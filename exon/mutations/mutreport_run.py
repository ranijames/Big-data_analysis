import os
import argparse
from logging import basicConfig, DEBUG
import logging
import sys
import subprocess
from os.path import join

import pandas as pd
import numpy as np

from allpy.SimpleTask import PatientTask
from allpy.config import pipeline_config
from allpy.config.config_constants import GENOME, NOTIFY_EMAIL, CORES, BCBIO_BACKUP, MUTREPORT_CANCERGENES, \
    MUTREPORT_CGC
from allpy.exon.mutations.mpileup_commands import MPILEUP_CMD
from allpy.exon.mutations.gemini_extraction import GeminiExtractor
from allpy.exon.mutations.multivcf import MultiVCFLoader
from allpy.exon.mutations.utils import get_mut_batch_folder, write_df, chromosome_order, pileup_read_details, read_df

ALL_MUTS = "ALLMUTS"

__author__ = 'ARJ'


# logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)
logger = logging.getLogger()


def unify_called_by(df):
    if df.shape[0] == 1:
        return df
    called_by = set()
    for x in df.called_by.tolist():
        for caller in x.split(','):
            called_by.add(caller)

    df.called_by = ','.join(sorted(list(called_by)))
    df.caller_count = len(called_by)
    return df.drop_duplicates()


class MutationReporter(PatientTask):
    def __init__(self, patientid=None, sample_dir=None, prefix=""):
        super().__init__("MutationReporter")

        self._new_master_file = False
        self.conf_dict = self.load_configuration()

        if patientid is None and sample_dir is None:
            raise RuntimeError('Either sample_dir or patientid must be specified')
        if patientid is None:
            patientid = os.path.basename(sample_dir)
        if sample_dir is None:
            sample_dir = '/home/mpschr/mounts/all_data/vbp_all/processed/exon/{}'.format(patientid)
        self.debug = False
        self.patientid = patientid
        try:
            self.in_folder = get_mut_batch_folder(self.patientid, sample_dir, self.conf_dict[BCBIO_BACKUP])
        except FileNotFoundError:
            logging.warning('No bcbio folder found')
        self.out_folder = join(sample_dir, 'mutations')
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)
        self.filename_stem = join(self.out_folder, '{}{}_{}'.format(prefix, patientid, "{}"))
        self.filename = self.filename_stem.format("{}-{}.tsv")

        self._cancer_genes = None
        self._cgc = None
        self.load_cancer_genes()

    def get_mut_folder(self):
        return self.in_folder

    def get_master_report_file(self, time_points):
        time_points = sorted(time_points)
        tpstring = "-".join(time_points)
        return self.filename.format(tpstring, ALL_MUTS)

    def get_df(self, time_points, min_reads=1, min_alt_reads=1, min_frequency=0,
               pam_only=False, force_recreate=False, avoid_master_creation=True,
               debug=False, time_point_normal=None, max_normal_frequency=0.05,
               inject_file = None) -> pd.DataFrame:

        self.debug = debug

        self.metadata = []
        df = self.get_master_report_df(avoid_master_creation, force_recreate, time_points, time_point_normal)

        if df is None:
            raise RuntimeError("No mutation dataframe could be retrieved!")

        obsolete_cols = ['gene_gemimini', 'exon_gemini', 'aa_change_gemini', 'impact_so']
        df.drop([x for x in obsolete_cols if x in df.columns], axis=1, inplace=True)

        df = self._anotate_cancer_genes(df)

        rows = df.shape[0]
        orig_rows = rows

        for tp in time_points:
            reads = 'reads_{}'.format(tp)
            freq = 'alt_freq_{}'.format(tp)
            cond = (df[reads] < min_reads) & (df[freq] > 0)
            if sum(cond) > 0:
                self.add_metadata("Set unreliable freq (reads < {}) to 0.0 for {} {}-variants".format(min_reads, sum(cond), tp))
                df.loc[cond, freq] = 0.0

        pam_col = 'pam_genes_count'

        if pam_only:
            query = "{} > 0 and biotype == 'protein_coding'".format(pam_col)
            df, rows = self.df_filter(df, query, rows)

        # filter tumor samples - min_reads
        query = ' or '.join(["reads_{} >= {}".format(tp, min_reads) for tp in time_points])
        df, rows = self.df_filter(df, query, rows)

        # filter tumor samples - min_alt_reads
        query = ' or '.join(["alt_reads_{} >= {}".format(tp, min_alt_reads) for tp in time_points])
        df, rows = self.df_filter(df, query, rows)

        # filter tumor samples min alt freqency
        query = ' or '.join(["alt_freq_{} >= {}".format(tp, min_frequency) for tp in time_points])
        df, rows = self.df_filter(df, query, rows)

        # filter normal samples - min_reads
        query = "reads_{} >= {}".format(time_point_normal, min_reads)
        df, rows = self.df_filter(df, query, rows)

        # filter normal samples max alt freqency
        query = 'alt_freq_{} < {}'.format(time_point_normal, max_normal_frequency)
        df, rows = self.df_filter(df, query, rows)

        query = 'chrom != "MT"'
        df, rows = self.df_filter(df, query, rows)

        inject_file_name = os.path.join(self.out_folder, "report_filter.py")
        if os.path.exists(inject_file_name):
            inject_file = inject_file_name
        if not inject_file is None:

            inject_file = os.path.abspath(inject_file)
            scope = {}
            exec(open(inject_file).read(), scope)

            df = scope["inject_func"](df, self.patientid)
            if df.shape[0] < rows:
                self.add_metadata("{} variants discarded with: {}".format(df.shape[0] - rows, "custom patient filter: {}".format(os.path.basename(inject_file))))
                rows = df.shape[0]
            for q in scope["patient_queries"]:
                df, rows = self.df_filter(df, q, rows)

            df = scope["inject_posthook"](df, self.patientid)


        self.add_metadata('{} variants returned ({:.2f}%)'.format(rows, rows / orig_rows * 100))

        for tp in time_points:
            tp_variants = df.query("alt_freq_{} > 0".format(tp)).shape[0]
            self.add_metadata('    {} variants returned for {}_{}'.format(tp_variants, self.patientid, tp))
        for tp in time_points:
            back_tracked = df.query("alt_freq_{0} > 0 and alt_freq_{0} < 0.05".format(tp)).shape[0]
            self.add_metadata(
                '    {} variants backtracked from other sample in {}_{}'.format(back_tracked, self.patientid, tp)
            )


        return df

    def load(self, inject_file):
        inject_file = os.path.abspath(inject_file)
        exec(open(inject_file).read(), globals(), locals())


    def add_metadata(self, metastring, do_log=True):
        self.metadata.append(metastring)
        if do_log:
            logging.info(metastring)

    def df_filter(self, df, query, rows):
        df_new = df.query(query)
        rows_new = df_new.shape[0]
        if rows_new < rows:
            self.add_metadata("{} variants discarded with: {}".format(rows_new - rows, query))
        return df_new, rows_new

    def get_standard_df(self, timepoints, time_point_normal, coding_only, force_recreate_master=False):

        min_reads = 10
        min_alt_reads = 3
        min_freq = 0.1

        ext = get_filename_extension(coding_only, min_alt_reads, min_freq, min_reads)
        timepoints = sorted(timepoints)
        tpstring = "-".join(timepoints)
        standard_filename = os.path.join(self.out_folder,
                                         self.filename.format(tpstring, ext))

        if os.path.exists(standard_filename) and not force_recreate_master:

            metadata,df = read_df(standard_filename, comment="#", index_col=None,dtype={"chrom": object})
            print('METADATA')
            for m in metadata:
                print(m.replace('\n',''))
            return df
            #return pd.read_table(standard_filename, index_col=None, comment='#', dtype={"chrom": object})



        df = self.get_df(timepoints, min_reads, min_alt_reads, min_freq,
                         avoid_master_creation=False, pam_only=coding_only, time_point_normal=time_point_normal,
                         force_recreate=force_recreate_master)


        if not os.path.exists(standard_filename) or self._new_master_file:
            logging.info("Written file: {}".format(standard_filename))
            write_df(df, standard_filename, self.metadata, False)
        return df

    def get_master_report_df(self, avoid_master_creation, force_recreate, time_points,
                             time_point_normal) -> pd.DataFrame:
        report_file = self.get_master_report_file(time_points)
        df = None
        if os.path.exists(report_file) and not force_recreate:
            df = pd.read_table(report_file, sep='\t', index_col=None, comment='#', dtype={"chrom": object})
            self.add_metadata('{} variants loaded from from {}'.format(df.shape[0], report_file))
        elif not avoid_master_creation:
            # create master report
            df = self.create_master_report(time_points, report_file, time_point_normal)
            self.add_metadata("Master report with {} variants created at: {}".format(df.shape[0], report_file))
            self._new_master_file = True
        else:
            logging.warning('Master report not created yet (and not allowed to create)!')

        if 'start' in df.columns and 'end' in df.columns:
            df.insert(1, 'pos', df['start'] + 1)
            df.drop(['start', 'end'], axis=1, inplace=True)
        return df

    def create_master_report(self, time_points, report_file, time_point_normal):
        master_metadata = []
        if self.debug:
            logger.debug('DEBUG MODE ON')

        limit = 20 if self.debug else 0

        # get gemini DFs
        sample_dfs = []
        for tp in time_points:
            sampleid = self.patientid + "_" + tp
            geminidf = GeminiExtractor(self.in_folder, sampleid).get_df(min_reads=1, min_alt_reads=1, coding_only=False)
            geminidf.chrom = geminidf.chrom.str.replace('chr', '')
            if limit > 0:
                geminidf = geminidf.sample(n=round(limit / len(time_points)))
            metadata_instance = "Gemini shape for {}: {}".format(sampleid, geminidf.shape)
            logger.info(metadata_instance)
            master_metadata.append(metadata_instance)
            if self.debug:
                logger.debug(geminidf.head())

            logger.info("Reading VCFs for sample")
            vcf_df = MultiVCFLoader(self.in_folder, sampleid).get_df(self.debug)
            metadata_instance = "VCF df shape for {}: {}".format(sampleid, vcf_df.shape)
            logger.info(metadata_instance)
            master_metadata.append(metadata_instance)

            geminidf.rename(columns={'aa_change': 'aa_change_gemini', 'exon': 'exon_gemini', 'gene': 'gene_gemini'},
                            inplace=True)
            sdf = vcf_df.reset_index().rename(columns=lambda x: x.lower()).merge(geminidf)

            metadata_instance = "Sample-DF shape for {}: {}".format(sampleid, sdf.shape)
            logger.info(metadata_instance)
            master_metadata.append(metadata_instance)

            diffgenes = (
                sdf[['chrom', 'pos', 'alt', 'ref', 'pam_genes_count', 'gene_gemini', 'genes', 'ct_vcf']]
                    .query('pam_genes_count > 0')
                    .replace("None", np.NaN)
                    .fillna("")
                    .query('gene_gemini != genes and genes != ""')
            )
            pd.set_option('display.width', 1000)

            if diffgenes.shape[0] > 0:
                metadata_instance = "{} Variants mapped to different genes in  {}".format(diffgenes.shape[0], sampleid)
                logger.info(metadata_instance)
                master_metadata.append(metadata_instance)
            sample_dfs.append(sdf)

        # create bed files for .bam files
        positions = (
            pd.concat(sample_dfs).drop_duplicates()
                .assign(chrom=lambda df: df.chrom.astype("category", categories=chromosome_order(df.chrom)))
                .sort_values(["chrom", "pos"])
        )

        # remove duplicate entries and aggregate called_by columns between samples
        positions = positions.groupby(['chrom', 'pos', 'ref', 'alt']).apply(unify_called_by)

        assert positions[positions.duplicated(['chrom', 'pos', 'ref', 'alt'], keep=False)].shape[0] == 0

        metadata_instance = "Merged DF samples for {}:  {}".format(self.patientid, positions.shape)
        logger.info(metadata_instance)
        master_metadata.append(metadata_instance)

        bed = self.filename_stem.format("variant-sites.bed")
        vcf = self.filename_stem.format("variant-sites.vcf")
        if (os.path.exists(bed) and
                os.path.exists(vcf) and
                    pd.read_table(bed, header=None).shape[0] == positions.shape[0]):
            logger.info("Using existing VCF file")
        else:
            self.run_mpileup(positions, bed, vcf, time_point_normal, time_points)

        samples = time_points + [time_point_normal]
        header = 'chrom,pos,id,ref,alt,QUAL,FILTER,INFO,FORMAT,{}'.format(','.join(samples)).split(',')
        clean_pileup_df = (
            pd.read_table(vcf, comment='#', header=-1, dtype={0: object})
                .rename(columns=lambda x: header[x])
                .drop(['id', 'QUAL', 'FILTER', 'INFO'], axis=1)
                .apply(lambda r: pileup_read_details(r, samples), axis=1)
                .drop(['FORMAT'] + samples, axis=1)
        )

        print("clean_pileup_df\n", clean_pileup_df.head())
        print("positions\n", clean_pileup_df.head())

        result_all_df = (
            clean_pileup_df.merge(positions, how='right')
                .drop_duplicates(['chrom', 'pos', 'ref', 'alt'])
                .assign(chrom=lambda df: df.chrom.astype("category", categories=chromosome_order(df.chrom)))
                .sort_values(['chrom', 'pos', 'ref', 'alt'])
        )
        print("result_all_df\n", result_all_df.head())

        null_variants = result_all_df[pd.isnull(result_all_df['reads_{}'.format(time_points[0])])].shape[0]

        metadata_instance = "{} variants have NaN values - discarded by mpileup in {}".format(null_variants,
                                                                                              self.patientid)
        logger.info(metadata_instance)
        master_metadata.append(metadata_instance)

        intcols = ['pos'] + [x for x in result_all_df.columns if 'reads_' in x]
        for intcol in intcols:
            result_all_df[intcol] = result_all_df[intcol].fillna(0).astype(int)

        metadata_instance = "final shape for {}: {}".format('merged', result_all_df.shape)
        logger.info(metadata_instance)
        master_metadata.append(metadata_instance)

        result_all_df['patient'] = self.patientid

        write_df(result_all_df, report_file, metadata=master_metadata, index=False)

        return result_all_df

    def run_mpileup(self, positions, bed, vcf, time_point_normal, time_points):
        bed_df = positions[['chrom', 'pos']]

        bed_df.to_csv(bed, sep="\t", index=False, header=False)
        samples = time_points + [time_point_normal]
        bams = {tp: os.path.join(os.path.dirname(self.in_folder),
                                 self.patientid + "_" + tp,
                                 self.patientid + "_" + tp + "-ready.bam")
                for tp in samples}

        # pileup all variants
        mpileup_cmd_format = MPILEUP_CMD.format(**{'GENOME': self.conf_dict[GENOME], 'BED': bed,
                                                   'TUMOR_SAMPLES': ' '.join(bams[tp] for tp in time_points),
                                                   'NORMAL_SAMPLE': bams[time_point_normal], 'OUTPUT_VCF': vcf})
        logging.debug(mpileup_cmd_format)
        p2 = subprocess.Popen(mpileup_cmd_format,
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=True)
        out2, err2 = p2.communicate()

    def get_metadata(self):
        return self.metadata

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            CORES,
                                            GENOME,
                                            BCBIO_BACKUP,
                                            MUTREPORT_CANCERGENES,
                                            MUTREPORT_CGC])
        return conf_dict

    def _anotate_cancer_genes(self, df):


        idx = df.columns.tolist().index('pam_genes_count') + 1

        df.insert(idx,
                           'known_cancer_gene',
                           df.genes.apply(lambda genes: any([x in self._cancer_genes.index for x in str(genes).split(',')]))
                          )
        df.insert(idx+1,
                          'cancer_gene_source',
                           df.genes.apply(lambda genes: max(['']+[self._cancer_genes.source.loc[x] for x in str(genes).split(',') if x in self._cancer_genes.index]))
                          )
        df.insert(idx+2,
                          'cancer_gene_source_count',
                          df.cancer_gene_source.apply(len))
        df.insert(idx+3,
                           'cgc_mol_genetics',
                            df.genes.apply(lambda genes: max([''] + [self._cgc.molecular_genetics.loc[x] for x in str(genes).split(',') if x in self._cgc.index], key=lambda x: len(str(x))))
                           )
        df.insert(idx+4,
                           'cgc_mut_type',
                            df.genes.apply(lambda genes: ",".join([str(self._cgc.mutation_types.loc[x]) for x in str(genes).split(',') if x in self._cgc.index]))
                           )

        return df

    def load_cancer_genes(self):

        cancer_gene_dfs = []
        gene_lists = self.conf_dict[MUTREPORT_CANCERGENES]
        self._cgc = pd.read_table(self.conf_dict[MUTREPORT_CGC]).set_index('gene_symbol')

        blacklist = []
        for source, filename in gene_lists.items():
            cgdf = (pd.read_table(filename, header=None, names=[source])
                    .assign(source=source.upper()[0])
                    .set_index(source)
                    [['source']]
                    )
            if source == "blacklist":
                blacklist = cgdf.index.tolist()
            else:
                cancer_gene_dfs.append(cgdf)

        cancer_gene_dfs.append(pd.DataFrame({'source': 'C'}, index=self._cgc.index))

        self._cancer_genes = (
            pd.concat(cancer_gene_dfs)
                .reset_index()
                .rename(columns={'index': 'symbol'})
                .drop_duplicates()
                .pivot_table(index=['symbol'], aggfunc=lambda l: ''.join(sorted(l)))
                .assign(list_count=lambda df: df.source.apply(len))
        ).query('index != "{}"'.format(blacklist))
        self._cancer_genes = self._cancer_genes[~self._cancer_genes.index.isin(blacklist)]


def annotate_protein_coding(df, sequence_ontology_colname, target_colname):
    if sequence_ontology_colname in df:
        df[target_colname] = df[sequence_ontology_colname].map(protein_coding_map)
        unmapped_so = df[sequence_ontology_colname][df[target_colname] == '?'].unique().tolist()
        if len(unmapped_so) > 0:
            logger.warning('unmapped SOs: {}'.format(unmapped_so))
    else:
        logging.info('Column "{}" not present'.format(sequence_ontology_colname))


def protein_coding_map(so):
    if so in ['5_prime_UTR_premature_start_codon_gain_variant',
              'downstream_gene_variant',
              'exon_variant',
              'intergenic_variant',
              'intron_variant',
              'upstream_gene_variant',
              'downstream_gene_variant', '3_prime_UTR_variant', 'synonymous_variant', '5_prime_UTR_variant']:
        return 'not_coding'
    elif so in ['frameshift_variant',
                'missense_variant',
                'splice_region_variant',
                'stop_gained',
                'stop_retained_variant',
                'inframe_deletion', 'disruptive_inframe_deletion', 'inframe_insertion',
                'start_lost', 'stop_lost', 'initiator_codon_variant', 'disruptive_inframe_insertion']:
        return 'is_coding'
    elif so == 'None':
        return so
    else:
        return '?'


def get_filename_extension(coding_only, min_alt_reads, min_frequency, min_reads):
    filename_extensions = []
    if coding_only:
        filename_extensions.append('coding')
    if min_reads > 1:
        filename_extensions.append("r{}".format(min_reads))
    if min_alt_reads > 1:
        filename_extensions.append("a{}".format(min_alt_reads))
    if min_frequency > 0:
        filename_extensions.append("f{}".format(min_frequency))
    ext = ".".join(filename_extensions)
    if ext == "":
        ext = ALL_MUTS
    return ext


class MutationReportRun:
    def __init__(self):
        pass

    def run(self, patientid, patient_dir, min_reads, min_alt_reads, min_frequency, time_points=None,
            time_point_normal='CR',
            coding_only=False, debug=False, standard_report=False, force_recreate_master=False):

        if time_points == None or time_points == []:
            time_points = ["ID", "REL"]
        logging.debug(time_points)

        if patientid is None:
            patientid = os.path.basename(patient_dir)

        ext = get_filename_extension(coding_only, min_alt_reads, min_frequency, min_reads)

        logging.info("Running mutreport for {} ({} and {})".format(patientid, ",".join(time_points), time_point_normal))
        mutreporter = MutationReporter(patientid, patient_dir)

        if standard_report:
            df_mutations = mutreporter.get_standard_df(time_points, time_point_normal, coding_only,
                                                       force_recreate_master=force_recreate_master)
        else:
            df_mutations = mutreporter.get_df(time_points, min_reads, min_alt_reads, min_frequency,
                                              pam_only=coding_only, debug=debug, time_point_normal=time_point_normal,
                                              avoid_master_creation=False, force_recreate=force_recreate_master)
        if ext != ALL_MUTS and not standard_report:
            time_points = sorted(time_points)
            tpstring = "-".join(time_points)
            out_file_name = join(mutreporter.out_folder, mutreporter.filename.format(tpstring, ext))
            write_df(df_mutations, out_file_name, mutreporter.metadata)
            logging.info("Written file: {}".format(out_file_name))


def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-D', dest='patient_dir', required=True,
                        help='Directory with sample name. The master file should be '
                             'found at the following path: sample/mutations/sample_ID-ALLMUTS.tsv')
    parser.add_argument('-r', dest='min_reads', default=1, help="Minimum depths for variants to be considered",
                        type=int)
    parser.add_argument('-a', dest='min_alt_reads', default=1, help="Minimum alt depths for variants to be considered",
                        type=int)
    parser.add_argument('-f', dest='min_frequency', default=0, help="Minimum variant frequency (VAF) to be considered",
                        type=float)
    parser.add_argument('-c', dest='coding_only', action='store_true', default=False,
                        help="Protein coding mutations only")
    parser.add_argument('-t', dest='time_points', action='append')
    parser.add_argument('-T', dest='time_point_normal', default='CR')
    parser.add_argument('--get-standard', dest='standard', action='store_true', default=True)
    parser.add_argument('--debug', dest='debug', action='store_true', default=False)
    parser.add_argument('-M', dest='force_recreate_master', action='store_true', default=False,
                        help="Force the re-creation of the Master mutation table")
    options = parser.parse_args()
    task = MutationReportRun()
    task.run(options.patientid, os.path.abspath(options.patient_dir), options.min_reads, options.min_alt_reads,
             options.min_frequency,
             options.time_points, options.time_point_normal, options.coding_only, options.debug, options.standard,
             options.force_recreate_master)


if __name__ == "__main__":  # detects if called from system
    cmdline()
