import os
import argparse
from logging import basicConfig, DEBUG
import logging
import sys
from matplotlib.pyplot import savefig

import pandas as pd
import math
from pybedtools import BedTool
import numpy as np
import time

from cnvlib.commands import do_call as cnvkit_call, do_scatter as cnvkit_scatterplot, \
     do_breaks as cnvkit_breaks, do_gainloss as cnvkit_gainloss
from cnvlib.vary import VariantArray
from cnvlib.cnary import CopyNumArray
import pylab
from skgenome import tabio
from allpy.exon.mutations.patient_snps import PatientSNPer, EXAC, BED_FILES

from allpy.exon.mutations.utils import Timer, write_df, chromosome_order
from os.path import join
from allpy.notify.notify import send_email
from allpy.pipelineutils.sample_info import get_sample_info

ALLCNAS = 'ALLCNAS'

REL = "REL"

ID = "ID"

__author__ = 'ARJ'

GENES_BED_FILE = '/home/mpschr/bin/bcbionextgen/data/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.gtf'


popvar_sites = "/home/mpschr/Data/snp_sites/exac_variant_sites.snp.tsv"
popvar_sites_df = pd.read_table(popvar_sites, header=None, names=['#CHROM','POS', 'REF', 'ALT', 'MAF', 'AN', 'AC', 'AC_Hom', 'FILTER'], dtype={'#CHROM': str})
popvar_sites_df.head()

#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)
logger = logging.getLogger()

corr_plot_centering = {
    'AE05_ID': 0.08,
    'AE05_REL': 0.05,
    'AE29_ID': 0.4,
    'AE29_REL': 0.3,
    'AL03_ID': 0.18,
    'AL03_REL': 0.07,
    'AL10_ID': 0.25,
    'AL10_REL': 0.15,
    'AL27_ID': 0.1,
    'AL27_REL': 0.1,
    'AL28_ID': 0.1,
    'AL28_REL': 0.06,
    'AL32_ID': 0.3,
    'AL32_REL': 0.2,
    'PL07_ID': 0.1,
    'PL07_REL': 0.1
}

corr_purity = {
    'AE04_REL': 68,
    'AL03_REL': 55,
    'AL10_ID': 63,
    'AL10_REL': 76,
    'AL28_ID': 89,
    'AL32_REL': 79,
    'AL04_ID': 43,
    'AL04_REL': 50,
}

def get_correction(samplename, sample_info):
    if samplename in sample_info.columns and 'cna_log2_center_correction' in sample_info.index:
        return sample_info.ix['cna_log2_center_correction', samplename]
    if samplename in corr_plot_centering:
        return corr_plot_centering[samplename]

    return 0


def get_segment_file(patient_dir):
    return join(patient_dir, 'copywriteR', 'CNAprofiles', 'log2_CNA.segmented.tsv')

def get_log2_ratio_file(patient_dir, not_bcbio, sample=None):
    if not_bcbio:
        return join(patient_dir, 'copywriteR', 'CNAprofiles', 'log2_CNA.igv')
    else:
        return join(patient_dir, 'bcbio',sample,sample+"-cnvkit.cnr")

def get_bin_kb(patient_dir):
        return int([x for x in os.listdir('{}/copywriteR/CNAprofiles'.format(patient_dir)) if x.endswith('kb.info')][0].split('kb')[0])


def add_genes(row: pd.Series, cnvkit_cnr: CopyNumArray):
    genes = set()
    for ix, probe in cnvkit_cnr.in_range(row.chromosome, row.start, row.end, mode='trim').data.iterrows():
        if pd.isnull(probe.gene):
            continue
        for gene in probe.gene.split(','):
            if gene != "":
                genes.add(gene)

    row.gene = ",".join(sorted(genes))
    return row


class CNAReporter:

    def __init__(self, patientid, patient_dir, genes_bed_file=GENES_BED_FILE, not_bcbio=False):
        self.debug = False
        self.genes_bed = BedTool(genes_bed_file)
        self.patientid = patientid
        self.patient_dir = patient_dir
        self.out_folder = join(patient_dir, 'cna')
        self.not_bcbio = not_bcbio
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)

    def get_master_report_file(self, time_point):
        return os.path.join(self.out_folder, "{}_{}.{}.tsv".format(self.patientid, time_point, ALLCNAS))


    def get_df(self, time_point, normal_time_point, popvar_vcf_file,
               different_only=False, coding_only=False, from_list=None) -> pd.DataFrame:

        self.metadata = []
        self.popvar_vcf_file = popvar_vcf_file
        df = self.get_master_report_df(time_point, normal_time_point)

        if df is None:
            raise RuntimeError("No CNA dataframe could be retrieved!")

        return df

        rows = df.shape[0]
        orig_rows = rows

        if different_only:
            query = "CNA < 2 or CNA > 3"
            df, rows = self.df_filter(df, query, rows)

        if coding_only:
            query = "biotype == 'protein_coding'"
            df, rows = self.df_filter(df, query, rows)

        if from_list is not None and len(from_list) > 0:
            query = "gene == {}".format(from_list)
            df, rows = self.df_filter(df, query, rows, 'in required list')


        segdf = self.get_segments_df(time_point)
        selected_segs = df.sort('segments').segments.unique().tolist()
        selected_kbs = segdf[segdf.segment.isin(selected_segs)][['seg.length.kb']]

        self.add_metadata('{} ({:.2f}%) genes from {} segments along {} kbs returned'.format(
            rows, rows/orig_rows*100, len(df.segments.unique()), selected_kbs.sum()[0] ))

        selected_kbs_per_CN = segdf[segdf.segment.isin(selected_segs)].groupby(['chrom','copynumber']).agg(
            {'seg.length.kb': np.sum}
        ).sort()

        for ids, seglen in selected_kbs_per_CN.iterrows():
            self.add_metadata("chrom {}, CN {}: {} kbs".format(ids[0], ids[1], seglen[0]))
        return df

    def add_metadata(self, metastring, do_log=True):
        self.metadata.append(metastring)
        if do_log:
            logging.info(metastring)

    def df_filter(self, df, query, rows, alt_text=None):
        df_new = df.query(query)
        rows_new = df_new.shape[0]
        if rows_new < rows:
            text = query if alt_text is None else alt_text
            self.add_metadata("{} genes discarded with: {}".format(rows_new - rows, text))
        return df_new, rows_new

    def get_master_report_df(self, time_point, normal_time_point) -> pd.DataFrame:
        report_file = self.get_master_report_file(time_point)
        df = None
        if os.path.exists(report_file):
            df = pd.read_table(report_file, sep='\t', comment='#').reset_index()
        else:
            # create master report
            df = self.create_master_report(time_point, normal_time_point, report_file)

        self.add_metadata('{} genes loaded from from {}'.format(df.shape[0], report_file))
        return df


    def get_metadata(self):
        return self.metadata

    def create_master_report(self, time_point, normal_time_point, report_file):
        master_metadata = []

        sample = self.patientid + "_" + time_point
        normal_sample = self.patientid + "_" + normal_time_point

        self.annodf = get_sample_info(self.patient_dir)
        self.annodf['gender'] = self.annodf.gender.replace('W', 'female').replace('M', 'male')

        if sample in corr_purity:
            logger.info("correcting purity from {} to {}".format(annodf.loc[sample].purity, corr_purity[sample]))
            self.annodf.set_value(sample, 'purity', corr_purity[sample])

        cnr_filename = join(self.out_folder, "{}.cnr".format(sample))

        logger.debug("Getting log2 ratio df for sample {}".format(sample))
        cnvkit_cnr = self.get_log2_ratio_df(sample)
        cnvkit_cns = self.get_segments_df(sample, cnvkit_cnr)
        cnvkit_vaf = self.get_cnvkit_vaf(sample, normal_sample)

        logger.debug('pipeline provance set to (not_bcbio): {}'.format(self.not_bcbio))
        if self.not_bcbio:
            try:
                gender = self.annodf.loc[sample.replace('CR', 'REL').replace('REL2','REL')].gender
            except KeyError:
                gender = None

            try:
                purity = self.annodf.loc[sample].purity
                logger.info('purity: {}'.format(purity))
            except KeyError:
                purity = 100
            logging.info({"gender": gender, "purity": purity})

            calling_method = 'clonal' if purity > 90 else 'threshold'

            cnvkit_called = cnvkit_call(cnvkit_cns, variants=cnvkit_vaf,
                             is_sample_female= gender == 'female',
                             is_reference_male= gender == 'male',
                             purity=purity/100, method=calling_method)

            if time_point != normal_time_point:
                breaks = pd.DataFrame(cnvkit_breaks(cnvkit_cnr, cnvkit_called)).fillna("").replace('nan', '')
                breaks.columns = ['Gene', 'Chrom.', 'Location', 'Change',
                                 'ProbesLeft', 'ProbesRight']
                #print(breaks.query("Gene != ''"))

                gainloss = pd.DataFrame(cnvkit_gainloss(cnvkit_cnr, segments=cnvkit_called, male_reference=gender == 'male'))
                #gainloss = pd.DataFrame(cnvkit_gainloss(cnvkit_cnr,min_probes=1, male_reference=gender == 'male'))
                print(cnr_filename.replace('.cnr', '.gainloss'))
                gainloss.to_csv(cnr_filename.replace('.cnr', '.gainloss'), sep="\t", index=None)
        else:
            cnvkit_called = pd.read_table(get_log2_ratio_file(self.patient_dir, self.not_bcbio, sample).replace('.cnr','-call.cns'),
                                    dtype={'chromosome': 'str'}).loc[lambda df: df.chromosome.str.startswith('GL') == False]
            cnvkit_called = CopyNumArray(cnvkit_called)


        #cnvkit_cnr.write(cnr_filename)
        tabio.write(cnvkit_cnr,cnr_filename)
        #cnvkit_called.write(cnr_filename.replace(".cnr", ".called.cns"))
        tabio.write(cnvkit_called,cnr_filename.replace(".cnr", ".called.cns"))

        metadata_instance = "called"
        logger.info(metadata_instance)
        master_metadata.append(metadata_instance)

        do_plots = True
        if do_plots:
            pylab.rcParams['figure.figsize'] = (25, 8)
            cnvkit_scatterplot(
                cnarr=cnvkit_cnr, segments=cnvkit_called, variants=cnvkit_vaf, do_trend=True,
                title=sample
            )
            savefig(join(self.out_folder, '{}.karyotype.png'.format(sample)))
            pylab.clf()

        #write_df(cnvkit_called.reset_index(), report_file, metadata=master_metadata, index=False)

        return cnvkit_called.data

    def get_log2_ratio_df(self, sample) -> CopyNumArray:
        log2ratios_file = get_log2_ratio_file(self.patient_dir, self.not_bcbio, sample)

        if self.not_bcbio:


            ## log2 ratios
            log2ratios = pd.read_table(log2ratios_file, comment="#",dtype={'Chromosome':str})

            log2ratio_cols = [x for x in log2ratios.columns if '.vs.' in x]
            log2ratios = (
                log2ratios[log2ratios.columns[0:4].tolist() + log2ratio_cols]
                .rename(columns=lambda x: x.lower() if not "_" in x else x.split('.vs.')[0])
                .assign(chromosome = lambda df: df.chromosome.astype("category", categories=chromosome_order(df.chromosome)))
                .sort_values(['chromosome', 'start'])
            )

            log2ratios['log2'] = log2ratios[[sample]] + get_correction(sample, self.annodf)

            log2ratios_df_cnvkit = log2ratios[['chromosome', 'start', 'end', 'gene', 'log2']]
            return CopyNumArray(log2ratios_df_cnvkit)
        else:
            return CopyNumArray(pd.read_table(log2ratios_file, dtype={'chromosome': 'str'}).loc[lambda df: df.chromosome.str.startswith('GL') == False])



    def get_segments_df(self, sample, cnvkit_cnr: CopyNumArray) -> CopyNumArray:

        segment_file_copywriter = get_segment_file(self.patient_dir)
        if self.not_bcbio:
            segments_cpwr = pd.read_table(segment_file_copywriter)

            # R converts '-' to . in ID name.
            if '-' in self.patientid:
                id_in_R = self.patientid.replace('-', '.')
                segments_cpwr['ID'] = segments_cpwr.ID.apply(lambda v: v.replace(id_in_R, self.patientid))


            if 'CR' in sample:
                querystring = '{}.vs.none'.format(sample)
                segments_cpwr = segments_cpwr[segments_cpwr.ID.str.contains('_CR.vs.none')]
            else:
                querystring = '{}.vs'.format(sample)
                segments_cpwr = segments_cpwr[~segments_cpwr.ID.str.contains("none") & segments_cpwr.ID.str.contains(querystring)]
            logging.info(segments_cpwr.ID.drop_duplicates().tolist())

            kb = get_bin_kb(self.patient_dir)

            correction = get_correction(sample, self.annodf)
            segments_cnvkit = (
                segments_cpwr
                    .assign(ID = lambda df: df.ID.str.replace('\..*', ''))
                    .query("ID == '{}'".format(sample))
                    .assign(start = lambda df: (df['loc.start'] - (kb*1000/2) + 1).astype(int))
                    .assign(end = lambda df: (df['loc.end'] + (kb*1000/2)).astype(int))
                    .assign(log2 = lambda df: df['seg.mean'] + correction)
                    .rename(columns={'chrom': 'chromosome',
                            'num.mark': 'probes',
                           })
                    .assign(gene = '-')
            )[['chromosome', 'start', 'end', 'gene', 'log2', 'probes']]

            segments_cnvkit = segments_cnvkit.apply(lambda row: add_genes(row, cnvkit_cnr), axis=1)
            return CopyNumArray(segments_cnvkit)
        else:
            mat = pd.read_table(get_log2_ratio_file(self.patient_dir, self.not_bcbio, sample).replace('cnr','cns'),
                                dtype={'chromosome': 'str'})\
                    .loc[lambda df: df.chromosome.str.startswith('GL') == False]
            return CopyNumArray(mat)

    def get_cnvkit_vaf(self, sample, normal_sample) -> VariantArray:
        """
        :param sample:
        :param normal_sample:
        :rtype : VariantArray
        """
        variants = tabio.read(self.popvar_vcf_file, fmt='vcf', sample_id=sample,
                              normal_id=normal_sample,
                              min_depth=40,
                              skip_somatic=False)
        variants_filtered = VariantArray(
            variants.data
                .merge(popvar_sites_df.rename(columns={
                        "#CHROM": 'chromosome',
                        "POS": 'end',
                        "REF": 'ref',
                        "ALT": 'alt'
                    }
                ))
                .query("n_alt_freq > 0.4 and n_alt_freq < 0.6")
        )
        return variants_filtered

    def logtime(self, row, start_time):
        current_time = time.time()
        elapsed = current_time - start_time
        mins = math.floor(elapsed / 60)
        secs = math.ceil(elapsed % 60)
        print('At segment {} - {}m{}s'.format(row['segment'], mins, secs), end='\r')



def get_ext(different_only, coding_only):
    filename_extensions = []
    if different_only:
        filename_extensions.append('altered')

    if coding_only:
        filename_extensions.append('coding')

    ext = ".".join(filename_extensions)
    if ext == "":
        ext = ALLCNAS
    return ext


class CNAReportRun:

    def __init__(self):
        pass

    def run(self, patient_dir, time_points=None, time_point_normal="CR", not_bcbio=False):

        patientid = os.path.basename(patient_dir)

        snper = PatientSNPer(time_points, time_point_normal, patient_dir, BED_FILES[EXAC], EXAC)
        vcf_file = snper.get_popvar_vcf()

        if time_points is None or time_points == []:
            time_points = [ID, REL]
        logging.debug(time_points)

        reporter = CNAReporter(patientid, patient_dir, not_bcbio=not_bcbio)

        for tt in time_points:

            sample = patientid + "_" + tt

            df = reporter.get_df(tt, time_point_normal, popvar_vcf_file=vcf_file, different_only=True, coding_only=True)

            out_file_name = join(reporter.out_folder, '{}-cna.tsv'.format(sample))
            #write_df(df, out_file_name, reporter.metadata, index=False)
            logging.info("Written file: {}".format(out_file_name))



def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-D', dest='sample_dir', required=True, help='Where the output folder from bcbio is, the root '
                                                                       'for the sample results')
    parser.add_argument('-t', dest='time_points', action='append')
    parser.add_argument('-T', dest='time_point_normal', default='CR')
    parser.add_argument('--copywriter', dest='copywriter',action='store_true',default=False,
                        help="Specify if CopywriteR has been run instead of the cnvkit within the bcbio pipeline.")

    options = parser.parse_args()
    print(options)
    task = CNAReportRun()
    task.run(os.path.abspath(options.sample_dir), time_points=options.time_points, time_point_normal=options.time_point_normal, not_bcbio=options.copywriter)

if __name__ == "__main__":  # detects if called from system
    cmdline()
