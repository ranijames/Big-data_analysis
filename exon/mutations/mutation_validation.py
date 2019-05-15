import os
import pandas as pd
import argparse
from logging import basicConfig, DEBUG
import logging
import sys
import subprocess
from allpy.SimpleTask import PatientTask
from allpy.config import pipeline_config
from allpy.config.config_constants import GENOME, NOTIFY_EMAIL, CORES, BCBIO_BACKUP
from allpy.exon.mutations.mpileup_commands import MPILEUP_CMD
from allpy.exon.mutations.mutreport_run import MutationReporter
from allpy.exon.mutations.utils import get_mut_batch_folder, write_df, chromosome_order, pileup_read_details
from os.path import join
import numpy as np
from allpy.externalcommands import ExternalCommand

NO_PANEL_READS = "NO_PANEL_READS"


ALL_MUTS = "ALLMUTS"

__author__ = 'ARJ'


# logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)
logger = logging.getLogger()


class MutValidator(PatientTask):
    def __init__(self, patientid=None, patient_dir_exon=None, patient_dir_panel=None, prefix="", time_points=['ID','REL'], time_point_normal='CR'):
        super().__init__("MutationValidator")

        self.conf_dict = self.load_configuration()

        if patientid is None and patient_dir_panel is None:
            raise RuntimeError('Either sample_dir or patientid must be specified')
        if patientid is None:
            patientid = os.path.basename(patient_dir_panel)

        self.debug = False
        self.patientid = patientid

        self.out_folder = join(patient_dir_panel, 'panel', 'mutations')
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)
        self.filename_stem = join(self.out_folder, '{}{}_{}'.format(prefix, patientid, "{}"))
        self.filename = self.filename_stem.format("{}-{}.tsv")


        reporter = MutationReporter(patientid, patient_dir_exon)
        bed = reporter.filename_stem.format("variant-sites.bed")
        #exon_vcf = self.filename_stem.format("variant-sites.vcf")
        panel_vcf = self.filename_stem.format("variant-sites-panel.vcf")
        panel_tsv = panel_vcf.replace('vcf', 'tsv')
        if os.path.exists(panel_tsv):
            logger.warn("OUTPUT file exists - exiting: {}".format(panel_tsv))
            return

        exon_mutations = reporter.get_standard_df(time_points, time_point_normal=time_point_normal, coding_only=True)

        time_points = [''.join([i for i in x if not i.isdigit()]) for x in time_points]

        for t in [time_point_normal] + time_points:
            if "alt_freq_{}".format(t) not in exon_mutations.columns:
                renamed_sample = False
                for x in [1, 2, 3, 4]:
                    if "alt_freq_{}{}".format(t, x) in exon_mutations.columns:
                        alt_t = t + str(x)
                        print("renaming {} to {}".format(alt_t, t))
                        rename_dict = {c : c.replace(alt_t, t) for c in exon_mutations.columns if alt_t in c}
                        print(rename_dict)
                        exon_mutations.rename(columns=rename_dict,inplace=True)
                        renamed_sample=True
                if not renamed_sample:
                    print("{} not in EXON muts".format(t))

        bed_stem = join(patient_dir_panel, "panel/{0}_{1}/{0}_{1}-ready.bam").format(patientid,"{0}")
        panel_bams = {tp:bed_stem.format(tp) for tp in time_points + [time_point_normal]}

        self.run_mpileup(bed, panel_vcf, time_point_normal, time_points, panel_bams)

        mpileup_df = self.read_mpileup_vcf(exon_mutations, panel_vcf, time_point_normal, time_points)
        df = mpileup_df.merge(exon_mutations, how='right', on=['chrom', 'pos', 'ref', 'alt'],
                           suffixes=['', '_WES'])

        df = df.apply(self.correct_empty_data, axis=1, mpileup_df=mpileup_df)
        df['panel_flag_ID'] = df.apply(self.panel_flags,axis=1, timepoint='ID')
        df['panel_flag_REL'] = df.apply(self.panel_flags, axis=1, timepoint='REL')
        df['panel_flag_COV'] = df.apply(self.panel_coverage_flags,axis=1)

        #df.query('panel_flag_ID == {}'.format(NO_PANEL_READS)).apply(lambda r: correct_no_panel, axis=1)

        metadata = ['']

        id_df = df.query("alt_freq_ID_WES > 0.05")
        rel_df = df.query("alt_freq_REL_WES> 0.05")

        metadata.append("{}/{} variants confirmed for {} ({}%) ".format(id_df.panel_flag_ID.value_counts().ix['OK'], id_df.shape[0], patientid + "_ID",id_df.panel_flag_ID.value_counts().ix['OK']/id_df.shape[0]))
        if "OK" in rel_df.panel_flag_REL.values:
            metadata.append("{}/{} variants confirmed for {} ({}%)".format(rel_df.panel_flag_REL.value_counts().ix['OK'], rel_df.shape[0], patientid + "_REL", rel_df.panel_flag_REL.value_counts().ix['OK']/rel_df.shape[0]))
        else:
            metadata.append("{}/{} variants confirmed for {} ({}%)".format(0, rel_df.shape[0], patientid + "_REL", 0))

        if NO_PANEL_READS in df.panel_flag_ID.values:
            metadata.append("no panel reads: {}".format(df.panel_flag_ID.value_counts().ix[NO_PANEL_READS]))

        for m in metadata:
            logger.info(m)

        write_df(df , panel_tsv, index=False, metadata=metadata)
        #df.to_csv(panel_vcf.replace('vcf', 'tsv'), sep='\t', index=None)

    def read_mpileup_vcf(self, exon_mutations, panel_vcf, time_point_normal, time_points):
        samples = time_points + [time_point_normal]
        header = 'chrom,pos,id,ref,alt,QUAL,FILTER,INFO,FORMAT,{}'.format(','.join(samples)).split(',')
        clean_pileup_df = (
            pd.read_table(panel_vcf, comment='#', header=-1, dtype={0: object})
                .rename(columns=lambda x: header[x])
                .drop(['id', 'QUAL', 'FILTER', 'INFO'], axis=1)
                .apply(lambda r: pileup_read_details(r, samples), axis=1)
                .drop(['FORMAT'] + samples, axis=1)
        )

        return clean_pileup_df

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            CORES,
                                            GENOME,
                                            BCBIO_BACKUP])
        return conf_dict

    def run_mpileup(self, bed, vcf, time_point_normal, time_points, bams):

        # pileup all variants
        mpileup_cmd_format = MPILEUP_CMD.format(**{'GENOME': self.conf_dict[GENOME], 'BED': bed,
                                                   'TUMOR_SAMPLES': ' '.join(bams[tp] for tp in time_points),
                                                   'NORMAL_SAMPLE': bams[time_point_normal], 'OUTPUT_VCF': vcf})

        cmd = ExternalCommand(mpileup_cmd_format,output_file=vcf)
        ret_code = cmd.run()


    @staticmethod
    def correct_empty_data(r,mpileup_df):
        if all(pd.isnull(r[['reads_ID', 'reads_REL', 'reads_CR']])):
            for tp in ['ID', 'REL', 'CR']:
                real_reads = mpileup_df.query('pos == {} and chrom == "{}"'.format(r.pos,r.chrom))['reads_{}'.format(tp)].max()
                if real_reads > 0:
                    r['reads_{}'.format(tp)] = real_reads
                    r['alt_reads_{}'.format(tp)] = 0
                    r['alt_freq_{}'.format(tp)] = 0.0

        return r

    @staticmethod
    def panel_flags(r,timepoint):

        flag = ["OK"]

        if all(pd.isnull(r[['reads_ID', 'reads_REL', 'reads_CR']])):
            flag = [NO_PANEL_READS]
        elif any(pd.isnull(r[['alt_freq_ID', 'alt_freq_REL', 'alt_freq_CR']])):
                flag = ["SOME_MISSING_PANEL_READS" ]
        elif r.alt_freq_CR > 0.05 and r.alt_freq_CR_WES < 0.05:
            flag =['PANEL_GERMLINE']

        else:
            panel = 'alt_freq_{}'.format(timepoint)
            wes = panel + "_WES"
            altflag = []
            if r[panel] == 0.0 and r[wes] > 0.05:
                altflag.append('{}_DISAPPEARED'.format(timepoint))
            elif r[panel] < 0.05 and r[wes] > 0.05:
                 altflag.append('{}_BELOW_THRESHOLD'.format(timepoint))

            if len(altflag) > 0:
                flag = altflag
        return ';'.join(sorted(flag))

    @staticmethod
    def panel_coverage_flags(r):
        cov_insufficient = []
        cov_low = []
        for tp in ['ID','CR','REL']:
            if r['reads_{}'.format(tp)] < r['reads_{}_WES'.format(tp)] or r['reads_{}'.format(tp)] < 50:
                cov_insufficient.append(tp)
            elif r['reads_{}'.format(tp)] < 100:
                cov_low.append(tp)
        i = None
        if len(cov_insufficient) > 0:
            i = 'INSUFFICIENT:{}'.format(','.join(sorted(cov_insufficient)))
        l = None
        if len(cov_low) > 0:
            l = 'LOW:{}'.format(','.join(sorted(cov_low)))

        return ';'.join(sorted([x for x in [i,l] if x is not None]))





class MutationReportRun:
    def __init__(self):
        pass

    def run(self, patientid, patient_dir_exon, patient_dir_panel, time_points=None,
            time_point_normal='CR'):

        if time_points == None or time_points == []:
            time_points = ["ID", "REL"]
        logging.debug(time_points)

        if patientid is None:
            patientid = os.path.basename(patient_dir_exon)

        logging.info("Running mutreport for {} ({} and {})".format(patientid, ",".join(time_points), time_point_normal))
        mutreporter = MutValidator(patientid, patient_dir_exon, patient_dir_panel,
                                   time_points=time_points, time_point_normal=time_point_normal)








def cmdline():

    print(sys.argv)
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-D', dest='patient_dir_exon', required=True,
                        help='Directory with sample name. The master file should be '
                             'found at the following path: sample/mutations/sample_ID-ALLMUTS.tsv')
    parser.add_argument('-d', dest='patient_dir_panel', required=True,
                    help='Directory with sample name. The master file should be '
                         'found at the following path: sample/mutations/sample_ID-ALLMUTS.tsv')

    parser.add_argument('-t', dest='time_points', action='append')
    parser.add_argument('-T', dest='time_point_normal', default='CR')
    parser.add_argument('--debug', dest='debug', action='store_true', default=False)
    options = parser.parse_args()
    task = MutationReportRun()
    task.run(options.patientid, os.path.abspath(options.patient_dir_exon), os.path.abspath(options.patient_dir_panel),
             options.time_points, options.time_point_normal)


if __name__ == "__main__":  # detects if called from system
    cmdline()
