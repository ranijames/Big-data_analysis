import os
import argparse
from logging import basicConfig, DEBUG
import logging
import sys
from allpy.SimpleTask import PatientTask
from allpy.config import pipeline_config
from allpy.config.config_constants import GENOME, NOTIFY_EMAIL, BCBIO_BACKUP, \
    POPVAR_EXACSNP_BED, POPVAR_COSMIC_BED,  MUTREPORT_CGC
from allpy.exon.mutations.mpileup_commands import MPILEUP_CMD, MPILEUP_DELETIONS_CMD
from allpy.exon.mutations.utils import get_mut_batch_folder, write_df, chromosome_order
from os.path import join
from allpy.externalcommands import ExternalCommand

COSMIC = 'cosmic-rec5'
EXAC = 'exac-snps'
BED_FILES = {
    COSMIC: POPVAR_COSMIC_BED,
    EXAC: POPVAR_EXACSNP_BED
}

__author__ = 'ARJ'


# logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)
logger = logging.getLogger()


class PatientSNPer(PatientTask):
    def __init__(self, time_points, time_point_normal, patient_dir, bed_file, bed_name):
        super().__init__("PatientSNPer")

        self.patient_dir = patient_dir
        self.derive_patient_name(patient_dir, None)
        self._bed_file = bed_file
        self._bed_name = bed_name

        if time_points is None:
            time_points = ['ID', 'REL']
        self.time_points = time_points

        if time_point_normal is None:
            time_point_normal = "CR"
        self.time_point_normal = time_point_normal

        self._new_master_file = False
        self.conf_dict = self.load_configuration()
        self.folder_and_file_names(self.patient_dir)

    def derive_patient_name(self, patient_dir, patientid):
        if patientid is None:
            patientid = os.path.basename(patient_dir)
        if patient_dir is None:
            patient_dir = '/home/mpschr/mounts/all_data/vbp_all/processed/exon/{}'.format(patientid)
        self.patientid = patientid
        return self.patientid

    def get_popvar_vcf(self):
        if not os.path.exists(self._popvar_vcf):
            self.run()
        return self._popvar_vcf

    def run(self):

        super().run(task_input=self.patient_dir, cmd_patient_name=None, exec_dir=".")
        logging.info(self.patientid)
        self.debug = False


        samples = self.time_points + [self.time_point_normal]
        bams = {tp: os.path.join(os.path.dirname(self.in_folder),
                                 self.patientid + "_" + tp,
                                 self.patientid + "_" + tp + "-ready.bam")
                for tp in samples}

        for bam in bams.values():
            assert os.path.isfile(bam)

        # pileup all variants
        mpileup_cmd_format = MPILEUP_CMD.format(**{'GENOME': self.conf_dict[GENOME],
                                                   'BED': self.conf_dict[BED_FILES[self._bed_name]],
                                                   'TUMOR_SAMPLES': ' '.join(bams[tp] for tp in self.time_points),
                                                   'NORMAL_SAMPLE': bams[self.time_point_normal],
                                                   'OUTPUT_VCF': self._popvar_vcf})

        return_code = ExternalCommand(mpileup_cmd_format, output_file=self._popvar_vcf).run()

        if return_code == 0:
            logging.info("Written file {}".format(self._popvar_vcf))
            logging.info("Task finished sucessfully in {}".format(self.timer.get_elapsed_time()))
        else:
            logging.error("Task finished with error after {}".format(self.timer.get_elapsed_time()))


        return self._popvar_vcf

    def folder_and_file_names(self, patient_dir):
        try:
            self.in_folder = get_mut_batch_folder(self.patientid, patient_dir, self.conf_dict[BCBIO_BACKUP])
        except FileNotFoundError:
            logging.warning('No bcbio folder found')
        self.out_folder = join(patient_dir, 'mutations')
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)
        filename_stem = join(self.out_folder, '{}_{}'.format(self.patientid, "{}"))
        self._popvar_vcf = filename_stem.format("popvar.{}.vcf".format(self._bed_name))
        logging.info("The output file will be: {}".format(self._popvar_vcf))

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            GENOME,
                                            BCBIO_BACKUP,
                                            POPVAR_EXACSNP_BED,
                                            POPVAR_COSMIC_BED,
                                            MUTREPORT_CGC])
        return conf_dict




def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-D', dest='patient_dir', required=True,
                        help='Directory with sample name. The master file should be '
                             'found at the following path: sample/mutations/sample_ID-ALLMUTS.tsv')


    parser.add_argument('-t', dest='time_points', action='append')
    parser.add_argument('-T', dest='time_point_normal', default='CR')
    parser.add_argument('--bed', dest='bed', choices=[EXAC, COSMIC])
    parser.add_argument('--debug', dest='debug', action='store_true', default=False)
    parser.add_argument('-M', dest='force_recreate_master', action='store_true', default=False)
    options = parser.parse_args()


    task = PatientSNPer(time_points=options.time_points, time_point_normal=options.time_point_normal,
                        patient_dir=os.path.abspath(options.patient_dir),
                        bed_file=BED_FILES[options.bed], bed_name=options.bed)
    task.run()


if __name__ == "__main__":  # detects if called from system
    cmdline()
