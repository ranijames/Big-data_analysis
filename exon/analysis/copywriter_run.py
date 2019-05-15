import argparse
from logging import basicConfig, DEBUG, info, debug
import logging
import os
import sys
import subprocess
import pkg_resources
from allpy.SimpleTask import PatientTask
from allpy.config import pipeline_config
from allpy.config.config_constants import BCBIO_BACKUP, NOTIFY_EMAIL

from allpy.exon.mutations.utils import Timer, get_mut_batch_folder
from allpy.externalcommands import ExternalCommand
from allpy.notify.notify import send_email

REL = 'REL'

ID = 'ID'

CMD_OUTER = "screen -S {1}_cpwr -dm sh -c '{0} ; " \
            " selfnotify \"{1} COPYWRITER done\" ; echo $?; exec bash -i' "
# " echo $?; exec bash -i' "

# retrieve absolute paths of R scripts from the source package
CMD_CPWR = "Rscript " + pkg_resources.resource_filename('allpy', 'copywriteRun.R') + " -p {0} -o {1} -i {2} -r {3} -c {4} -t {5} -g {6} -b {7}"
CMD_CPWR_one_tumor_sample = "Rscript  " + pkg_resources.resource_filename('allpy', 'copywriteRun_idonly.R') + " -p {0} -o {1} -i {2} -c {3} -t {4} -g {5} -b {6}"

CMD_MV = "cd {1} ; touch {2}/{0}.info ; mv {2} {3} ; ln -s {3} {2} "

COVERED_BED = '/home/mpschr/Data/agilent_sure_select_v5+utr/S04380219_Covered.bed'
GENES_BED = '/home/mpschr/bin/bcbionextgen/data/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.bed'

# logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


class CopywriteRunner(PatientTask):
    def __init__(self, patient_dir, patientid=None):
        self.conf_dict = pipeline_config.Config().get_config_dict([BCBIO_BACKUP, NOTIFY_EMAIL])
        self.patient_dir = patient_dir
        self.patient = self.derive_patient_name(patient_dir, patientid)

    def derive_patient_name(self, task_input, cmd_patient_name):
        if cmd_patient_name is None:
            return os.path.basename(task_input)
        else:
            return cmd_patient_name

    def run(self, out_dir=None, covered_regions_bed=COVERED_BED, genes_bed=GENES_BED, bin_size='20kb',
            normal_sample="CR", time_points=None):

        timer = Timer(start_now=True)

        if out_dir is None:
            out_dir = self.patient_dir

        if len(self.patient) < 4:
            raise RuntimeError("too short patientid: {}".format(self.patient))

        if time_points is None:
            self.time_points = [ID, REL]
        else:
            self.time_points = time_points


        bcbio_folder = os.path.dirname(get_mut_batch_folder(self.patient, self.patient_dir, self.conf_dict[BCBIO_BACKUP]))

        id_file = os.path.join(bcbio_folder, self.patient + '_' + self.time_points[0], self.patient + '_' + '{}-ready.bam'.format(self.time_points[0]))
        if len(self.time_points) == 2:
            rel_file = os.path.join(bcbio_folder, self.patient + '_' + self.time_points[1], self.patient + '_' + '{}-ready.bam'.format(self.time_points[1]))
        cr_file = os.path.join(bcbio_folder, self.patient + '_' + normal_sample, self.patient + '_{}-ready.bam'.format(normal_sample))

        out_dir = os.path.join(out_dir, 'copywriteR')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        elif os.path.exists(os.path.join(out_dir, "CNAprofiles")) and os.path.islink(os.path.join(out_dir, "CNAprofiles")):
            logging.info("Unlinking CNAprofiles")
            os.unlink(os.path.join(out_dir, "CNAprofiles"))

        bin_numeric = bin_size.replace('kb', '000')
        cmd = ''
        if len(self.time_points) == 2:
            cmd = CMD_CPWR.format(self.patient, out_dir, id_file, rel_file, cr_file, covered_regions_bed, genes_bed, bin_numeric)
        else:
            cmd = CMD_CPWR_one_tumor_sample.format(self.patient, out_dir, id_file, cr_file, covered_regions_bed, genes_bed, bin_numeric)
        oldprofiles = 'CNAprofiles'
        newprofiles = 'CNAprofiles_' + bin_size
        cmd_cleaning = CMD_MV.format(bin_size, out_dir, oldprofiles, newprofiles)
        returncode = ExternalCommand(" ; ".join(["set -e",cmd,cmd_cleaning]), output_file=os.path.join(out_dir, newprofiles)).run()

        elapsed_time = timer.get_elapsed_time()
        logging.info(elapsed_time)
        loginfo = "\n"
        try:
            logfile = os.path.join(out_dir, newprofiles, "CopywriteR.log")
            logproc = subprocess.Popen("grep -i mad {} |  head -3 | sed 's/^INFO.\+]//'".format(logfile), stdout=subprocess.PIPE, shell=True)
            loginfo += logproc.communicate()[0].decode()
        except Exception as e:
            logging.debug("log-error",e)

        if returncode == 0:
            send_email('michael.p.schroeder@gmail.com', "{} COPYWRITER {} done".format(self.patient, bin_size), elapsed_time+loginfo)
        else:
            send_email('michael.p.schroeder@gmail.com', "FAILURE of {} COPYWRITER {}".format(self.patient, bin_size), elapsed_time)
        return returncode


def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-D', dest='sample_dir', required=True, help='Where the output folder from bcbio is, the root '
                                                                     'for the sample results')
    parser.add_argument('-t', dest='time_points', action='append')

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-b', dest='bin_size', default='20kb', choices=['20kb', '30kb', '40kb', '50kb', '80kb','100kb'],
                        help='The bin size for the read count data. Default 20kb')
    parser.add_argument('-o', dest='out_dir', required=False, help='Where the program results should stored. Specify'
                                                                   ' if different than bam_root_dir.')
    parser.add_argument('--normal-sample', dest='normal_sample', default='CR', help='indicate if different than CR')

    options = parser.parse_args()
    task = CopywriteRunner(os.path.abspath(options.sample_dir), options.patientid)
    task.run(options.out_dir,  bin_size=options.bin_size, normal_sample=options.normal_sample, time_points=options.time_points)


if __name__ == "__main__":  # detects if called from system
    cmdline()
