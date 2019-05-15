import argparse
import logging
import os
from os.path import join
from allpy.SimpleTask import SampleTask
from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL, DEFUSE_CORES, DEFUSE_CONFIG, DEFUSE_DATA_DIR
from allpy.externalcommands import ExternalCommand

from allpy.notify.notify import send_email

CR = "CR"


CMD_FASTA =  "gunzip -c {INPUT_DIR}/{SAMPLE}_R1.fastq.gz > {SCRATCH_DIR}/{SAMPLE}_R1.fastq; gunzip -c {INPUT_DIR}/{SAMPLE}_R2.fastq.gz > {SCRATCH_DIR}/{SAMPLE}_R2.fastq "
CMD_DEFUSE = "defuse_run.pl -d {DEFUSE_DATA_DIR} -c {DEFUSE_CONFIG} -1 {SCRATCH_DIR}/{SAMPLE}_R1.fastq -2 {SCRATCH_DIR}/{SAMPLE}_R2.fastq -o {SCRATCH_DIR} --rescla --resfil --parallel {CORES} "
CMD_CLEANUP = "rm -r {SCRATCH_DIR}/jobs  {SCRATCH_DIR}/reads.*  {SCRATCH_DIR}/cdna.pair.sam {SCRATCH_DIR}/{SAMPLE}_R1.fastq {SCRATCH_DIR}/{SAMPLE}_R2.fastq; mv {SCRATCH_DIR}/* {OUTPUT_FOLDER}/ ; rmdir --ignore-fail-on-non-empty -p {SCRATCH_DIR}"

logger = logging.getLogger()


class DefuseRunner(SampleTask):

    def __init__(self):
        super().__init__("DEFUSE")

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            DEFUSE_CORES,
                                            DEFUSE_DATA_DIR,
                                            DEFUSE_CONFIG
					])
        return conf_dict

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            sample = cmd_sample_name
        else:
            sample = os.path.basename(os.path.abspath(task_input))
        return sample

    def run(self, patient_dir=None, patient_id=None, sample_dir=None, scratch_dir=None,
            output_dir=None):
        assert patient_dir != sample_dir

        if output_dir is None:
            output_dir = os.getcwd()

        super().run(sample_dir,
                    output_dir,
                    patient_id=patient_id)

        conf_dict = self.load_configuration()


        del output_dir

        res_file = self.output_files.generate('defuse', 'results.tsv')
        if os.path.isfile(res_file):
            logger.warning("FILE already exists {}. \n\t EXITING".format(res_file))
            return

        sample_out_dir = self.work_dirs.get('defuse')

        if scratch_dir is None:
            scratch_dir = self.work_dirs.get('defuse', scratch=True)

        commands = []

        commands.append(
            ExternalCommand(
                CMD_FASTA.format(**{
                    "INPUT_DIR": sample_dir,
                    "SCRATCH_DIR": scratch_dir,
                    "SAMPLE": self.sample
                    }),
                output_file=os.path.join(scratch_dir, self.sample + "_R2.fastq")
            )
        )

        commands.append(
            ExternalCommand(
                CMD_DEFUSE.format(**{
                    "CORES": conf_dict[DEFUSE_CORES],
                    "DEFUSE_CONFIG": conf_dict[DEFUSE_CONFIG],
                    "DEFUSE_DATA_DIR": conf_dict[DEFUSE_DATA_DIR],
                    "SCRATCH_DIR": scratch_dir,
                    "SAMPLE": self.sample
                    })
            )
        )


        commands.append(
            ExternalCommand(
                CMD_CLEANUP.format(**{
                    "SCRATCH_DIR": scratch_dir,
                    "SAMPLE": self.sample,
                    "OUTPUT_FOLDER": sample_out_dir
                    })
            )
        )

        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

        for external_command in commands:
            c = external_command.run()
            elapsed_time = self.timer.get_elapsed_time()
            if c != 0:
                send_email('michael.p.schroeder@gmail.com',"FAILURE of {}  DEFUSE".format(self.sample), elapsed_time)
                return

        logging.info(elapsed_time)

        send_email('michael.p.schroeder@gmail.com',"{} DEFUSE done".format(self.sample), elapsed_time)


def absolute_path(path):
    if path is None:
        return None

    return os.path.abspath(path)

def cmdline_rna():
    cmdline(seqtype='rna')

def cmdline(seqtype = 'exon'):

    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-D', dest='patient_dir', default=None,
                        help="The output folder of the patient, containing one folder for each sample")
    group.add_argument('-S', dest='sample_dir', default=None,
                        help="Sample folder, containing the 'align' folder with the aligned .bam file")

    parser.add_argument('-p', dest='patientid', required=False)

    parser.add_argument('--output-dir', dest='output_dir', required=True,
                        help = "base dir of where the final outputs should be stored")

    parser.add_argument('--scratch-dir', dest='scratch_dir', default=None, required=False)
    #parser.add_argument('-S', dest='sequencetype', choices=['exon', 'rna'], default='exon')

    options = parser.parse_args()

    task = DefuseRunner()

    if not options.patient_dir is None:
        sampledirs = []
        patientdir = absolute_path(options.patient_dir)
        for dir in os.listdir(patientdir):
            if os.path.isdir(join(patientdir, dir)):
                sampledirs.append(join(patientdir,dir))

        if options.patientid is None:
            options.patientid = os.path.basename(options.patient_dir)
    else:
        sampledirs = [absolute_path(options.sample_dir)]

    logger.info("Found {} sample dirs: \n\t\t - {}".format(len(sampledirs), "\n\t\t - ".join(sampledirs)))

    for sdir in sampledirs:
        task.run(absolute_path(options.patient_dir), options.patientid, sdir, scratch_dir=options.scratch_dir,
                 output_dir=options.output_dir)

if __name__ == "__main__":  # detects if called from system
    cmdline()
