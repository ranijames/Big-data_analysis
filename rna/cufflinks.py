import argparse
from logging import basicConfig, DEBUG, info
import logging
import os
import sys
from allpy.SimpleTask import SampleTask

from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL, CORES, \
    GENOME, TRANSCRIPT_GTF, TRANSCRIPT_GTF_MASK, CUFFLINKS_BIN
from allpy.externalcommands import ExternalCommand
from allpy.notify.notify import send_email
from allpy.pipelineutils.filepaths import CUFFLINKS, CUFFLINKS_NOVEL

CMD_CUFFLINKS = """{CUFFLINKS_BIN} \
 --num-threads {CORES} \
 --quiet \
 --no-update-check \
 --max-bundle-frags 2000000 \
 --library-type fr-unstranded \
 --frag-bias-correct {GENOME} \
 --multi-read-correct \
 --GTF{GUIDE_OR_NO} {TRANSCRIPT_GTF} \
 --mask-file {TRANSCRIPT_GTF_MASK} \
 --output-dir {WORK_DIR} \
 {BAM_FILE}
"""
__author__ = 'ARJ'

#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


class Cufflinks(SampleTask):

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            sample = cmd_sample_name
        else:
            sample = os.path.basename(os.path.abspath(task_input)).split('.bam')[0]
        return sample

    def __init__(self):
        super().__init__('CUFFLINKS')

    def run(self, input_bam, exec_dir='.', sample_name=None, detect_unannotated=False):
        super().run(input_bam, exec_dir, sample_name)

        if detect_unannotated:
            workdir_cufflinks = self.work_dirs.get_popvar_vcf(CUFFLINKS_NOVEL)
        else:
            workdir_cufflinks = self.work_dirs.get_popvar_vcf(CUFFLINKS)

        conf_dict = self.load_configuration()

        os.chdir(exec_dir)

        c = self.run_cufflinks(conf_dict, self.sample, input_bam, workdir_cufflinks, detect_unannotated)

        elapsed_time = self.timer.get_elapsed_time()
        logging.debug(self.timer.get_elapsed_time())

        runtype = "normal" if not detect_unannotated else "guide"
        if c == 0:
            send_email(conf_dict['DEFAULT.notify_email'],"{} {} () done".format(self.sample, self.task_name, runtype), elapsed_time)
        else:
            send_email(conf_dict['DEFAULT.notify_email'],"FAILURE of {} CUFFLINKS {}".format(self.sample, self.task_name, runtype), elapsed_time)

        return c

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            CORES,
                                            GENOME,
                                            TRANSCRIPT_GTF, TRANSCRIPT_GTF_MASK,
                                            CUFFLINKS_BIN
                                            ])
        return conf_dict

    def run_cufflinks(self, conf_dict, sample, input_bam, workdir_cufflinks, detect_unannotated):
        cmd = CMD_CUFFLINKS.format(**{'SAMPLE': sample,
                                        'BAM_FILE': input_bam,
                                        'CUFFLINKS_BIN': conf_dict[CUFFLINKS_BIN],
                                        'WORK_DIR': workdir_cufflinks,
                                        'GENOME': conf_dict[GENOME],
                                        'GUIDE_OR_NO': "" if not detect_unannotated else "-guide",
                                        'TRANSCRIPT_GTF': conf_dict[TRANSCRIPT_GTF],
                                        'TRANSCRIPT_GTF_MASK': conf_dict[TRANSCRIPT_GTF_MASK],
                                      'CORES': conf_dict[CORES]
                                      })

        c = ExternalCommand(cmd, output_file="{}/genes.fpkm_tracking".format(workdir_cufflinks), touch_output=True).run()
        if c != 0:
            os.remove(workdir_cufflinks)
            logging.warning("Removed output dir due to error")
        return c



#bcbio_nextgen.py config/trios/AL02.yaml --workdir scratch_AL02 -t local -n 20
#screen -S testy -dm sh -c 'ls -al; echo $?; exec bash -i'

def cmdline():
    parser = argparse.ArgumentParser(description='Run cufflinks for one sample\n\n'
                                                 'The program will use/create a scratch dir for the intermediate files'
                                                 'and the final cufflink files inside.')

    parser.add_argument('-i', dest='input_bam', required=True, help='The bam file to be used')
    parser.add_argument('--novel', dest='novel', action='store_true', default=False,
                        required=False, help='Detect novel transcripts annotation (intensive)')
    parser.add_argument('-x', dest='exec_dir', default='.', required=False, help='Where the program should be run and results stored')
    parser.add_argument('-s', dest='sample_name', required=False, help='Only needed to specify if the samplename should be '
                                                                       'different than the bam file name')

    options = parser.parse_args()
    task = Cufflinks()
    task.run(options.input_bam, options.exec_dir, options.sample_name, detect_unannotated=options.novel)

if __name__ == "__main__":  # detects if called from system
    cmdline()
