import argparse
from logging import basicConfig, DEBUG, info
import logging
import os
import sys
import re
from allpy.SimpleTask import SimpleTask, SampleTask
import pandas as pd

from allpy.config import pipeline_config
from allpy.config.config_constants import STAR_CORES, NOTIFY_EMAIL, STAR_BIN, CUTADAPT_BIN, SAMBAMBA_BIN, CORES, \
    CUTADAPT_ADAPTERS
from allpy.config.config_constants import STAR_GENOME_DIR
from allpy.externalcommands import ExternalCommand
from allpy.notify.notify import send_email
from allpy.exon.mutations.utils import Timer
from allpy.pipelineutils.filepaths import WorkDirs, ALIGN, TRIM, OutFiles

STAR_OVERHANG = 'stargenome.LOG.out_overhang'

CMD_CUTADAPT_PE = "{CUTADAPT_BIN} " \
      "--times=2 --quality-base=33 --quality-cutoff=5 --format=fastq " \
      "{ADAPTERS} " \
      " --minimum-length=25  " \
      "-o {WORK_DIR}/{SAMPLE}_R1.trimmed.tmp.fastq.gz " \
      "-p {WORK_DIR}/{SAMPLE}_R2.trimmed.tmp.fastq.gz " \
      "{INPUT_PATH}/{SAMPLE}_R1.fastq.gz {INPUT_PATH}/{SAMPLE}_R2.fastq.gz;" \
      "{CUTADAPT_BIN} " \
      "--times=2 --quality-base=33 --quality-cutoff=5 --format=fastq " \
      "{ADAPTERS} " \
      " --minimum-length=25  " \
      "-o {WORK_DIR}/{SAMPLE}_R2.trimmed.fastq.gz " \
      "-p {WORK_DIR}/{SAMPLE}_R1.trimmed.fastq.gz " \
      "{WORK_DIR}/{SAMPLE}_R2.trimmed.tmp.fastq.gz " \
      "{WORK_DIR}/{SAMPLE}_R1.trimmed.tmp.fastq.gz; " \
      "rm {WORK_DIR}/{SAMPLE}_R1.trimmed.tmp.fastq.gz {WORK_DIR}/{SAMPLE}_R2.trimmed.tmp.fastq.gz "


CMD_STAR_2PASS ="{STAR_BIN} " \
          "--genomeDir {GENOME_DIR} " \
          "--twopassMode Basic --sjdbOverhang {STAR_OVERHANG} " \
          "--readFilesIn {WORK_DIR_TRIM}/{SAMPLE}_R1.trimmed.fastq.gz {WORK_DIR_TRIM}/{SAMPLE}_R2.trimmed.fastq.gz " \
          "--runThreadN {CORES} --outFileNamePrefix {WORK_DIR_ALIGN}/2PASS_ " \
          "--outReadsUnmapped Fastx --outFilterMultimapNmax 10 " \
          "--outStd SAM --outSAMunmapped Within " \
          "--outSAMattributes NH HI NM MD AS --readFilesCommand zcat  " \
          "--outSAMattrRGline ID:3 PL:illumina PU:2PASS_ SM:{SAMPLE} " \
          "--outSAMstrandField intronMotif " \
          " | samtools sort -@ {CORES} -m 1G  " \
          "-T {WORK_DIR_ALIGN}/{SAMPLE}-sorttmp " \
          "-o {WORK_DIR_ALIGN}/{SAMPLE}.bam " \
          "/dev/stdin"


CMD_STAR_1PASS ="{STAR_BIN} " \
          "--genomeDir {GENOME_DIR} " \
          "--twopassMode Basic --sjdbOverhang {STAR_OVERHANG} " \
          "--readFilesIn {WORK_DIR_TRIM}/{SAMPLE}_R1.trimmed.fastq.gz {WORK_DIR_TRIM}/{SAMPLE}_R2.trimmed.fastq.gz " \
          "--runThreadN {CORES} --outFileNamePrefix {WORK_DIR_ALIGN}/1PASS_ " \
          "--outReadsUnmapped Fastx --outFilterMultimapNmax 10 " \
          "--outStd SAM --outSAMunmapped Within " \
          "--outSAMattributes NH HI NM MD AS --readFilesCommand zcat  " \
          "--outSAMattrRGline ID:3 PL:illumina PU:1PASS_ SM:{SAMPLE} " \
          "--outSAMstrandField intronMotif " \
          " | samtools sort -@ {CORES} -m 1G  " \
          "-T {WORK_DIR_ALIGN}/{SAMPLE}-sorttmp " \
          "-o {WORK_DIR_ALIGN}/{SAMPLE}.bam " \
          "/dev/stdin"

CMD_INDEX = "{SAMBAMBA_BIN} index -t {CORES} {WORK_DIR_ALIGN}/{SAMPLE}.bam"


__author__ = 'mpschr'


#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


class TrimAndStar(SampleTask):

    def __init__(self):
        super().__init__('TRIM_AND_ALIGN')

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            sample = cmd_sample_name
        else:
            sample = os.path.basename(os.path.abspath(task_input))
        return sample

    def run(self, input_dir, exec_dir='.', sample_name=None, passes=2, trim=True):
        super().run(input_dir, exec_dir, sample_name)

        self.check_fastqs_exist(input_dir, self.sample)

        workdir_trim = self.work_dirs.get(TRIM)
        workdir_align = self.work_dirs.get(ALIGN)

        conf_dict = self.load_configuration()

        os.chdir(exec_dir)


        if trim:
            self.run_cutadapt(conf_dict, input_dir, self.sample, workdir_trim)
        else:
            workdir_trim = input_dir

        c = self.run_star(conf_dict, self.sample, workdir_align, workdir_trim, passes, trim)

        elapsed_time = self.timer.get_elapsed_time()
        logging.debug(self.timer.get_elapsed_time())

        if c == 0:
            send_email(conf_dict['DEFAULT.notify_email'],"{} STAR-ALIGN done".format(self.sample), elapsed_time)
        else:
            send_email(conf_dict['DEFAULT.notify_email'],"FAILURE of {} STAR-ALIGN".format(self.sample), elapsed_time)

        return c

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([STAR_CORES,
                                            STAR_GENOME_DIR,
                                            NOTIFY_EMAIL,
                                            STAR_BIN,
                                            CUTADAPT_BIN,
                                            CUTADAPT_ADAPTERS,
                                            SAMBAMBA_BIN,
                                            CORES])
        conf_dict[STAR_OVERHANG] = self.read_sjdboverhang(conf_dict)
        return conf_dict

    def read_sjdboverhang(self, conf_dict):
        data_file = open(os.path.join(conf_dict[STAR_GENOME_DIR], 'Log.out'))

        params=None
        for line in data_file:
            if 'sjdbOverhang' in line:
                params =  line.rstrip().split()
                print(params)
                break

        data_file.close()

        return params[1]

    def run_star(self, conf_dict, sample, workdir_align, workdir_trim, passes, trim):
        info('working in dir: {}'.format(os.path.abspath(workdir_align)))
        info('RUNNING STAR\n=========================')
        cmd = self.get_star_command(conf_dict, sample, workdir_align, workdir_trim, passes, trim)
        starcommand = ExternalCommand(cmd, output_file=self.output_files.get(ALIGN))
        c = starcommand.run()

        cmd2 = CMD_INDEX.format(**{'SAMPLE': sample,
                                    'WORK_DIR_ALIGN': workdir_align,
                                    'SAMBAMBA_BIN': conf_dict[SAMBAMBA_BIN],
                                    'CORES': conf_dict[CORES]})
        c2 = ExternalCommand(cmd2, output_file="{}/{}.bam.bai".format(workdir_align, sample)).run()

        return c

    def get_star_command(self, conf_dict, sample, workdir_align, workdir_trim, passes, trim):

        if passes == 2:

            cmd = CMD_STAR_2PASS.format(**{'SAMPLE': sample,
                                           'STAR_BIN': conf_dict[STAR_BIN],
                                           'STAR_OVERHANG': conf_dict[STAR_OVERHANG],
                                           'GENOME_DIR': conf_dict[STAR_GENOME_DIR],
                                           'WORK_DIR_TRIM': workdir_trim,
                                           'WORK_DIR_ALIGN': workdir_align,
                                           'CORES': conf_dict[STAR_CORES]})

        elif passes == 1:

            cmd = CMD_STAR_1PASS.format(**{'SAMPLE': sample,
                                           'STAR_BIN': conf_dict[STAR_BIN],
                                           'GENOME_DIR': conf_dict[STAR_GENOME_DIR],
                                           'WORK_DIR_TRIM': workdir_trim,
                                           'WORK_DIR_ALIGN': workdir_align,
                                           'CORES': conf_dict[STAR_CORES]})
        else:
            raise RuntimeError('Incorrect number of passes specified for STAR')

        if not trim:
            cmd = cmd.replace('.trimmed', '')

        return cmd

    def run_cutadapt(self, conf_dict, input_dir, sample, workdir_trim):
        info('working in dir: {}'.format(os.path.abspath(workdir_trim)))
        info('RUNNING CUTADAPT\n=========================')
        cmd = CMD_CUTADAPT_PE.format(**{'SAMPLE': sample,
                                        'ADAPTERS': ' '.join( "--adapter={}".format(x) for x in conf_dict[CUTADAPT_ADAPTERS] ),
                                        'CUTADAPT_BIN': conf_dict[CUTADAPT_BIN],
                                        'WORK_DIR': workdir_trim,
                                        'INPUT_PATH': input_dir})

        c = ExternalCommand(cmd, output_file="{}/{}_R1.trimmed.fastq.gz".format(workdir_trim, sample)).run()
        assert c == 0

    def check_fastqs_exist(self, input_dir, sample):
        fastq1 = "{0}/{1}_R1.fastq.gz".format(input_dir, sample)
        fastq2 = "{0}/{1}_R2.fastq.gz".format(input_dir, sample)
        missing_fastqs = []
        if not os.path.exists(fastq1):
            missing_fastqs.append(fastq1)
        if not os.path.exists(fastq2):
            missing_fastqs.append(fastq2)
        if len(missing_fastqs) > 0:
            raise RuntimeError('Following fastqs not found: {}'.format(missing_fastqs))


#bcbio_nextgen.py config/trios/AL02.yaml --workdir scratch_AL02 -t local -n 20
#screen -S testy -dm sh -c 'ls -al; echo $?; exec bash -i'

def cmdline():
    parser = argparse.ArgumentParser(description='Run cutadatpt and STAR aligner for one sample\n\n'
                                                 'The program will create a scratch dir with the intermediate files'
                                                 'and the final .bam file inside.')

    parser.add_argument('-i', dest='input_dir', required=True, help='The directory with the two fastq files of the sample')
    parser.add_argument('-x', dest='exec_dir', default='.', required=False, help='Where the program should be run and results stored')
    parser.add_argument('-s', dest='sample_name', required=False, help='Only needed to specify if the samplename should be '
                                                                       'different than the folder name containing the fastq files')
    parser.add_argument('--passes', dest='passes', default=2, help='Run STAR aligner with 1 pass or 2 pass - default 2',
                            type=int, choices=[1, 2])
    parser.add_argument('--no-trim', dest='no_trim', action='store_true')
    #parser.add_argument('-adapters')
#    parser.add_argument('-n', dest='parallel_tasks', type=int, default=20)

    options = parser.parse_args()
    task = TrimAndStar()

    task.run(options.input_dir, options.exec_dir, options.sample_name, options.passes, trim=not options.no_trim)

if __name__ == "__main__":  # detects if called from system
    cmdline()
