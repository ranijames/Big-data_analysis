import sys
import argparse
import logging
import os
import subprocess
from os.path import join
from allpy.SimpleTask import SampleTask
from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL, CORES, QUALIMAP_BIN, TRANSCRIPT_GTF, TARGET_REGIONS
from allpy.exon.cna.cna_report import ID, REL

from allpy.notify.notify import send_email

CR = "CR"



CMD_EXON = "unset DISPLAY && {QUALIMAP_BIN} bamqc -bam {BAM} -gff {TARGET_BED} -nt 8 -c -outdir {OUTDIR}  -nr 5000  --java-mem-size=15G ;"
CMD_RNA = "unset DISPLAY && {QUALIMAP_BIN} rnaseq -outdir {OUTDIR} -a proportional -bam {BAM} -gtf {GTF} --paired --java-mem-size=15G -oc {COUNTS_FILE};"
    #unset DISPLAY && /home/mpschr/bin/bcbionextgen/tools/bin/qualimap rnaseq -outdir {OUTDIR}
    # -a proportional -bam {bam} -gtf /home/mpschr/bin/bcbionextgen/data/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf --java-mem-size=110g


logger = logging.getLogger()



class QualimapRunner(SampleTask):

    def __init__(self):
        super().__init__("QUALIMAP")

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([QUALIMAP_BIN,
                                            TRANSCRIPT_GTF,
                                            TARGET_REGIONS,
                                            NOTIFY_EMAIL,
                                            CORES])
        return conf_dict

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            sample = cmd_sample_name
        else:
            sample = os.path.basename(os.path.abspath(task_input))
        return sample

    def run(self, patient_dir, patientid=None, time_points=None, exonseq=False, rnaseq=False, panel=False):
        super().run(patient_dir, patient_dir, patientid)

        if time_points == None or time_points == []:
            time_points = [ID, REL, CR]
            if rnaseq:
                time_points = [ID, REL]

        jobs = {}

        conf_dict = self.load_configuration()

        if patientid is None:
            patientid = os.path.basename(patient_dir)


        for t in time_points:
            sample_bam_dir = None
            bamfile = None
            outfile = None
            if exonseq or panel:
                folder = 'bcbio' if exonseq else 'panel'
                target_regions_bed = conf_dict[TARGET_REGIONS]['exon' if exonseq else 'panel']
                sample_bam_dir = join(patient_dir, folder, '{}_{}'.format(patientid, t))
                bamfile = "{0}/{1}_{2}-ready.bam".format(sample_bam_dir, patientid, t)
                cmd = CMD_EXON.format(**{
                    'QUALIMAP_BIN': conf_dict[QUALIMAP_BIN],
                    'OUTDIR': join(sample_bam_dir, 'qc', 'qualimap'),
                    'TARGET_BED': target_regions_bed,
                    'BAM': bamfile
                })
            elif rnaseq:
                sample_bam_dir = join(patient_dir, '{}_{}'.format(patientid, t), 'align')
                bamfile = "{0}/{1}_{2}.bam".format(sample_bam_dir, patientid, t)
                countfile = bamfile.replace('.bam', 'qualimap-genecounts.txt')
                cmd = CMD_RNA.format(**{
                    'QUALIMAP_BIN': conf_dict[QUALIMAP_BIN],
                    'OUTDIR': join(sample_bam_dir, 'qc', 'qualimap'),
                    'BAM': bamfile,
                    'GTF': conf_dict[TRANSCRIPT_GTF],
                    'COUNTS_FILE': countfile
                })
                outfile = countfile

            if not os.path.exists(bamfile):
                logging.warning("BAM file not found: {}".format(bamfile))
                continue
            else:
                logging.info("BAM file: {}".format(bamfile))

            #logging:
            if not os.path.exists(join(sample_bam_dir, 'qc')):
                os.mkdir(join(sample_bam_dir, 'qc'))

            logfilename = join(sample_bam_dir, 'qc', '{}_{}-qualimap.log'.format(patientid, t))
            log_file = open(logfilename, 'w+')
            logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

            logging.info('Full command: "{}"'.format(cmd))
            logging.info('See log file: {}'.format(logfilename))
            c = subprocess.Popen(cmd, shell=True, stderr=log_file, stdout=log_file)
            jobs["{}_{}".format(patientid, t)] = c


        returncodes = {}
        for p, j in jobs.items():
            j.wait()
            logging.debug('Return code of {} is: {}'.format(p, j.returncode))
            returncodes[p] = j.returncode

        elapsed_time = self.timer.get_elapsed_time()
        logging.info(elapsed_time)

        if max(returncodes.values()) == 0 and min(returncodes.values()) == 0:
            send_email('michael.p.schroeder@gmail.com',"{} Qualimap done".format(patientid), elapsed_time)
        else:
            send_email('michael.p.schroeder@gmail.com',"FAILURE of {} Qualimap".format(patientid), "{}\n\n{}".format(elapsed_time,returncodes))


def cmdline_exon():
    cmdline(seqtype='exon')

def cmdline_panel():
    cmdline(seqtype='panel')

def cmdline_rna():
    cmdline(seqtype='rna')

def cmdline(seqtype = 'exon'):

    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-D', dest='patient_dir', required=True, help='Where the output folder from bcbio is, the root '
                                                                       'for the sample results')
    parser.add_argument('-t', dest='time_points', action='append')
    #parser.add_argument('-S', dest='sequencetype', choices=['exon', 'rna'], default='exon')

    options = parser.parse_args()
    task = QualimapRunner()
    task.run(os.path.abspath(options.patient_dir), options.patientid, options.time_points,
             seqtype == 'exon', seqtype == 'rna', seqtype == 'panel')

if __name__ == "__main__":  # detects if called from system
    cmdline()
