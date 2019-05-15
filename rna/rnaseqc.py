import argparse
import logging
import os
import tempfile
from os.path import join
from allpy.SimpleTask import PatientTask
from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL, TRANSCRIPT_GTF, RNASEQC_JAR, GENOME
from allpy.externalcommands import ExternalCommand

from allpy.notify.notify import send_email

CR = "CR"

CMD_RNA = "{PIPEFAIL} java -Xmx6g -Xms256m -XX:-UseParallelGC -jar {RNASEQC_JAR} rnaseq -o {OUTDIR} -r {GENOME} -s '{SAMPLE}|{BAM}|Notes' -t {GTF} -ttype 2 {MV_SCRATCH}"

logger = logging.getLogger()


class RNAseqcRunner(PatientTask):

    def __init__(self):
        super().__init__("RNASEQ-C")

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([RNASEQC_JAR,
                                            TRANSCRIPT_GTF,
                                            NOTIFY_EMAIL,
                                            GENOME])
        return conf_dict

    def derive_patient_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            patient = cmd_sample_name
        else:
            patient = os.path.basename(os.path.abspath(task_input))
        return patient

    def run(self, patient_dir=None, patientid=None, sample_dir=None, scratch_dir=None,alt_gtf=None):
        assert patient_dir != sample_dir

        super().run(patient_dir if patient_dir is not None else sample_dir,
                    patient_dir if patient_dir is not None else sample_dir,
                    patientid)

        conf_dict = self.load_configuration()
        if not alt_gtf is None:
            conf_dict[TRANSCRIPT_GTF] = alt_gtf

        samples = {}
        if sample_dir == None:
            for sdir in os.listdir(patient_dir):
                samples[sdir] = os.path.join(patient_dir, sdir)
        else:
            samples[os.path.basename(sample_dir)] = sample_dir

        jobs = {}

        for sampleid, sdir in samples.items():
            print(sdir)
            sample_bam_dir = join(sdir, 'align')
            bamfile = "{0}/{1}.bam".format(sample_bam_dir, sampleid)

            if not os.path.exists(bamfile):
                logging.warning("BAM file not found: {}\n-------------------------\n".format(bamfile))
                continue
            else:
                logging.info("BAM file: {}\n-------------------------\n".format(bamfile))

            outdir = join(sample_bam_dir, 'qc', 'rnaseqc')
            if scratch_dir is None:
                scratch_dir_temp = outdir
            else:
                os.makedirs(scratch_dir, exist_ok=True)
                prefix = join(os.path.expanduser(scratch_dir), 'TMP.seqc-' + sampleid)
                scratch_dir_temp = tempfile.mkdtemp(prefix = prefix)

            cmd = CMD_RNA.format(**{
                'PIPEFAIL': "set -euxo pipefail; ",
                'RNASEQC_JAR': conf_dict[RNASEQC_JAR],
                'OUTDIR': scratch_dir_temp,
                'BAM': bamfile,
                #'GTF': conf_dict[TRANSCRIPT_GTF].replace('.gtf','.RNASEQC.gtf'),
                'GTF': conf_dict[TRANSCRIPT_GTF],
                'GENOME': conf_dict[GENOME],
                'SAMPLE': sampleid,
                'MV_SCRATCH': '; mv {} {}'.format(scratch_dir_temp, outdir) if scratch_dir_temp != outdir else "; "
                })

            if not os.path.exists(join(sample_bam_dir, 'qc')):
                os.mkdir(join(sample_bam_dir, 'qc'))

            logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

            print(cmd)
            c = ExternalCommand(cmd, output_file=outdir).run()
            jobs[sampleid] = c
            if c != 0:
                break


        returncodes = {}
        for p, c in jobs.items():
            logging.debug('Return code of {} is: {}'.format(p, c))
            returncodes[p] = c

        elapsed_time = self.timer.get_elapsed_time()
        logging.info(elapsed_time)

        reportid = None
        if patientid is None and patient_dir is not None:
            reportid = os.path.basename(patient_dir)
        else:
            reportid = list(samples.keys())[0]

        if max(returncodes.values()) == 0 and min(returncodes.values()) == 0:
            send_email('michael.p.schroeder@gmail.com',"{} RNA-SEQ-C done".format(reportid), elapsed_time)
        else:
            send_email('michael.p.schroeder@gmail.com',"FAILURE of {}  RNA-SEQ-C".format(reportid), "{}\n\n{}".format(elapsed_time,returncodes))

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
    parser.add_argument('-G', dest='alt_gtf', required=False,
                        help='alternative GTF file, as RNASEQc can be picky about lines without "transcript_id"')

    parser.add_argument('--scratch-dir', dest='scratch_dir', default=None, required=False)
    #parser.add_argument('-S', dest='sequencetype', choices=['exon', 'rna'], default='exon')

    options = parser.parse_args()
    print(options)


    task = RNAseqcRunner()
    task.run(absolute_path(options.patient_dir), options.patientid, absolute_path(options.sample_dir),
             scratch_dir=options.scratch_dir,alt_gtf=options.alt_gtf)

if __name__ == "__main__":  # detects if called from system
    cmdline()
