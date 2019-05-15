import argparse
import logging
import os
import sys
import subprocess
import re
from os.path import join

#from allpy.exon.mutations.gemini_extraction import get_sample_muts
from allpy.exon.analysis.bcbio_run import BcbioRunner
from allpy.exon.analysis.copywriter_run import CopywriteRunner
from allpy.exon.analysis.sciclone_run import SciCloneRunner
from allpy.exon.mutations.plots import do_plots

logger = logging.getLogger()

class ExonAnalysisRun():

    def __init__(self):
        pass

    def do_run_task(self, task_name, last_return_code, verification_file):
        if os.path.exists(verification_file):
            return False

        if last_return_code == 0:
            return True

        if last_return_code != 0:
            logger.error("Exiting {} because of job failure".format(task_name))
            exit(-1)

    def run(self, patientid):

        bcbio = BcbioRunner()
        copywriter = CopywriteRunner(patientid=patientid)
        sciclone = SciCloneRunner()

        ptasks = 20
        exec_dir = '/home/mpschr/Documents/projects/exon'
        yaml_file = '/home/mpschr/Documents/projects/exon/config/trios/{}.yaml'.format(patientid)

        if self.do_run_task('bcbio', 0, yaml_file):
            bcbio_ret = bcbio.run(yaml_file, patientid, parallel_tasks=ptasks, exec_dir=exec_dir)

        verification_file = os.path.join(exec_dir, 'out_trios', patientid, 'bcbio', "{}_ID".format(patientid), "_ID-ready.bam".format(patientid))
        sample_dir = '/home/mpschr/Documents/projects/exon/out_trios/{}'.format(patientid)
        if self.do_run_task('cpwr', bcbio_ret, verification_file):
            copywriter_ret = copywriter.run(patient_dir=sample_dir)

        verification_file =  os.path.join(exec_dir, 'out_trios', patientid, 'copywriteR', 'CNAprofiles', 'log2_CNA.segmented.tsvs')
        if self.do_run_task('sciclone', copywriter_ret, verification_file):
            sciclone_ret = sciclone.run(sample_dir)

            if sciclone_ret == 0:
                do_plots(sample_dir)

def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=True)

    options = parser.parse_args()
    task = ExonAnalysisRun()
    task.run(options.patientid)

if __name__ == "__main__":  # detects if called from system
    cmdline()
