import argparse
from logging import basicConfig, DEBUG, info, debug
import logging
import os
import sys
import subprocess
from allpy.exon.mutations.utils import Timer

CMD = "screen -S {1} -dm sh -c 'bcbio_nextgen.py {0} --workdir scratch/scratch_{1} -t local -n {2};" \
      "selfnotify \"{1} BCBIO done\" ;" \
      "echo $?; exec bash -i'"

      # after selfnotify : " mv scratch/scratch_{1}/log out_trios/{1}/bcbio/ ;" \


__author__ = 'ARJ'


#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)

class BcbioRunner():

    def __init__(self):
        pass

    def run(self, yaml_file, patientid, parallel_tasks, exec_dir):

        timer = Timer(start_now=True)

        os.chdir(exec_dir)
        info('working in dir: {}'.format(os.path.abspath(exec_dir)))
        info('RUNNING BCBIO IN SCREEN AND EXIT\n=========================')
        cmd = CMD.format(yaml_file, patientid, parallel_tasks)
        debug(cmd)
        c = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        logging.debug('Return code {}'.format(c.returncode))

        return c.returncode

#bcbio_nextgen.py config/trios/AL02.yaml --workdir scratch_AL02 -t local -n 20
#screen -S testy -dm sh -c 'ls -al; echo $?; exec bash -i'

def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=True)
    parser.add_argument('-x', dest='exec_dir', required=True, help='Where the program should be run and results stored',
                        default='/home/mpschr/Documents/projects/exon')
    parser.add_argument('-y', dest='yaml_file', required=True)
    parser.add_argument('-n', dest='parallel_tasks', type=int, default=20, help="default 20")
    options = parser.parse_args()
    task = BcbioRunner()
    task.run(options.yaml_file, options.patientid, options.parallel_tasks, options.exec_dir)

if __name__ == "__main__":  # detects if called from system
    cmdline()
