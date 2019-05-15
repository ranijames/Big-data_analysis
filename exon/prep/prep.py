import argparse
import subprocess
from logging import info,warning,error,debug, basicConfig, DEBUG, INFO
import os
import sys
from allpy.exon.prep.find_trio import TrioFinder
from allpy.exon.prep.join_lanes import LaneJoiner
from allpy.exon.prep.create_config import ConfigCreator

__author__ = 'ARJ'
CMD_NOTIFY = 'selfnotify \"{} prepared\"'


#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


def extract_outputs(lane_joiner: LaneJoiner):
    outfiles = []
    for v in lane_joiner.joined_lanes.values():
        outfiles += (v['output'])
    for v in lane_joiner.previously_joined_lanes.values():
            outfiles += (v['output'])
    return outfiles

class SamplesPrep:

    def __init__(self, link_folder, yaml_file, config_dir, test_run=False, time_points=[]):
        self.config_dir = config_dir
        self.yaml_file = yaml_file
        self.test_run = test_run
        self.link_folder = link_folder
        self.time_points = time_points

    def run(self, inputpath, patientids):
        self.find_candidate_trios(inputpath, patientids, self.time_points)

    def find_candidate_trios(self, inputpath, patientids, time_points):
        trio_finder = TrioFinder()
        trio_finder.run(inputpath, patientids, time_points)
        if len(trio_finder.matches) > 0:
            self.join_lanes(trio_finder.matches, inputpath)

    def join_lanes(self, matches, inputpath):
        lane_joiner = LaneJoiner()
        for patient in matches:
            samples, sample_files = zip(*matches[patient])
            lane_joiner.run(inputpath, sampleids=samples, sample_files=sample_files, test_run=self.test_run)
            self.create_datalinks(patient, extract_outputs(lane_joiner))

    def create_datalinks(self, patient, sample_files):
        linkfolder = os.path.join(self.link_folder, patient)
        if os.path.exists(linkfolder):
            warning("Datalink folder for {} already exists: {}".format(patient, linkfolder))
        else:
            os.mkdir(linkfolder)

        for f in sample_files:
            linkpath = os.path.join(linkfolder, os.path.basename(f))
            if not (os.path.lexists(linkpath) or os.path.exists(linkpath)):
                os.symlink(src=f, dst=linkpath)
                info('LINKED: {}'.format(os.path.basename(f)))
            elif os.path.exists(linkpath) and os.readlink(linkpath) == f:
                debug('already linked: {}'.format(os.path.basename(f)))
                continue
            elif os.path.lexists(linkpath) and os.readlink(linkpath) != f:
                raise Exception('Link exists, seems to point to a bad destination: {} ---> {}'.format(linkpath, os.path.realpath(linkpath)))

        self.create_config(linkfolder, patient)

    def create_config(self, linkfolder, patient):
        config_creator = ConfigCreator(self.yaml_file, linkfolder, config_dir=self.config_dir, time_points = self.time_points)
        config_creator.run()

        c = subprocess.Popen(CMD_NOTIFY.format(patient), shell=True, stdout=subprocess.PIPE)
        streamdata = c.communicate()[0]
        debug("return code {}:".format(c))



#subprocess.call("screen -dmS test ./server", shell=True)

def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='append', dest='inputpaths', required=True)
    parser.add_argument('-p', dest='patientids', action="append", required=True)
    parser.add_argument('-t', dest='test_run', action='store_true',  default=False)
    parser.add_argument('-d', dest='link_folder', required=True)
    parser.add_argument('-y', dest='yaml_file', required=True)
    parser.add_argument('-c', '--config-dir', dest='config_dir', required=True, help='where the config file will be written to')
    parser.add_argument('-T', dest='time_points', action="append", default=[])
    options = parser.parse_args()
    task = SamplesPrep(options.link_folder, options.yaml_file, options.config_dir, test_run=options.test_run,
                       time_points = options.time_points)
    task.run(options.inputpaths, options.patientids)

if __name__ == "__main__":  # detects if called from system
    info('Command line')
    cmdline()
else:
    info('Module import')
