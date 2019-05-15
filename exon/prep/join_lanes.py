import argparse
import fnmatch
from logging import basicConfig, DEBUG, info
import subprocess
import os
import tempfile

#logging:
import sys
import pandas as pd

basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


R1 = '_R1'
R2 = '_R2'
CMD = 'cat {1} > {0}; mv {0} {2}'
CMD_LN = 'ln {0} {1}'

__author__ = 'ARJ'


class LaneJoiner():

    def __init__(self):
        self.joined_lanes = {}
        self.previously_joined_lanes = {}

    def join_lanes(self, sample, files, test_run):
        lanes_r1 = [f for f in files if R1 in f]
        lanes_r2 = [f for f in files if R2 in f]

        summary_dict = {'files': files, 'output': []}

        if len(lanes_r1) > 1:
            self._join_lanes_exec(lanes_r1, sample, summary_dict, test_run, R1)
        else:
            self._link_file(lanes_r1, sample, summary_dict, test_run, R1)

        if len(lanes_r2) > 1:
            self._join_lanes_exec(lanes_r2, sample, summary_dict, test_run, R2)
        else:
            self._link_file(lanes_r2, sample, summary_dict, test_run, R2)

        if len(summary_dict['output']) > 0:
            self.joined_lanes[sample] = summary_dict

    def _link_file(self, lanes, sample, summary_dict, test_run, r_suffix):

        outfile = joined_file_name(lanes[0], sample, r_suffix)

        cmd = CMD_LN.format(
            lanes[0], outfile
        )
        info('COPYING 1LANE-SAMPLE {} for {}: {}'.format(r_suffix, sample, cmd))
        if test_run:
            info('TEST RUN! COMMAND NOT SENT TO OS')
        else:
            c = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).wait()
            summary_dict['output'].append(outfile)



    def _join_lanes_exec(self, lanes, sample, summary_dict, test_run, r_suffix):

        #ramdisk = '/media/ramdisk/'
        ramdisk = '.'
        outfile = joined_file_name(lanes[0], sample, r_suffix)
        tmpfile = tempfile.mkstemp(dir=ramdisk, prefix=os.path.basename(outfile) + '_')[1]

        cmd = CMD.format(
            tmpfile, " ".join(lanes), outfile
        )

        info('JOINING {} for {}: {}'.format(r_suffix, sample, cmd))
        if test_run:
            info('TEST RUN! COMMAND NOT SENT TO OS')
        else:
            c = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).wait()
            summary_dict['output'].append(outfile)

    def run(self, inputpath, sampleids, sample_files=None, test_run=False):
        info('FINDING FILES TO JOIN\n===================================')

        filename_cache = []
        for i in inputpath:
            for root, dirnames, filenames in os.walk(i):
                #if os.path.basename(root) == 'raw_sequence_data':
                filename_cache += [os.path.join(root, f) for f in filenames if f.endswith('.fastq.gz')]

        previously_done = []
        matches = {}
        for s_id, s_files in zip(sampleids, sample_files):
            sample_matches = []
            print(s_files)
            if len(s_files) > 0:
                series = pd.Series(filename_cache)
                series.index = series.apply(lambda s: s.split('/')[-1])
                sample_matches = series[series.index.isin(s_files)].tolist()
            else:
                for filename in fnmatch.filter(filename_cache, '*{}*_L*.fastq.gz'.format(s_id)):
                    sample_matches.append(filename)

            if len(sample_matches) > 0:

                sample_matches = self.remove_existing(s_id, sample_matches)

                if len(sample_matches) == 0:
                    previously_done.append(s_id)

                elif len(sample_matches) % 2 == 0:
                    matches[s_id] = sample_matches
                    info('FOUND {} files for {}{}'.format(len(sample_matches), s_id, R1 if R1 in sample_matches[0] else R2))

                else:
                    raise RuntimeError('ERROR: Wrong number of files found for sample with id {}: {} files\n\t{}'.format(
                        s_id, len(sample_matches), '\n\t'.join(sample_matches)
                    ))
            else:
                    print(sample_matches)
                    raise RuntimeError('ERROR: No files found for id {}'.format(s_id))

        info('JOINING FILES\n=========================')
        for sample, files in matches.items():
            self.join_lanes(sample, files, test_run)

    def remove_existing(self, sample, fastq_files):

        summary = {'output': []}

        r1_out = joined_file_name(fastq_files[0], sample, R1)
        if os.path.exists(r1_out):
            fastq_files = [f for f in fastq_files if R1 not in f]
            info('SKIPPING: R1 for sample {} already joined'.format(sample))
            summary['output'].append(r1_out)

        r2_out = joined_file_name(fastq_files[0], sample, R2)
        if os.path.exists(r2_out):
            fastq_files = [f for f in fastq_files if R2 not in f]
            info('SKIPPING: R2 for sample {} already joined'.format(sample))
            summary['output'].append(r2_out)

        self.previously_joined_lanes[sample] = summary

        return fastq_files

def joined_file_name(orig_file, sample, readgroup):
        return os.path.join(os.path.dirname(orig_file), sample + readgroup + '.fastq.gz')

def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='append', dest='inputpath', required=True)
    parser.add_argument('-s', dest='sampleids', action="append", required=True)
    parser.add_argument('-t', '--test-run', action='store_true', dest='test',  default=False)
    options = parser.parse_args()
    task = LaneJoiner()
    task.run(options.inputpath, options.sampleids, test_run=options.test)

if __name__ == "__main__":  # detects if called from system
    cmdline()
