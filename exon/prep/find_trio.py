import argparse
import fnmatch
from logging import info,warning,error,debug, basicConfig, DEBUG, INFO
import os
import re
import sys
import pandas as pd

TIME_POINT = 'TIME_POINT'

PATIENT = 'PATIENT'

FASTQ_FILE = 'FASTQ_FILE'

R1 = '_R1'
R2 = '_R2'
CMD = 'cat {} {} > {}'
SAMPLE_PATTERN = "(PL|PE|AL|AE)[0-9]{2}"

__author__ = 'ARJ'

#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)

class TrioFinder:

    def __init__(self):
        self.matches = {}
        self.meta = None

    def run(self, inputpath, patientid, required_time_points=[]):
        info('FINDING TRIOS\n=========================')

        for p_id in patientid:
            if not re.compile("^{}$".format(SAMPLE_PATTERN)).match(p_id):
                raise Exception("Not a valid pattern {}".format(p_id))

        filename_cache = []

        for i in inputpath:
            sample_fastq_pattern = re.compile("^{}.+_L.+.fastq.gz$".format(SAMPLE_PATTERN))
            meta_pattern = re.compile(".*meta.tsv")
            for root, dirnames, filenames in os.walk(i):
                for filename in filenames:
                    if self.meta is None and meta_pattern.match(filename):
                        self._load_meta(os.path.join(root, filename))
                    if sample_fastq_pattern.match(filename):
                        patient, timepoint = os.path.basename(filename).split("_")[0:2]
                        filename_cache.append([filename, patient, timepoint])
                    elif self.meta is not None and filename in self.meta.FASTQ_FILE.tolist():
                        patient, timepoint = self.meta.set_index(FASTQ_FILE).loc[filename].SAMPLE_ID.split("_")
                        filename_cache.append([filename, patient, timepoint])

        if len(filename_cache) == 0:
            raise Exception("no fastq.gz files in path!")
        else:
            filename_cache_df = pd.DataFrame(filename_cache)
            filename_cache_df.columns = [FASTQ_FILE, PATIENT, TIME_POINT]

        for p_id in patientid:

            patient_matches = filename_cache_df.query("{} == '{}'".format(PATIENT, p_id))
            if required_time_points != []:
                patient_matches = patient_matches.query("{} == {}".format(TIME_POINT, required_time_points))
            patient_TPs = patient_matches.TIME_POINT.drop_duplicates()
            if len(patient_TPs) == 3:
                self.matches[p_id] = []
                for TP, group in patient_matches.groupby(TIME_POINT):
                    self.matches[p_id].append([p_id + '_' + TP, group.FASTQ_FILE.tolist()])
                    #self.matches[p_id] = [p_id + '_' + TP for TP in patient_TPs]
            else:
                raise RuntimeError('ERROR: No files or no trio found for id {},{}'.format(p_id, patient_matches))

        for patient, timepoints in self.matches.items():
            info("sample trio found for {}: {}".format(patient, '\n\t\t\t\t\t' + '\n\t\t\t\t\t'.join([(x + ": " + str(len(y)) + " fastqs") for x,y in timepoints])))

    def _load_meta(self, sequencing_meta_file):
        print("loading meta.tsv!")
        self.meta = (
            pd.read_table(sequencing_meta_file)
            .query('SAMPLE_ID != "undetermined"')
            .assign(SAMPLE_ID = lambda df: df.SAMPLE_ID.apply(lambda s: s.replace("_2", "2")))
            .assign(PATIENT = lambda df: df.SAMPLE_ID.apply(lambda s: s.split('_')[0]))
        )



def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='append', dest='inputpath', required=True)
    parser.add_argument('-s', dest='sampleids', action="append", required=True)
    options = parser.parse_args()
    task = TrioFinder()
    task.run(options.inputpath, options.sampleids)

if __name__ == "__main__":  # detects if called from system
    cmdline()
