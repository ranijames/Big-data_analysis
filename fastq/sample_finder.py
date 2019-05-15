import logging
import traceback
from allpy.SimpleTask import SimpleTask
from allpy.fastq.cohorts_regex import FASTQ_COHORTS_REGEX, FastqRegexCohort

__author__ = 'arj'

import argparse
from logging import info
import os
import itertools
import re
import sys

R1 = '_R1'
R2 = '_R2'

class SampleFastQ():
    def __init__(self, sample, patient, time_point, lane, read_pair, filename):
        self.filename = filename
        self.lane = lane
        self.time_point = time_point.upper()
        self.sample = sample.upper()
        self.patient = patient.upper()
        self.read_pair = read_pair

    def __repr__(self):
        return " ".join(["+++\n",self.sample,"\n",self.patient, self.time_point, self.lane, self.read_pair,"\n",self.filename,"\n---"])



class Patient():

    def __init__(self):
        self.patient = None
        self.time_points = {}
        self.samples = []


    def add_fastq(self, sample_fastq: SampleFastQ):
        if self.samples == []:
            self.patient = sample_fastq.patient
        else:
            assert self.patient == sample_fastq.patient

        self.samples.append(sample_fastq)
        idx = len(self.samples) - 1

        if sample_fastq.time_point not in self.time_points:
            self.time_points[sample_fastq.time_point] = []
        self.time_points[sample_fastq.time_point].append(idx)

    def __repr__(self):
        return "Patient: {} -> {} fastqs\n\t\t\t{}".format(self.patient, len(self.samples),
                                                     "\n\t\t\t".join(["{}: {} fastqs".format(tp, len(l)) for tp, l in self.time_points.items()]))

class SampleFastQGenerator():

    def __init__(self, patient_extractor_pattern, time_point_extractor_pattern):
        self.time_point_pattern = time_point_extractor_pattern
        self.patient_pattern = patient_extractor_pattern
        self.lane_pattern = "_(L[0-9]{3})_"

    def extract_sample_fastq(self, filenames, sampleids):
        sample_fastqs = []
        for f in filenames:
            patient = re.search(self.patient_pattern, f, re.IGNORECASE).group(1)
            time_point = re.search(self.time_point_pattern, f, re.IGNORECASE).group(1)
            lane = re.search(self.lane_pattern, f, re.IGNORECASE).group(1)
            sample = "{}_{}".format(patient, time_point)
            read_pair = R1.replace('_', '') if R1 in f else R2.replace('_', '')

            if sample in sampleids or "*" in sampleids:
                sfastq = SampleFastQ(sample, patient, time_point, lane, read_pair, f)
                sample_fastqs.append(sfastq)
        return sample_fastqs

class SampleFinder(SimpleTask):

    def __init__(self):
        super().__init__("SAMPLE FINDER")
        self.patients = {}

    def run(self, inputpath, cohort_key, sampleids):
        super().run()

        logging.debug(cohort_key)
        sample_regex = FASTQ_COHORTS_REGEX[cohort_key]
        assert isinstance(sample_regex, FastqRegexCohort)

        if sampleids is None:
            sampleids = ["*"]

        filename_cache = []
        for i in inputpath:
            i = os.path.abspath(i)
            sample_fastq_pattern = re.compile(sample_regex.get_fastq_pattern(), re.IGNORECASE)
            for root, dirnames, filenames in os.walk(i):
                logging.debug("Scanning {}".format(root))
                for filename in filenames:
                    if sample_fastq_pattern.match(filename):
                        filename = os.path.join( root, filename)
                        filename_cache.append(filename )
        if len(filename_cache) == 0:
            raise Exception("no fastq.gz files in path!")

        sample_generator = SampleFastQGenerator(sample_regex.get_patient_extractor(), sample_regex.get_time_point_extractor())
        sample_fastqs = sample_generator.extract_sample_fastq(filename_cache, sampleids)

        key = lambda x: x.patient
        patients = 0
        samples = 0
        fastqs = 0
        for key, g in itertools.groupby(sorted(sample_fastqs, key=key), key=key):
            p = Patient()
            for s in g:
                p.add_fastq(s)
            patients += 1
            samples += len(p.time_points)
            fastqs += len(p.samples)
            self.patients[p.patient] = p
        logging.info("Found {} fastqs for {} samples of {} patients".format(fastqs, samples, patients))



        # for s in sampleids:
        #     sample_matches = []
        #     for filename in filename_cache:
        #         patient, timepoint = os.path.basename(filename).split("_")[0:2]
        #         if timepoint not in sample_matches:
        #             sample_matches.append(timepoint)
        #
        #     if len(sample_matches) == 3:
        #         self.matches[s] = [s + '_' + m for m in sample_matches]
        #
        #     else:
        #             raise RuntimeError('ERROR: No files or no trio found for id {},{}'.format(s, sample_matches))
        #             #warning('ERROR: No files or no trio found for id {},{}'.format(s, sample_matches))
        #
        # for sample, timepoints in self.matches.items():
        #     info("sample trio found {}".format(timepoints, sample))



def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='append', dest='input_fastq_path', required=True)
    parser.add_argument('-C', dest='cohort_key', choices=list(FASTQ_COHORTS_REGEX.keys()), required=True)
    parser.add_argument('-s', dest='sampleids', action="append", required=False)
    options = parser.parse_args()
    task = SampleFinder()
    task.run(options.input_fastq_path, options.cohort_key, options.sampleids)

if __name__ == "__main__":  # detects if called from system
    cmdline()
