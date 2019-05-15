import argparse
from logging import info, debug
from allpy.SimpleTask import SimpleTask
from allpy.fastq.cohorts_regex import FASTQ_COHORTS_REGEX
from allpy.fastq.sample_finder import SampleFinder, Patient, SampleFastQ
from os.path import join, exists, abspath
import os
import pandas as pd

__author__ = 'ARJ'

class SampleLinker(SimpleTask):

    def __init__(self):
        super().__init__('SAMPLE LINKER')
        self.meta = None

    def run(self, sample_finder: SampleFinder, output_directory: str):
        super().run()

        for p in sample_finder.patients.values():
            assert isinstance(p, Patient)
            info(p)
            patient_dir = join(abspath(output_directory), p.patient)
            if not exists(patient_dir):
                os.makedirs(patient_dir)

            for sample in p.samples:
                assert isinstance(sample, SampleFastQ)
                sample_dir = join(patient_dir, sample.sample)
                if not exists(sample_dir):
                    os.makedirs(sample_dir)

                f = sample.filename
                ext = 'fastq.gz' if f.endswith('fastq.gz') else 'fastq'
                linkpath = os.path.join(sample_dir, "{}_{}.{}".format(sample.sample, sample.read_pair, ext))
                if not (os.path.lexists(linkpath) or os.path.exists(linkpath)):
                    os.symlink(src=f, dst=linkpath)
                    info('LINKED: {} as {}'.format(os.path.basename(f), os.path.basename(linkpath)))
                elif os.path.exists(linkpath) and os.readlink(linkpath) == f:
                    debug('already linked: {}'.format(os.path.basename(f)))
                    continue
                elif os.path.lexists(linkpath) and os.readlink(linkpath) != f:
                    raise Exception('Link exists, seems to point to a bad destination: {} ---> {}'.format(linkpath, os.path.realpath(linkpath)))

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
    parser.add_argument('-i', action='append', dest='input_fastq_path', required=True)
    parser.add_argument('-o', dest='output_directory', required=True)
    parser.add_argument('-C', dest='cohort_key', choices=list(FASTQ_COHORTS_REGEX.keys()), required=True)
    parser.add_argument('-s', dest='sampleids', action="append", required=False)
    parser.add_argument('-m', dest='metafile')
    options = parser.parse_args()

    finder = SampleFinder()
    finder.run(options.input_fastq_path, options.cohort_key, options.sampleids, options.metafile)

    task = SampleLinker()
    task.run(finder, options.output_directory)
