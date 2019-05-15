import gzip
import math
import os
import re
import time
import pandas as pd
from os.path import join
from itertools import takewhile


__author__ = 'ARJ'

def get_mut_batch_folder(patientid, sample_dir, bcbio_backup):
    filter_string = "" + "_.+[0-9]{6}$"
    bcbio_dir = join(sample_dir,'bcbio')
    muts_dir = None
    try:
        muts_dir = [match for match in filter(lambda x: re.search(filter_string, x), os.listdir(bcbio_dir))][0]
    except (FileNotFoundError, IndexError):
        bcbio_dir = join(bcbio_backup, os.path.basename(sample_dir), 'bcbio')
        muts_dir = [match for match in filter(lambda x: re.search(filter_string, x), os.listdir(bcbio_dir))][0]
    return join(bcbio_dir, muts_dir)


def flatten_dict(d):
    def expand(key, value):
        if isinstance(value, dict):
            return [ (key + '.' + k, v) for k, v in flatten_dict(value).items() ]
        else:
            return [ (key, value) ]

    items = [ item for k, v in d.items() for item in expand(k, v) ]

    return dict(items)


def nulls_in_dict(d):
    for x in list(d.values()):
        if pd.isnull(x):
            return True
    return False


def chromosome_order(chrom_series):
    return sorted(chrom_series.unique().tolist(), key=lambda x: (-x.isnumeric(), len(x), x))

def write_df(df, out_file_name, metadata=[], index=True):
    out_file = open(out_file_name, 'w')
    if len(metadata) > 0:
        for md in metadata:
            out_file.writelines(["# {}\n".format(md)])
    df.to_csv(out_file, sep="\t", index=index)


def read_df(filename, comment='#', **kwargs):
    with open(filename, 'r') as fobj:
        # takewhile returns an iterator over all the lines
        # that start with the comment string
        headiter = takewhile(lambda s: s.startswith(comment), fobj)
        # you may want to process the headers differently,
        # but here we just convert it to a list
        header = list(headiter)
    df = pd.read_table(filename, comment=comment, **kwargs)
    return header,df

class Timer:
    def __init__(self, start_now=True):
        if start_now:
            self.start = time.time()
        else:
            self.start = -1

    def start(self):
        self.start = time.time()

    def get_elapsed_time(self):
        if self.start < 0:
            return 'Not started yet'

        elapsed = time.time() - self.start
        hour_factor = (60 * 60)
        hours = math.floor(elapsed / hour_factor)
        elapsed -= hours * hour_factor
        mins = math.floor(elapsed / 60)
        secs = math.ceil(elapsed % 60)
        elapsed_time = 'Elapsed time: {}h{}m{}s'.format(hours, mins, secs)
        return elapsed_time


def get_commented_rows(f, commentstr="##"):
    skiprows = 0

    fin = None
    if f.endswith('.gz'):
        fin = gzip.open(f, 'r')
    else:
        fin = open(f, 'r')

    for l in fin:
        if type(l) != str:
            l = l.decode()

        if l.startswith(commentstr):
            skiprows += 1
        else:
            break
    return skiprows


def pileup_read_details(row ,samples):
    try:
        DPidx = row.FORMAT.split(':').index('DP')
        ADidx = row.FORMAT.split(':').index('AD')
    except AttributeError as e:
        print(row)
        raise e

    for s in samples:
        sample_details = row[s].split(":")
        row['reads_{}'.format(s)] = int(sample_details[DPidx])
        row['alt_reads_{}'.format(s)] = int(sample_details[ADidx].split(',')[1])
        try:
            row['alt_freq_{}'.format(s)] = int(sample_details[ADidx].split(',')[1]) / int(sample_details[DPidx])
        except ZeroDivisionError:
            row['alt_freq_{}'.format(s)] = 0.0
    return row
