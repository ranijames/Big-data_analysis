import logging
from os.path import join
import pandas as pd
import subprocess
from io import BytesIO

__author__ = 'ARJ'
logger = logging.getLogger()


CMD_GEMINI_VARIANTS = '''gemini query --header -q \
             "SELECT chrom, start, end, gene, exon, ref, alt, aa_change, rs_ids, cosmic_ids, \
                aaf_adj_exac_all, type, sub_type, impact_so
               FROM variants " \
               {}'''


# select the sample which is not 0!
CMD_GEMINI_SAMPLES = '''gemini query -q \
            "SELECT name \
            FROM samples \
            WHERE instr(name, \'CR\') == 0" \
{}'''

class GeminiExtractor:

    def __init__(self, muts_dir, sample, calling='ensemble'):
        self.sample = sample
        self.sample_muts_file = join(muts_dir, 'b{}-{}.db'.format(sample, calling))

    @staticmethod
    def remove_multispaces(somestring):
        return " ".join(somestring.split())

    @staticmethod
    def calc_alt_allele_frequency(row, sample_name):
        alt_depth, total_depth = row[['gt_alt_depths.' + sample_name, 'gt_depths.' + sample_name]]
        if int(alt_depth) < 0 or int(total_depth) < 0:
            return 0

        return int(alt_depth)/int(total_depth)

    def get_df(self, min_alt_reads=1, min_reads=20, coding_only=True, debug_run=False) -> pd.DataFrame:

        logging.info('reading ' + self.sample_muts_file)

        gemini_variants_cmd = self.remove_multispaces(
            CMD_GEMINI_VARIANTS.format(self.sample_muts_file)
        )
        logger.debug(gemini_variants_cmd)

        p2 = subprocess.Popen(gemini_variants_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out2, err = p2.communicate()

        #create df
        df = pd.read_csv(BytesIO(out2), sep='\t')


        if 'start' in df.columns and 'end' in df.columns:
            df.insert(1, 'pos', df['start'] + 1)
            df.drop(['start', 'end'], axis=1, inplace=True)

        return df



def relapsed_muts(ID, REL):
    on_cols = ['chrom', 'start', 'end', 'ref', 'alt', 'gene', 'sub_type']
    cols_of_interest = on_cols + ['gt_alt_depths', 'alt_freq']
    df = pd.merge(ID[cols_of_interest], REL[cols_of_interest], on=on_cols, how='inner', suffixes=['_ID','_REL'])
    df['delta_alt_freq'] = df.alt_freq_REL - df.alt_freq_ID
    return df
