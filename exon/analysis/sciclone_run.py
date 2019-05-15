import sys
import argparse
import logging
import os
import subprocess
from os.path import join

from allpy.exon.mutations.mutreport_run import MutationReporter
from allpy.exon.mutations.utils import write_df, Timer
from allpy.notify.notify import send_email

MIN_FREQ = 0

MIN_READS = 40
MIN_ALT_READS = 5


CMD_OUTER = "{0} ; notify-send \"Analysis DONE\" \"sciclone done for sample {1}\" -t 3000  "

CMD_SCICLONE = "Rscript /home/mpschr/Documents/projects/exon/exonR/scicloneRun.R -p {0} -o {1} -i {2} -c {3} -m {4} -f {5} --highlighted.genes /home/mpschr/Documents/projects/exon/exonR/highlightedgenes.txt"

logger = logging.getLogger()

def create_sciclone_muts_file(patient, reporter: MutationReporter, out_dir, calling='ensemble', min_alt_reads=10, min_reads=40, min_frequency=0.05):

    time_points = ['ID', 'REL']
    sample_muts = reporter.get_df(time_points=time_points, min_alt_reads=min_alt_reads, min_reads=min_reads, min_frequency=min_frequency, time_point_normal='CR')
    filter_metadata = reporter.get_metadata()
    #muts_dir = reporter.get_mut_folder()


    read_cols = []
    for t in time_points:

        reads_col = 'reads_{}'.format(t)
        alt_reads_col = 'alt_reads_{}'.format(t)
        ref_reads_col = 'ref_reads_{}'.format(t)
        alt_freq_col = 'alt_freq_{}'.format(t)

        read_cols += reads_col, alt_reads_col, alt_freq_col

        sample_muts[ref_reads_col] = sample_muts[reads_col] - sample_muts[alt_reads_col]
        sample_muts[alt_freq_col] = sample_muts[alt_freq_col] * 100

    #sample_muts = sample_muts.reset_index()
    #reduced_df = sample_muts[['chrom', 'end', reads_col, alt_reads_col, alt_freq_col, 'gene', 'base_qual', 'mapping_qual']]

    reduced_df = sample_muts[['chrom', 'pos'] + read_cols + [ 'genes', 'pam_genes_count', 'known_cancer_gene']].rename(columns  = {'genes': 'gene'})
    reduced_df.gene = reduced_df.gene.fillna('-')

    print(reduced_df.head())


    filterstring =  'depth{}-alt{}-freq{}'.format(min_reads, min_alt_reads, min_frequency)
    out_file_name = join(out_dir,
                         'b{}-{}.sciclone-muts.{}.tsv'.format(patient, calling, filterstring))
    write_df(reduced_df, out_file_name, filter_metadata, index=False)
    #reduced_df.to_csv(out_file_name, sep="\t", index=False)

    logging.info("Written into: " + out_file_name)

    return out_file_name, filterstring


class SciCloneRunner():

    def __init__(self):
        pass

    def run(self, sample_dir, patientid=None, min_reads=MIN_READS,
            min_alt_reads=MIN_ALT_READS, min_frequency=MIN_FREQ, allmuts=True):

        timer = Timer(start_now=True)
        out_dir = os.path.join(sample_dir, 'sciClone', 'allmuts' if allmuts else 'coding')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        if patientid is None:
            patientid = os.path.basename(sample_dir)

        #logging:
        log_file = os.path.join(out_dir, 'sciclone-run.{}'.format(min_reads))
        logging.basicConfig(filename=log_file+'.pylog', filemode="w+", format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)
        logging.debug(log_file)

        reporter = MutationReporter(patientid, sample_dir)
        muts_file, filterstring = create_sciclone_muts_file(reporter=reporter, out_dir=out_dir, patient=patientid,
                                                 min_reads=min_reads, min_alt_reads=min_alt_reads, min_frequency=min_frequency)

        cna_segments_file = os.path.join(sample_dir, 'copywriteR', 'CNAprofiles', 'log2_CNA.segmented.tsv')

        cmd = CMD_SCICLONE.format(patientid, out_dir, muts_file, cna_segments_file, min_reads, filterstring)
        logging.info('Full command: "{}"'.format(cmd) )
        rlogfile = open(log_file+'{}.Rlog'.format(filterstring), 'w+')
        c = subprocess.Popen(cmd, shell=True, stderr=rlogfile, stdout=rlogfile)
        streamdata = c.communicate()[0]
        rlogfile.flush()
        rlogfile.close()
        logging.debug('R Return code {}'.format(c.returncode))

        elapsed_time = timer.get_elapsed_time()
        logging.info(elapsed_time) 

        if c.returncode == 0:
            send_email('michael.p.schroeder@gmail.com',"{} SCICLONE done".format(patientid), elapsed_time)
        else:
            send_email('michael.p.schroeder@gmail.com',"FAILURE of {} SCICLONE".format(patientid), elapsed_time)
        return c.returncode


def cmdline():
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-D', dest='sample_dir', required=True, help='Where the output folder from bcbio is, the root '
                                                                       'for the sample results')
    parser.add_argument('-m', dest='min_reads', default=MIN_READS, help="Minimum depths for variants to be considered")
    parser.add_argument('-a', dest='min_alt_reads', default=MIN_ALT_READS, help="Minimum depths for variants to be considered")
    parser.add_argument('-f', dest='min_frequency', default=MIN_FREQ, help="Minimum variant frequency (VAF) to be considered")

    options = parser.parse_args()
    task = SciCloneRunner()
    task.run(os.path.abspath(options.sample_dir), options.patientid, options.min_reads, options.min_alt_reads,
             options.min_frequency)

if __name__ == "__main__":  # detects if called from system
    cmdline()
