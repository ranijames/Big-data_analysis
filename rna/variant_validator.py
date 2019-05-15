import argparse
from logging import basicConfig, DEBUG, info
import logging
import os
import sys
import vcf
import pandas
from allpy.SimpleTask import SampleTask

from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL, CORES, \
    GENOME, SAMTOOLS_BIN, BCFTOOLS_BIN, SAMBAMBA_BIN
from allpy.exon.mutations.mutreport_run import annotate_protein_coding, ALL_MUTS
from allpy.exon.mutations.utils import write_df
from allpy.externalcommands import ExternalCommand, ExternalCommandOutput
from allpy.notify.notify import send_email

CMD_WHATEVER = """{CUFFLINKS_BIN} \
 --num-threads {CORES} \
 --quiet \
 --no-update-check \
 --max-bundle-frags 2000000 \
 --library-type fr-unstranded \
 --frag-bias-correct {GENOME} \
 --multi-read-correct \
 --GTF{GUIDE_OR_NO} {TRANSCRIPT_GTF} \
 --mask-file {TRANSCRIPT_GTF_MASK} \
 --output-dir {WORK_DIR} \
 {BAM_FILE}
"""
__author__ = 'ARJ'

#logging:
basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


def reads(sample_call: vcf.model._Call, name):
    d = {}
    d['reads_{}'.format(name)] = sample_call.data.DP
    d['alt_reads_{}'.format(name)] = sample_call.data.DV
    try:
        d['alt_freq_{}'.format(name)] = sample_call.data.DV / sample_call.data.DP
    except ZeroDivisionError:
        d['alt_freq_{}'.format(name)] = 0.0
    return d

class VariantValidator(SampleTask):

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            sample = cmd_sample_name
        else:
            sample = os.path.basename(os.path.abspath(task_input[0][0])).split('.bam')[0]
        return sample

    def add_metadata(self, metastring, do_log=True):
        self.metadata.append(metastring)
        if do_log:
            logging.info(metastring)

    def __init__(self):
        super().__init__('EXON_RNA_VARIANT_VALIDATOR')
        self.metadata = []

    def run(self, input_bams, input_variants, exec_dir='.', sample_name=None, detect_unannotated=False, scratch_dir=None):
        super().run(input_bams, exec_dir, sample_name)

        os.chdir(exec_dir)
        if not os.path.exists('mutations'):
            os.mkdir('mutations')

        conf = self.load_configuration()

        if scratch_dir is None:
            scratch_dir = os.path.join(exec_dir, 'mutations')

        samples = []
        # get sample name EXON
        sample = self.bam_samplename(input_bams[0][0], conf[SAMTOOLS_BIN])
        input_bams[0].append(sample)
        samples.append(sample)

        # get sampmle name RNA
        sample = self.bam_samplename(input_bams[1][0], conf[SAMTOOLS_BIN])
        input_bams[1].append(sample)
        samples.append(sample)

        ## output files
        jointcall_vcf = 'mutations/{}-{}.vcf'.format(self.sample, ALL_MUTS)
        jointcall_tsv = 'mutations/{}-{}.tsv'.format(self.sample, ALL_MUTS)

        bams = []
        if len(set(samples)) < len(input_bams):
            logging.warning("reapeated sample names! Will be renamed")
            # create new bams if necessary
            counter = 0
            for bam, name, samplename in input_bams:
                counter += 1
                new_samplename = self.sample + '_' + name
                samples[counter-1] = new_samplename

                # avoid run if .vcf already exists
                if os.path.exists(jointcall_vcf):
                    info('SKIPPING bam reheader: .vcf already exists')
                    continue
                ExternalCommand(
                    cmd='{5} view -H {0} | sed "s/SM:{1}/SM:{2}/" | samtools reheader - {0} > {7}/tmp{3}-{4} ; {6} index {7}/tmp{3}-{4}'
                                .format(bam,
                                        samplename,
                                        new_samplename,
                                        counter,
                                        os.path.basename(bam),
                                        conf[SAMTOOLS_BIN],
                                        conf[SAMBAMBA_BIN],
                                        scratch_dir),
                    output_file="{}/tmp{}-{}".format(scratch_dir, counter, os.path.basename(bam))
                ).run()
                bams.append("{}/tmp{}-{}".format(scratch_dir, counter, os.path.basename(bam)))
        else:
            for bam, name, samplename in input_bams:
                bams.append(bam)

        logging.info("Samples: {}".format(', '.join(samples)))

        # create bed regions file
        variants_bed = 'mutations/{}-variant-sites.bed'.format(self.sample)
        known_variants = pandas.read_table(input_variants, sep='\t', comment='#')
        known_variants['chrom'] = known_variants.chrom.apply(str)
        self.add_metadata('Loaded {} variants from {}'.format(known_variants.shape[0],
                                                                                  input_variants))
        known_variants['pos'] = known_variants['start'] + 1
        known_variants[['chrom', 'pos']].to_csv(variants_bed, sep='\t', header=False, index=False)
        annotate_protein_coding(known_variants, 'impact_so', 'is_coding')

        # mpileup and bcftools call
        cmd = "{SAMTOOLS_BIN} mpileup -u -Q 0 -t DP4 -t DV -t DP -f {GENOME} --positions {VARIANTS} {BAMS} " \
        "| {BCFTOOLS_BIN} view | bcftools norm -m-both -f {GENOME} - | grep -v '<X>' > {OUTPUTFILE}".format(**{'SAMTOOLS_BIN': conf[SAMTOOLS_BIN],
                                                            'GENOME': conf[GENOME],
                                                            'VARIANTS': variants_bed,
                                                            'BAMS': " ".join(bams),
                                                            'BCFTOOLS_BIN': conf[BCFTOOLS_BIN],
                                                            'OUTPUTFILE': jointcall_vcf})
        c = ExternalCommand(cmd, output_file=jointcall_vcf).run()
        ### samtools mpileup -u -t DP4 -t DV -t DP -f {GENOME} --positions exon_variants.bed  {BAMS} | bcftools call -m  > jointcall.vcf

        #cmd=''
        #c = ExternalCommand(cmd).run()

        reader = vcf.Reader(open(jointcall_vcf))
        sites = []
        for v in reader:
            exon = v.genotype(samples[0])
            rna = v.genotype(samples[1])

            site = {}
            site['chrom'] = v.CHROM
            site['pos'] = v.POS
            site['ref'] = v.REF
            site['alt'] = str(v.ALT[0])

            if len(v.ALT) > 1:
                raise RuntimeError('More than one alternative found.')

            site.update(reads(exon, "_".join(exon.sample.split('_')[-2:])))
            site.update(reads(rna, "_".join(rna.sample.split('_')[-2:])))

            sites.append(site)

        sitesdf = pandas.DataFrame.from_dict(sites)
        self.add_metadata("Read {} variants at {} sites from {}".format(sitesdf.shape[0], len(sitesdf.pos.unique()), jointcall_vcf))

        #sitesdf.sort("alt_freq_ID_RNA", ascending=False).head()

        variant_idx = ['chrom', 'pos', 'alt', 'ref']
        joined = known_variants.set_index(variant_idx).join(sitesdf.set_index(variant_idx), how='left')

        field = "alt_reads_" + "_".join(samples[1].split("_")[1:])
        count = joined.query(field + ' > 1').shape[0]
        protein_coding_count = joined.query( field + ' > 1 and is_coding == "is_coding" ').shape[0]

        self.add_metadata('Confirmed {} variant sites ({} protein affecting) with at least 2 alternative read in RNA'.format(count, protein_coding_count))

        self.add_metadata('Writing {} variants to tsv'.format(joined.shape[0]))
        write_df(joined, jointcall_tsv, self.metadata, index=True)

        # delete newly created bams if created
        ExternalCommand(cmd='rm {0}/tmp*.bam {0}/tmp*.bai'.format(scratch_dir)).run()

        #return c
        elapsed = self.timer.get_elapsed_time()
        logging.info(elapsed)

    def bam_samplename(self, bam, samtools_bin):
        sample = None
        assert os.path.exists(bam), "{} cannot be found".format(bam)
        c = ExternalCommandOutput("{} view -H {}".format(samtools_bin, bam))
        out = c.run()
        for line in out.split("\n"):
            if line.startswith('@RG') and 'SM:' in line:
                sample = line.split('SM:')[1]
        return sample

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL,
                                            CORES,
                                            GENOME,
                                            SAMTOOLS_BIN,
                                            BCFTOOLS_BIN,
                                            SAMBAMBA_BIN
                                            ])
        return conf_dict


def cmdline():
    usage_example = '''seq-variant-validator -i exonseq.bam EXON -i rnaseq.bam RNA  -v variantlist.vcf -s SampleX '''
    parser = argparse.ArgumentParser(description='''Validate the variants found in e.g. exon''',
                                     usage=usage_example)

    parser.add_argument('-i', dest='input_bams', required=True, nargs=2, metavar=('BAM-FILE', 'NAME'), action='append',
                        help='The bam file(s) and and their names to be used. Names help to distinguish samples if two bam files have'
                                                                    'the same sample name. Minimum two')

    parser.add_argument('-v', dest='variants_file', required=True, help='The file to retrieve the variants from for which to check. Mutations MASTER file format.')

    parser.add_argument('-x', dest='exec_dir', default='.', required=False, help='Where the program should be run and results stored')
    parser.add_argument('-s', dest='sample_name', required=False, help='Only needed to specify if the samplename should be '
                                                                       'different than the bam file name')
    parser.add_argument('--scratch', default=None, required=False, help='A scratch dir to save a temporary files. Recommended to use'
                                                                        'a local storage location, if the output dir is remote')

    options = parser.parse_args()

    if len(options.input_bams) < 2:
        raise RuntimeError('At least two bam files and names should be specified! As in  \n\n{}'.format(usage_example))

    task = VariantValidator()
    task.run(options.input_bams, options.variants_file, options.exec_dir, options.sample_name, scratch_dir=options.scratch)

if __name__ == "__main__":  # detects if called from system
    cmdline()
