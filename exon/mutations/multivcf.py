import logging
import os
import time
import math
from os.path import join as pathjoin
from collections import Iterator

import pandas as pd
import vcf
import vcf.utils
import numpy as np

from allpy.exon.mutations.utils import flatten_dict, nulls_in_dict, get_commented_rows

__author__ = 'ARJ'

ENSEMBLE = 'ensemble'
VARDICT = 'vardict'
VARSCAN = 'varscan'
FREEBAYES = 'freebayes'
MUTECT = 'mutect'

READS = 'reads'
ALT_READS = 'alt_reads'
ALT_FREQ = 'alt_freq'

MINIMAL_REPORT_MAP = {
    'base_qual': ['mutect.tumor.BQ', 'vardict.tumor.QUAL'],
    'mapping_qual': ['freebayes.MQ', 'vardict.tumor.MQ', 'varscan.MQ']
}
COLUMNS_FIRST = ['chrom', 'start', 'end', 'ref', 'alt', 'called_by', 'caller_count', 'base_qual', 'mapping_qual']


allcts = set()


# protein affecting
pam_set_snpeff = [
    'FRAME_SHIFT',
    'NON_SYNONYMOUS_CODING',
    'STOP_GAINED',
    'NON_SYNONYMOUS_START',
    'START_LOST',
    'STOP_LOST',
    'START_GAINED',
    'CODON_DELETION',
    'CODON_INSERTION',
    'CODON_CHANGE_PLUS_CODON_INSERTION',
    'CODON_CHANGE_PLUS_CODON_DELETION',
    'PROTEIN_PROTEIN_CONTACT',
    'SPLICE_SITE_ACCEPTOR',
    'SPLICE_SITE_DONOR'
]


# protein affecting ensembl
pam_set_ensembl = [
    'FRAMESHIFT_VARIANT',
    'NON_SYNONYMOUS_CODING',
    'STOP_GAINED',
    'STOP_RETAINED_VARIANT',
    'NON_SYNONYMOUS_START',
    'START_LOST',
    'STOP_LOST',
    'START_GAINED',
    '5_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT',
    'MISSENSE_VARIANT',
    'DISRUPTIVE_INFRAME_DELETION',
    'DISRUPTIVE_INFRAME_INSERTION',
    'CONSERVATIVE_INFRAME_DELETION',
    'CONSERVATIVE_INFRAME_INSERTION',
    'STRUCTURAL_INTERACTION_VARIANT',
    'SPLICE_SITE_ACCEPTOR',
    'SPLICE_SITE_DONOR',
    'SPLICE_DONOR_VARIANT'
]

non_pam_set = [
    'SYNONYMOUS_CODING',
    'SYNONYMOUS_VARIANT',
    'SPLICE_SITE_REGION',
    'SPLICE_REGION_VARIANT',
    'SPLICE_ACCEPTOR_VARIANT',
    'SYNONYMOUS_START',
    'SYNONYMOUS_STOP',
    'NON_CODING_TRANSCRIPT_VARIANT',
    'TRANSCRIPT',
    'UTR_5_PRIME',
    'UTR_3_PRIME',
    '3_PRIME_UTR_VARIANT',
    '5_PRIME_UTR_VARIANT',
    'TF_BINDING_SITE_VARIANT',
    'UTR_5_DELETED',
    'SEQUENCE_FEATURE'
]

non_pc_set = [
    'EXON',
    'NON_CODING_TRANSCRIPT_EXON_VARIANT',
    'TRANSCRIBED_UNPROCESSED_PSEUDOGENE',
    'PROCESSED_TRANSCRIPT'
]

skip_set = [
    'DOWNSTREAM',
    'UPSTREAM',
    'INTERGENIC',
    'INTERGENIC_REGION',
    'DOWNSTREAM_GENE_VARIANT',
    'UPSTREAM_GENE_VARIANT',
    'INTRON',
    'INTRON_VARIANT',
    'DOWNSTREAM_GENE_VARIANT',
    'UPSTREAM_GENE_VARIANT'
]

pam_set = pam_set_ensembl + pam_set_snpeff
ct_set = pam_set + non_pam_set + non_pc_set

logger = logging.getLogger()
#logger.debug("CTSET: {}".format(ct_set))
#logger.debug("SKIP_SET: {}".format(skip_set))


TSL_df = (
    pd.read_table('/home/mpschr/Data/transcript_support_level/wgEncodeGencodeTranscriptionSupportLevelV23.txt.gz',header=None)
    .rename(columns={0: 'ENST', 1: 'TSL'})
    .assign(ENST = lambda df: df.ENST.apply(lambda enst: enst.split('.')[0]))
    .assign(TSL = lambda df: df.TSL.replace(-1,99))
    .set_index('ENST')
)


def is_pam_ct(ct):

    if ct == '' or ct is None:
        return False

    ct = ct.upper()

    if 'NEXT_PROT' in ct or 'MOTIF' in ct:
        return False

    if ct not in ct_set and ct not in skip_set:
        #print('unclassifiable CT: {}'.format(ct))
        raise RuntimeError('unclassifiable CT: {}. Please classify the consequence type in the (this) source code'.format(ct))
    return ct in pam_set



def consequence_types(cts):
    candidates = set()
    for ct in cts.split(','):
        if ct.split('|')[7] == 'protein_coding' or True:
            gene = ct.split('|')[3]
            ctstrings = ct.split('|')[1]
            if "NEXT_PROT" not in ctstrings:
                for ctstring in ctstrings.replace('+','&').split("&"):
                    #if is_protein_coding_ct(ctstring):
                    candidates.add(':'.join([gene,ctstring.lower()]))

    return ','.join(list(candidates))


def select_most_impacting_gene(ct_strings, coding_or_not):
    winners = {}
    winner_idx = -1
    winner_aa_change = None
    winner_consequence = ''
    winner_biotype = ''
    winner_has_warning = False
    winner_warning = ''
    winner_length = -1
    winner_transcript = ''
    winner_tsl = -1
    winner_sift = ''
    if ct_strings == "" or pd.isnull(ct_strings):
        return None

    for transcript in ct_strings.split(','):
        fields = transcript.replace(')', '').replace('(', '|').split('|')
        #Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank |
        # 10: HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
        if len(fields) == 17:
            #['Allele', Annotation,  'Impact'    , '_Name', '_ID ',         'Feature_Type', Feature_ID      Transcript_BioType   Rank   HGVS.c           HGVS.p     cDNA.pos     CDS.pos AA.pos
            #['G', 'intron_variant', 'MODIFIER', 'TTC34', 'ENSG00000215912', 'transcript', 'ENST00000401095', 'protein_coding', '4/6', 'c.688-12535G>C',  '',           '',        '',           '',      '', '']
            allele, ct, sift, gene, geneid, feature, transcript, biotype, rank, codon_change, aa_change, cdna_pos, cds_pos, aa_pos, dist , warning, kk = fields

        else:
        #    ct,sift,muttype,codon_change,aa_change,aa_length,gene,biotype,coding,transcript,exon,base = fields
            allele, ct, sift, gene, geneid, feature, transcript, biotype, rank, codon_change, aa_change, cdna_pos, cds_pos, aa_pos, dist , warning = fields

        if ct.upper() in skip_set:
            continue

        if feature != 'transcript':
            continue

        has_warning = warning is not None and warning != ''


        #is_non_coding_var = aa_change == "" or aa_change == None
        is_non_coding_var = all([is_pam_ct(ctsplit) == False for ctsplit in ct.split('&')])
        if coding_or_not == 'CODING' and is_non_coding_var or coding_or_not != 'CODING' and not is_non_coding_var:
            continue

        if 'NEXT_PROT' in ct or 'MOTIF' in ct:
            continue
        if '+' in ct:
            # select most impacting of multiple effects within 1 transcript
            ranks = [ct_set.index(ctplus.upper()) for ctplus in ct.split('+') if ctplus.upper() not in skip_set]
            ct = ct.split('+')[ranks.index(min(ranks))]
        if '&' in ct:
            # select most impacting of multiple effects within 1 transcript
            ranks = [ct_set.index(ctplus.upper()) if ctplus.upper() not in skip_set else 99 for ctplus in ct.split('&')]

            #print(ranks)

            if len(ranks) == 0:
                logging.debug('setting ct to "" {}'.format(ct))
                ct = ''
            else:
                ct = ct.split('&')[ranks.index(min(ranks))]
            #logging.debug('selected {}'.format(ct))


        idx = ct_set.index(ct.upper())

        try:
            aa_length = int(aa_pos.split('/'))
        except :
            aa_length = 0

        try:
            tsl = TSL_df.loc[transcript].TSL
        except KeyError:
            tsl = 99

        if idx == winner_idx and gene not in winners:
            winners.add(gene)
        elif (winner_idx < 0 or idx < winner_idx) or winner_idx == idx and (aa_length > winner_length or
                                                                            winner_has_warning and not has_warning or
                                                                            tsl < winner_tsl):
            winner_idx = idx
            winners = {gene}
            winner_consequence = ct_set[idx].lower()
            winner_aa_change = aa_change
            winner_biotype = biotype
            winner_has_warning = has_warning
            winner_warning = warning
            winner_transcript = transcript
            winner_length = aa_length
            winner_tsl = tsl
            winner_sift = sift.lower()

    if coding_or_not == 'CODING':
        return {'genes': ','.join(list(winners)),
                'aa_change': winner_aa_change,
                'consequence': winner_consequence,
                'impact': winner_sift,
                'biotype': winner_biotype,
                'transcript': winner_transcript,
                'warning': winner_warning,
                'tsl': winner_tsl}
    elif coding_or_not == 'NON_CODING':
        return {'genes_noncoding': ','.join(list(winners)),
                'consequence_noncoding': winner_consequence,
                'transcript_noncoding': winner_transcript,
                'biotype_noncoding': winner_biotype}
    else:
        raise RuntimeError('CODING OR NON_CODING?! ')

def longest_transcript_info(row):
    if (row.pam_genes_count == 0) and row.impact != 'synonymous_coding' and row.aa_change != '':
        #if coding == "CODING" and biotype != "protein_coding" and aa_change == "":
        return

    max_exon = -1
    max_aalength = -1
    max_exon_transcript = ''
    max_aachange = ''
    max_gene = ''
    max_transcript = ''
    for transcript in row.EFF.split(','):
        fields = transcript.replace(')', '').replace('(', '|').split('|')
        if len(fields) == 17:
            allele, ct, impact, gene, geneid, feature, transcript, biotype, rank, codon_change, aa_change, cdna_pos, cds_pos, aa_pos, dist , warning, error = fields

        else:
            allele, ct, impact, gene, geneid, feature, transcript, biotype, rank, codon_change, aa_change, cdna_pos, cds_pos, aa_pos, dist , warning = fields

        if biotype != 'protein_coding':
            continue


        if feature != 'transcript':
            continue

        if 'NEXT_PROT' in impact or 'MOTIF' in impact:
            continue

        if '+' in impact:
            # select most impacting of multiple effects within 1 transcript
            ranks = [ct_set.index(ctplus.upper()) for ctplus in impact.split('+') if ctplus.upper() not in skip_set]
            impact = impact.split('+')[ranks.index(min(ranks))]

        if '&' in impact:
            # select most impacting of multiple effects within 1 transcript
            ranks = [ct_set.index(ctplus.upper()) for ctplus in impact.split('&') if ctplus.upper() not in skip_set]
            impact = impact.split('+')[ranks.index(min(ranks))]

        if not impact in pam_set and impact.lower() != row.impact:
            continue

        if not gene in row.genes:
            continue

        exon = 0
        if exon == '' and aa_change == '':
            continue
        try:
            aa_length = int(aa_pos.split('/')[1])
        except:
            aa_length=0
        if exon > max_exon or exon == max_exon and aa_length > max_aalength:
            max_exon = exon
            max_aalength = aa_length
            max_exon_transcript = transcript
            max_aachange = aa_change
            max_gene = gene
            max_transcript = transcript

    if max_exon < 0:
        max_exon = None
    if max_aalength < 0:
        max_aalength = None

    return {'exon_longest': max_exon,
            'aa_length_longest': max_aalength,
            'aa_change_longest': max_aachange,
            'transcript_longest': max_transcript,
            'gene_longest': max_gene }


class MultiVCFLoader():
        def __init__(self, folder, sample):
            self.files = [
                pathjoin(folder, "b{}-ensemble.vcf.gz".format(sample)),
                pathjoin(folder, "b{}-mutect.vcf.gz".format(sample)),
                pathjoin(folder, "b{}-vardict.vcf.gz".format(sample)),
                pathjoin(folder, "b{}-varscan.vcf.gz".format(sample)),
                pathjoin(folder, "b{}-freebayes.vcf.gz".format(sample))
            ]

        def get_df(self, debug=False):
            dfs = {}
            self.unclassifiable_cts = set()

            # load the 5 vcf files
            files = self.files
            if debug:
                files = self.files[0:2]
            for f in files:
                skiprows = get_commented_rows(f)
                caller = os.path.basename(f).split("-")[-1].replace(".vcf.gz", "")
                status_update = "{}/{}: {}               ".format(self.files.index(f) + 1, len(self.files), caller)
                print(status_update, end='\r')
                d = (pd.read_table(f, skiprows=skiprows, dtype=str)
                     .rename(columns={"#CHROM": "CHROM"})
                     .query("FILTER == 'PASS'")
                     .assign(POS=lambda df: df.POS.astype(int))
                     .set_index(["CHROM", "POS", "REF", "ALT"])
                     )
                d['caller'] = caller
                dfs[caller] = d

            # merge the vcf files
            df_merge = None
            callerlist = ['ensemble', 'vardict', 'varscan', 'freebayes', 'mutect']
            if debug:
                callerlist = ['ensemble', 'mutect']
            for c in callerlist:
                if type(df_merge) != pd.core.frame.DataFrame:
                    d = dfs[c]
                    df_merge = d[['caller']].rename(columns={'caller': c}).reset_index()
                else:
                    d = dfs[c][dfs[c].index.isin(dfs['ensemble'].index)]
                    df_merge = df_merge.merge(d[['caller']].rename(columns={'caller': c}).reset_index(),
                                              how='left')
            # create the called df
            idx = ["CHROM", "POS", "REF", "ALT"]
            df_called = (
                df_merge
                    .fillna("")
                    .drop('ensemble', axis=1)
                    .set_index(idx)
                    .apply(lambda row: ",".join((sorted([x for x in row.tolist() if x != ""]))), axis=1)
                    .to_frame().rename(columns={0:'called_by'})
                    .assign(caller_count = lambda df: df.apply(lambda row: len(row['called_by'].split(",")), axis=1))
                    .reset_index()
                    .drop_duplicates(subset=idx)
                    .set_index(idx)
            )
            del df_merge

            try:
                assert sum(df_called.index == dfs['ensemble'].index)/dfs['ensemble'].shape[0] == 1
            except ValueError:
                id1 = df_called.index
                id2 = dfs['ensemble'].index
                print(df_called.index.shape, dfs['ensemble'].index.shape)
                print(df_called[id1.duplicated(keep=False)])
                exit(-1)

            assert sum(df_called.index == dfs['ensemble'].index)/dfs['ensemble'].shape[0]

            df_called['EFF'] = dfs['ensemble'].INFO.str.extract('ANN=(.+?);')
            #df_called['type'] = (dfs['ensemble'].INFO + ';').str.extract('TYPE=(.+?);')
            df_called['ct_vcf'] = df_called.EFF.apply(consequence_types)

            df_called.insert(2, 'genes_possibly_affected', df_called.ct_vcf.apply(lambda ctvcf: ','.join({x.split(':')[0] for x in ctvcf.split(',')})) )

            ct_pam = df_called.ct_vcf.apply(lambda ctstring: ','.join([x for x in ctstring.split(',') if len(x) > 0 and is_pam_ct(x.split(":")[1]) ]))
            df_called.insert(0, 'pam_genes_count', ct_pam.apply(lambda ctpam: len({x.split(':')[0] for x in ctpam.split(',')}) if ctpam != '' else 0))

            most_impacted_noncoding = df_called.EFF.apply(lambda x: pd.Series(select_most_impacting_gene(x, 'NON_CODING')))
            for s in reversed(['genes_noncoding', 'consequence_noncoding', 'biotype_noncoding', 'transcript_noncoding']):
                df_called.insert(3, s, most_impacted_noncoding[s])

            most_impacted_coding = df_called.EFF.apply(lambda x: pd.Series(select_most_impacting_gene(x, 'CODING')))
            for s in reversed(['genes', 'aa_change', 'impact', 'consequence', 'biotype', 'transcript', 'warning', 'tsl']):
                df_called.insert(0, s, most_impacted_coding[s])

            # update pam_genes_count
            df_called['pam_genes_count'] = df_called.apply(lambda row: len(row.genes.split(',')) if is_pam_ct(row.consequence) else 0, axis=1)

            if len(self.unclassifiable_cts) > 0:
                raise RuntimeError("Untranslatable CTs! abort!")

            pc_mut_info = df_called.apply(lambda x: pd.Series(longest_transcript_info(x)),axis=1).drop('gene_longest', axis=1)

            for s in pc_mut_info.columns:
                df_called.insert(7, s, pc_mut_info[s])

            df_called.drop('EFF',axis=1,inplace=True)

            return df_called



class MultiVCFCall(object):

    def __init__(self, call_dict):
        self.call_dict = call_dict
        self.count, self.called_by = self.callers()
        self.qual_info = self.quality_query()
        self.genotype_info = self.genotypes_query()

    def callers(self):
        pass_callers = []
        count = 0
        filters = []
        for c in [VARDICT, VARSCAN, FREEBAYES, MUTECT]:
            if self.call_dict[c] is not None and 'REJECT' not in self.call_dict[c].FILTER:
                count += 1
                pass_callers.append(c)
                pass_callers.sort()
        return count, ','.join(pass_callers)


    def genotypes_query(self):
        genotype_info = {}
        for sample_obj in self.call_dict[ENSEMBLE].samples:
            s = sample_obj.sample
            if s == '':
               continue
            genotype_info[s] = {}
            for c in self.called_by.split(','):
                try:
                    geno = self.call_dict[c].genotype(s)
                except KeyError:
                    print(c,s,self.call_dict[ENSEMBLE])
                    continue
                reads = -1
                alt_reads = -1
                if c == VARDICT:
                    reads = geno.data.DP
                    alt_reads = geno.data.AD[1]
                if c == FREEBAYES:
                    reads = geno.data.DP
                    ao = geno.data.AO
                    if type(ao) != list:
                        alt_reads = ao
                    elif len(ao) == 1:
                        alt_reads = geno.data.AO[0]
                if c == VARSCAN:
                    dp4 = [int(x) for x in geno.data.DP4.split(',')]
                    reads = np.sum(dp4)
                    alt_reads = np.sum([dp4[2:4]])
                if c == MUTECT:
                    reads = geno.data.DP
                    alt_reads = geno.data.AD[1]
                try:
                    genotype_info[s][c] = {'reads': int(reads), 'alt_reads': int(alt_reads), 'alt_freq': int(alt_reads)/int(reads)}
                except Exception:
                    if c != FREEBAYES and '_CR' not in s:
                        print(s,c,self.call_dict[c],geno.data)
        return genotype_info

    def quality_query(self):
        info_fields = ['MQ', 'MQRankSum', 'BQ', 'BaseQRankSum', 'QUAL', 'QD', 'QSTD', 'GQ']
        qualinfo = {}
        for caller, call in self.call_dict.items():
            if caller == ENSEMBLE:
                continue

            if call is not None:
                if call.QUAL is not None:
                    qualinfo[caller + '.QUALITY'] = call.QUAL

                for f in info_fields:
                    if f in call.INFO and call.INFO[f] is not None:
                        qualinfo[caller + '.' + f] = call.INFO[f]
                    if '_CR' in call.samples[0].sample:
                        tumor_sample = call.samples[1]
                        normal_sample = call.samples[0]
                    else:
                        tumor_sample = call.samples[0]
                        normal_sample = call.samples[1]
                    if f in tumor_sample.data._fields:
                        qualinfo[caller + '.tumor.' + f] = getattr(tumor_sample.data, f)

        return qualinfo


class MultiVCFWalker(Iterator):

    def __init__(self, folder, sample):

        self.current = None
        self.count = -1
        self.called_by = None
        self.qual_info = None
        self.sample = sample

        self.order = [ENSEMBLE, VARDICT, VARSCAN, FREEBAYES, MUTECT]

        ensembl_reader = vcf.Reader(filename=pathjoin(folder, "b{}-ensemble.vcf.gz".format(sample)), compressed=True)
        vardict_reader = vcf.Reader(filename=pathjoin(folder, "b{}-vardict.vcf.gz".format(sample)), compressed=True)
        varscan_reader = vcf.Reader(filename=pathjoin(folder, "b{}-varscan.vcf.gz".format(sample)), compressed=True)
        freebayes_reader = vcf.Reader(filename=pathjoin(folder, "b{}-freebayes.vcf.gz".format(sample)), compressed=True)
        mutect_reader = vcf.Reader(filename=pathjoin(folder, "b{}-mutect.vcf.gz".format(sample)), compressed=True)
        self.walker = vcf.utils.walk_together(ensembl_reader, vardict_reader, varscan_reader, freebayes_reader, mutect_reader)

    def __next__(self) -> MultiVCFCall:
        while True:
            r = next(self.walker)
            if not r[0] is None:
                self.current = MultiVCFCall({ ENSEMBLE: r[0], FREEBAYES: r[3], MUTECT: r[4], VARDICT: r[1], VARSCAN: r[2] })
                return self.current


class MultiVCFtoDF:

    def __init__(self, multi_vcf_walker: MultiVCFWalker):
        self.multi_vcf_walker = multi_vcf_walker
        self.sample = multi_vcf_walker.sample
        self.tumor_time_point = self.sample.split('_')[-1]
        self.patient = self.sample.replace(self.tumor_time_point, '')

    def choose_readcount_data(self, read_counts_df, time_points, priority=[VARDICT, FREEBAYES, VARSCAN, MUTECT]):
        if type(time_points) != list:
            time_points = [time_points]
        logging.debug("multivcf time points: {}".format(time_points))

        records = []
        for i, row in read_counts_df.iterrows():
            r = {}
            for tp in time_points:
                tp_reads = '{}.{}'.format(tp, READS)
                tp_alt_reads = '{}.{}'.format(tp, ALT_READS)
                tp_alt_freq = '{}.{}'.format(tp, ALT_FREQ)
                tp_source = '{}.source_readcounts'.format(tp)
                tpr = ({tp_reads: None,
                          tp_alt_reads: None,
                          tp_alt_freq: None,
                          #tp_source: None
                         })

                for source in priority:
                    tpr[tp_reads] = row['{}.{}.{}'.format(tp, source, READS)]
                    tpr[tp_alt_reads] = row['{}.{}.{}'.format(tp, source, ALT_READS)]
                    tpr[tp_alt_freq] = row['{}.{}.{}'.format(tp, source, ALT_FREQ)]
                    if not nulls_in_dict(tpr):
                        tpr[tp_source] = source
                        r.update(tpr)
                        break
            records.append(r)
        return pd.DataFrame.from_records(records)

    def get_df(self, limit=0, extended=False, normal_time_point='CR') -> pd.DataFrame:

        data = []
        genotype_data = []

        logging.info('reading all vcf files simultaneously')

        start_time = time.time()
        c = 0
        df = None
        additional_cols = []
        for current_calls in self.multi_vcf_walker:

            genotype_data.append(flatten_dict(current_calls.genotype_info))

            ensemble_call = current_calls.call_dict[ENSEMBLE]
            current_info = {'chrom': ensemble_call.CHROM, 'ref': ensemble_call.REF, 'alt': str(ensemble_call.ALT[0]),
                            'start': ensemble_call.start, 'end': ensemble_call.end,
                            'caller_count': current_calls.count, 'called_by': current_calls.called_by}

            q_info = {}
            for k, fields in MINIMAL_REPORT_MAP.items():
                try:
                    choose_from = [current_calls.qual_info[x] for x in fields
                                   if x in current_calls.qual_info and current_calls.qual_info[x] is not None]
                    q_info[k] = np.nanmax(choose_from)
                except ValueError:
                    pass

            current_info.update(q_info)

            if extended:
                current_info.update(current_calls.qual_info)
            data.append(current_info)

            c += 1
            if limit > 0 and c % limit == 0:
                logging.info(c)
                break
            elif c % 100 == 0:
                current_time = time.time()
                elapsed = current_time - start_time
                mins = math.floor(elapsed/60)
                secs = math.ceil(elapsed % 60)
                print('At variant {} - {}m{}s'.format(c, mins,secs),  end='\r')

        df_general = pd.DataFrame(data)
        df_genotype = pd.DataFrame(genotype_data)

        #remove sample specific string e.g. AE03
        df_genotype.columns = df_genotype.columns.str.replace(self.patient, '')

        #choose best genotype info
        df_genotype = self.choose_readcount_data(df_genotype, time_points=[self.tumor_time_point, normal_time_point])
        df_genotype.columns = [x.split('.')[1] + '_' + x.split('.')[0] for x in df_genotype.columns]

        # concat dataframes
        df = pd.concat([df_general, df_genotype], axis=1)
        #logger.debug(df.head())
        additional_cols = [x for x in df.columns if x not in COLUMNS_FIRST]
        additional_cols.sort()

        return df[COLUMNS_FIRST + additional_cols]
