import argparse
import itertools
import logging
import os
import random
from matplotlib import gridspec
from pandas.tools.plotting import parallel_coordinates


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from allpy.SimpleTask import SampleTask, PatientTask

from allpy.exon.mutations.mutreport_run import MutationReporter
from allpy.pipelineutils.sample_info import get_sample_info

PLOT_VAF2D = 'vaf-2d'

PLOT_VAF1D = 'vaf-1d'

PLOT_VIOLIN = 'violin'

PLOT_CALLERS = 'callers'

PLOT_ALL = 'all'




class ExonMutationPlotter(PatientTask):

    def derive_patient_name(self, task_input, cmd_patient_name):
        if cmd_patient_name is not None:
            patient = cmd_patient_name
        else:
            patient = os.path.basename(task_input)
        return patient


    def add_metadata(self, metastring, do_log=True):
        self.metadata.append(metastring)
        if do_log:
            logging.info(metastring)

    def __init__(self, patient_dir):
        super().__init__('EXON_MUTATION_PLOTTER')
        self.metadata = []
        self.blasten_col = "purity"
        self.blasten = get_sample_info(patient_dir)[[self.blasten_col]].fillna(100)
        self.blasten[self.blasten_col] = self.blasten[self.blasten_col].astype(int)


    def run(self, patientdir, patient=None, toplot=[PLOT_ALL], time_points=None, time_point_normal='CR'):
        super().run(patientdir, '.', patient)

        if time_points == None or time_points == []:
            time_points = ["ID", "REL"]
        logging.debug(time_points)

        if patient is None:
            patient = os.path.basename(patientdir)

        if toplot is None:
            toplot=[PLOT_ALL]

        reporter = MutationReporter(patient, patientdir)

        sns.set_style("whitegrid")

        df = read_data(patient, reporter, time_points, time_point_normal)

        # normalize time points:
        if 'alt_freq_REL2' in df:
            df['alt_freq_REL'] = df.alt_freq_REL2
        if time_point_normal != "CR":
            df['alt_freq_CR'] = df["alt_freq_{}".format(time_point_normal)]

        #violinplot
        if PLOT_VIOLIN in toplot or PLOT_ALL in toplot:
            logging.info('violin plot')
            #df = reporter.get_standard_df(time_points, time_point_normal, True)
            df, fixed = self.blasts_fix(df)
            plot_violin(df, reporter, self.patient)


            plot_scatter(df, reporter, self.patient)

            plot_parcoord(df, reporter, self.patient)

        exit()

        if PLOT_CALLERS in toplot or PLOT_ALL in toplot:
            logging.info('called mutations plot')
            plot_callerhist(id_df, patient, rel_df, reporter, time_points)

        try:

            bothm, idm, relm, to_highlight = sciclone_data_prep(id_df, patient, rel_df, reporter, time_points=time_points)

            if PLOT_VAF1D in toplot or PLOT_ALL in toplot:

                logging.info('vaf 1d plot')
                plot_vaf_1d(idm, patient, relm, reporter)

            if PLOT_VAF2D in toplot or PLOT_ALL in toplot:

                logging.info('vaf 2d plot')
                plot_vaf_2d(bothm, patient, reporter, to_highlight)

        except AssertionError:
            pass
        #except Exception:
        #    logging.error("Sciclone plots failed")

    def blasts_fix(self, df):

        id_sample = "{}_ID".format(self.patient)
        rel_sample = "{}_REL".format(self.patient)
        id_blasten = self.blasten.loc[id_sample, self.blasten_col] / 100 if id_sample in self.blasten.index else 1.0
        rel_blasten = self.blasten.loc[rel_sample,self.blasten_col] / 100 if rel_sample in self.blasten.index else 1.0

        df['alt_freq_ID_fix'] = df.alt_freq_ID / id_blasten
        df['alt_freq_REL_fix'] = df.alt_freq_REL / rel_blasten

        over1 = df[(df.alt_freq_ID_fix > 1) | (df.alt_freq_REL_fix > 1)][["patient","genes","alt_freq_ID_fix", 'alt_freq_CR',"alt_freq_REL_fix", "known_cancer_gene"]]
        df.loc[over1.index,("alt_freq_ID_fix","alt_freq_REL_fix")] = 1
        return df, over1


def complete_clone(both, ID, REL):
    individual_clones = []
    individual_clones.append("{:.0f}".format(int(both)))
    if int(ID) > 0:
        individual_clones.append("ID{:.0f}".format(int(ID)))
    if int(REL) > 0:
        individual_clones.append("REL{:.0f}".format(int(REL)))
    return "-".join(individual_clones)
    #return "{:.0f}-ID{:.0f}-REL{:.0f}".format(int(both), int(ID), int(REL))

def get_markers(clones):

    markers = []
    c = None
    for x,group in itertools.groupby(clones, key=lambda x: x.split('-')[0]):
        if c is None:
            c = itertools.cycle(['v','8','*', '<', '>'])
        else:
            c = itertools.cycle(['o','s','h','d', '*', '8'])
        for g in group:
            markers.append(next(c))

    return markers


def calledbysplit(calledby):
    freebayes = 1 if 'freebayes' in calledby else 0
    mutect = 1 if 'mutect' in calledby else 0
    vardict = 1 if 'vardict' in calledby else 0
    varscan = 1 if 'varscan' in calledby else 0
    three = 0
    four = 0
    if len(calledby.split(',')) == 3:
        three = 1
    elif len(calledby.split(',')) == 4:
        three = 1
        four = 1

    return pd.Series({'freebayes': freebayes, 'mutect': mutect,
                      'vardict': vardict, 'varscan': varscan,
                      'called by 2': 1, 'called by 3' : three, 'called by 4': four})

def get_colors(clones, palette):
    p = sns.color_palette(palette, len(clones))
    return [p[int(x.split('-')[0])] for x in clones]



def plot_callerhist(id_df, patient, rel_df, reporter, time_points):
    
    ID, REL = time_points #may be other time points than ID and REL    
    
    rel_df_calledby = rel_df.called_by.apply(calledbysplit)
    id_df_calledby = id_df.called_by.apply(calledbysplit)
    id_cmelt = pd.melt(pd.merge(id_df, id_df_calledby, left_index=True, right_index=True),
                       var_name='caller', value_name='called',
                       value_vars=id_df_calledby.columns.tolist(),
                       id_vars=['alt_freq_{}'.format(ID), 'reads_{}'.format(ID), 'gene', 'alt_reads_{}'.format(ID)])
    rel_cmelt = pd.melt(pd.merge(rel_df, rel_df_calledby, left_index=True, right_index=True),
                        var_name='caller', value_name='called',
                        value_vars=id_df_calledby.columns.tolist(),
                        id_vars=['alt_freq_{}'.format(REL), 'reads_{}'.format(REL), 'gene', 'alt_reads_{}'.format(REL)])

    #bin the daa
    binsize = 10
    id_cmelt['reads_binned'] = id_cmelt['alt_reads_{}'.format(ID)].apply(lambda x: binreads(x, binsize))
    rel_cmelt['reads_binned'] = rel_cmelt['alt_reads_{}'.format(REL)].apply(lambda x: binreads(x, binsize))
    id_cgrouped = id_cmelt.groupby(['caller', 'reads_binned']).agg({'called': np.sum})
    rel_cgrouped = rel_cmelt.groupby(['caller', 'reads_binned']).agg({'called': np.sum})

    #fill the gaps
    tools = ['freebayes', 'mutect', 'varscan', 'vardict', 'called by 2', 'called by 3', 'called by 4']
    start = min(id_cgrouped.index.levels[1].min(),rel_cgrouped.index.levels[1].min())
    end = max(id_cgrouped.index.levels[1].max(),rel_cgrouped.index.levels[1].max())

    for x in range(start, end+binsize, binsize):
        for t in tools:
            if x not in id_cgrouped.loc[t].index:
                id_cgrouped.loc[(t,x),'called'] = 0
            if x not in rel_cgrouped.loc[t].index:
                rel_cgrouped.loc[(t,x),'called'] = 0


    #add cummulative sum
    id_cgrouped['sample'] = '{}_{}'.format(patient, ID)
    rel_cgrouped['sample'] = '{}_{}'.format(patient, REL)
    cgrouped = pd.concat([id_cgrouped.reset_index(), rel_cgrouped.reset_index()], axis=0).sort(['sample', 'caller', 'reads_binned']).set_index(['sample', 'caller'])
    cgrouped['called_cum'] = cgrouped.groupby(level=[0,1])['called'].cumsum()


    #plot
    plt.clf()
    flatui = ["pink", "pink", "pink", "windows blue", "light orange", "dusty red", "green"]
    pallete = sns.xkcd_palette(flatui)
    sns.axes_style("darkgrid")
    fig = plt.figure(figsize=(20, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1])

    ax = plt.subplot(gs[0,0])
    g = sns.pointplot(x="reads_binned", y="called_cum", data=cgrouped.loc[(patient + '_' + ID)].reset_index(), hue="caller", linewidth=0, ax=ax, ci=0,
                      kind="point", palette=pallete, alpha=0.8,  dodge=0.3, markers=['2', '3', '4' ,'s','o','d', 'H', '8'],
                     linestyles=[':', ':', ':' ,'-','-','-', '-', '-'])
    g.set_ylabel("Called Somatic PC Mutations")
    ax = plt.subplot(gs[1,0])
    g = sns.pointplot(x="reads_binned", y="called_cum", data=cgrouped.loc[(patient + '_' + REL )].reset_index(), hue="caller", linewidth=0, ax=ax, ci=0,
                      kind="point", palette=pallete,  alpha=0.8,  dodge=0.3, markers=['2', '3', '4' ,'s','o','d', 'H', '8'],
                     linestyles=[':', ':', ':' ,'-','-','-', '-', '-'])
    g.set_ylabel("Called Somatic PC Mutations")

    #g.set_axis_labels("Alt. read count", "Called Somatic PC Mutations").set_titles("Detected mutations by caller and binned alt-reads for {row_var} {row_name}")



    ax = plt.subplot(gs[0,1])
    g = sns.barplot(x="caller", y="called", data=id_cmelt.groupby('caller').sum().reset_index(), palette=pallete)
    g.set_ylabel("")


    ax = plt.subplot(gs[1,1])
    g = sns.barplot(x="caller", y="called", data=rel_cmelt.groupby('caller').sum().reset_index(), palette=pallete)
    g.set_ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(reporter.out_folder, '{}_mutationcallers.coding.png'.format(patient)))


def binreads(reads, size=10):
    return (int(reads/size)+1)*size



def plot_vaf_1d(idm, patient, relm, reporter, time_points):

    ID, REL = time_points #may be other time points than ID and REL

    plt.clf()
    idm['sample'] = "{}_{}".format(patient, ID)
    relm['sample'] = "{}_{}".format(patient, REL)
    idm.rename(columns={'alt_freq_{}'.format(ID): 'alt_freq'}, inplace=True)
    relm.rename(columns={'alt_freq_{}'.format(REL): 'alt_freq'}, inplace=True)
    m = pd.concat([idm, relm])
    m.clone.fillna(0, inplace=True)
    m.clone = m.clone.astype(int).astype(str)
    m.clone[pd.isnull(m.clone)] = "none"
    m.sort(columns='alt_freq').head()
    # clones 1d plot
    f, axes = plt.subplots(1, 1, figsize=(12, 6), sharex=False)
    highlights = pd.read_table("/home/mpschr/Documents/projects/exon/exonR/highlightedgenes.txt")

    p = sns.stripplot(x="alt_freq", y='sample', data=m.sort(columns=['clone']),  hue='clone', split=True,
                       jitter=True, palette=sns.color_palette("Paired"), size=12, alpha=0.6, ax=axes)
    axes.set_xlim(0,1.2)

    to_highlight = m[m.gene.isin(highlights.GeneSymbol)].set_index('gene').sort(columns=['sample','clone', 'alt_freq'], ascending=[False, False, True])[['alt_freq', 'clone', 'sample']]
    label_y = 1
    maxclone = int(m.clone.max())
    clonespace = {6: 0.11, 5: 0.13, 4: 0.15, 3: 0.175, 2: 0.225, 1: 0.24, 0: 0.5}[maxclone]
    for gene, coords in to_highlight.iterrows():
        #print([gene]+coords.tolist())
        x,clone,sample = coords.tolist()
        s = 0 if ID in sample else 1
        plt.annotate(gene, xy=((x,s - 0.3 + int(clone)*clonespace)), xytext=(1.1, label_y+s/15), annotation_clip=True,
                     arrowprops=dict(facecolor='black', shrink=0.005, width=1, frac=0.01, headwidth=5, alpha=0.1))
        label_y -= 0.07


    plt.tight_layout()
    plt.savefig(os.path.join(reporter.out_folder, '{}_clones1d.png'.format(patient)))


def plot_vaf_2d(bothm, patient, reporter, to_highlight, time_points):

    ID, REL = time_points #may be other time points than ID and REL

    plt.clf()
    both_clusters = bothm.clone.unique().tolist()
    both_clusters.sort()
    hue_kws = {
        "marker": get_markers(both_clusters)
    }
    label_y = 1
    g = sns.FacetGrid(bothm.sort(columns='clone'), palette=get_colors(both_clusters, "Paired"),
                      hue="clone",size=7, aspect=1, hue_kws=hue_kws, xlim=(-0.2,1.2), ylim=(-0.05,1))

    if to_highlight.shape[0] > 30:
        p1 = to_highlight.query("alt_freq_{} == 0 and alt_freq_{} > 0".format(REL, ID)).head(10)
        p2 = to_highlight.query("alt_freq_{} > 0 and alt_freq_{} > 0".format(REL, ID)).head(10)
        p3 = to_highlight.query("alt_freq_{} > 0 and alt_freq_{} == 0".format(REL, ID)).head(30-p1.shape[0]-p2.shape[0])
        to_highlight = pd.concat([p1,p2,p3])

    for gene, coords in to_highlight.sort('alt_freq_{}'.format(REL), ascending=False).iterrows():
        x,y = coords.tolist()
        xpos = 1.05 if x != 0 else -0.15
        plt.annotate(gene, xy=((x, y)), xytext=(xpos, label_y), annotation_clip=True,
                     arrowprops=dict(facecolor='black', shrink=0.005, width=1, frac=0.01, headwidth=5, alpha=0.1))
        label_y -= 0.04
    g.map(plt.scatter, "alt_freq_{}".format(ID), "alt_freq_{}".format(REL), alpha=.45, s=90, linewidth=1.5)
    g.add_legend()
    plt.savefig(os.path.join(reporter.out_folder, '{}_clones2d.png'.format(patient)))


def sciclone_data_prep(id_df, patient, rel_df, reporter, allmuts=True, time_points=None):

    ID, REL = time_points #may be other time points than ID and REL

    # data prep sciclone
    sciclone_cols = ['chr', 'st', 'cluster']
    subfolder = 'allmuts' if allmuts else 'coding'
    try:
        sc_id_df = pd.read_table(os.path.join(os.path.dirname(reporter.out_folder),
                                              'sciClone',
                                              subfolder,
                                              'clusters.{}_{}.1d.depth40-alt5-freq0.tsv'.format(patient, ID)))[sciclone_cols]
        sc_rel_df = pd.read_table(os.path.join(os.path.dirname(reporter.out_folder),
                                               'sciClone',
                                               subfolder,
                                               'clusters.{}_{}.1d.depth40-alt5-freq0.tsv'.format(patient, REL)))[sciclone_cols]
        sc_id_df.columns = ['chrom', 'end', 'clone']
        sc_rel_df.columns = ['chrom', 'end', 'clone']
        idm = pd.merge(id_df[['chrom', 'end', 'alt_freq_ID', 'gene']], sc_id_df, how='left')
        relm = pd.merge(rel_df[['chrom', 'end', 'alt_freq_REL', 'gene']], sc_rel_df, how='left')
        relm.clone.fillna(0, inplace=True)
        idm.clone.fillna(0, inplace=True)

    except OSError:
        logging.error("No 1d sciclone files found for {}, {} or both".format(ID, REL))
        idm = id_df[['chrom', 'end', 'alt_freq_{}'.format(ID), 'gene']]
        idm['clone'] = 0
        relm = rel_df[['chrom', 'end', 'alt_freq_{}'.format(REL), 'gene']]
        relm['clone'] = 0

    clones2dfile = os.path.join(os.path.dirname(reporter.out_folder),
                                'sciClone',
                                subfolder,
                        'clusters.{}.2d.depth40-alt5-freq0.tsv'.format(patient))

    #print(sc_both_df.columns)
    bothm = pd.merge(idm, relm, how='outer', on=['chrom', 'end'])

    both_df_cols = ['chrom', 'end', 'clone', '{}_CN'.format(ID), '{}_CN'.format(REL)]
    if os.path.exists(clones2dfile):
        sciclone_cols += ['{}_{}.cleancn'.format(patient, ID), '{}_{}.cleancn'.format(patient, REL)]
        sc_both_df = pd.read_table(clones2dfile)[sciclone_cols]
        sc_both_df.columns = both_df_cols
    else:
        logging.error("No 2d sciclone file found")
        sc_both_df = pd.DataFrame(columns=both_df_cols)
        sc_both_df['clone'] = 0

    bothm = pd.merge(bothm, sc_both_df, how='left')
    bothm.fillna(0, inplace=True)
    bothm.clone = bothm.clone.astype(int).astype(str)
    bothm.clone = bothm.apply(lambda x: complete_clone(x['clone'], x['clone_x'], x['clone_y']), axis=1)
    bothm['gene'] = bothm.apply(lambda x: x['gene_x'] if type(x['gene_x']) is not int else x['gene_y'], axis=1)

    # clusters bo2d plot
    highlights = pd.read_table("/home/mpschr/Documents/projects/exon/exonR/highlightedgenes.txt")
    to_highlight = bothm[bothm.gene.isin(highlights.GeneSymbol)].set_index('gene').sort(columns=['alt_freq_{}'.format(REL)], ascending=False)[['alt_freq_{}'.format(ID), 'alt_freq_{}'.format(REL)]]

    return bothm, idm, relm, to_highlight


def plot_violin(df, reporter, patient):

    plt.clf()
    f, axarr = plt.subplots(ncols=2, sharey=True)

    freq_cols = ['alt_freq_ID', 'alt_freq_REL']
    sub_plot_violin(df, freq_cols, ax=axarr[0])
    freq_cols = ['alt_freq_ID_fix', 'alt_freq_REL_fix']
    sub_plot_violin(df, freq_cols, ax=axarr[1], value_suffix='fix')
    plt.savefig(os.path.join(reporter.out_folder, '{}_violins.png'.format(patient)))



def sub_plot_violin(df, freq_cols, no_points=False, ax=None, value_suffix=''):


    sns.set(font_scale=1.4)
    plt.rc("figure", figsize=(13, 6.5))

    #df['mutations'] = 'pam'

    value_name = 'alt_freq' + '_' + value_suffix

    varfreq_df = (
        pd.melt(pd.concat([df], ignore_index=True)[['mutset', 'mutations'] + freq_cols],
                id_vars=['mutset', 'mutations'], value_vars=freq_cols, var_name='timepoint', value_name=value_name)
            .assign(timepoint=lambda df: df.timepoint.apply(lambda s: s.split('_')[2]))
            .query('{} > 0'.format(value_name))
    )

    sns.set_style("whitegrid")

    ax = sns.violinplot(x='timepoint', y=value_name, data=varfreq_df, hue='mutations', cut=0.2,
                        palette='Set1', scale='area',ax=ax)
    if not no_points:
        ax = sns.swarmplot(x="timepoint", y=value_name, data=varfreq_df, hue='mutations', split=True,
                              size=6, alpha=.6, edgecolor='black',palette="pastel", ax=ax) #color='yellow',
        ax.set_ylim((0, 1))


def plot_scatter(df, reporter, patient):
    plt.clf()
    lm = sns.lmplot("alt_freq_ID_fix", "alt_freq_REL_fix", data=df,  hue='known_cancer_gene',
                    markers=['x','o'],
                    fit_reg=False, size=7, scatter_kws={"s": 80, "alpha":0.5, "linewidth": 1.5})
    lm.axes[0,0].set_xlim(-.05,1.05)
    lm.axes[0,0].set_ylim(-0.05,1.05)

    label_y = 1
    prefix = 1
    for gene, row in df.query('known_cancer_gene').sort('alt_freq_ID_fix', ascending=True).set_index('genes').iterrows():
        x,y = row[['alt_freq_ID_fix', 'alt_freq_REL_fix']]

        #xpos = 1.05 if x != 0 else -0.15
        prefix = prefix * -1
        xpos = x + random.uniform(0.01, 0.05)  * prefix
        label_y = y + random.uniform(0.02, 0.05) * prefix
        plt.annotate(gene, xy=((x, y)), xytext=(xpos, label_y), annotation_clip=True, fontsize=10,
                     arrowprops=dict(facecolor='black', shrink=0.005, width=1, frac=0.01, headwidth=5, alpha=0.1))
        label_y -= 0.04
    lm.add_legend()

    plt.savefig(os.path.join(reporter.out_folder, '{}_scatter.png'.format(patient)))


def plot_parcoord(df, reporter, patient):
    plt.clf()
    plt.rc("figure", figsize=(13, 6.5))
    parallel_coordinates(df.query('known_cancer_gene'), 'genes', ['alt_freq_ID_fix', 'alt_freq_REL_fix', 'alt_freq_CR'], alpha=0.8)
    plt.savefig(os.path.join(reporter.out_folder, '{}_parallel_coordinates.png'.format(patient)))


def read_data(patient, reporter, time_points, time_point_normal):
    patient = patient.replace('-', '')

    ID, REL = time_points #may be other time points than ID and REL

    df = reporter.get_standard_df(time_points, time_point_normal, coding_only=True)
    df['mutset'] = '{}_{}-coding'.format(patient, ID)
    df['chrom'] = df.chrom.apply(str)


    #df_all = reporter.get_standard_df(ID, time_point_normal, coding_only=False)
    #df_all['mutset'] = '{}_{}-all mutations'.format(patient, ID)
    #df_all['chrom'] = df_all.chrom.apply(str)

    df['mutations'] = df.mutset.apply(lambda x: x.split('-')[1])
    return df

def cmdline():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', dest='patientid', required=False)
    parser.add_argument('-D', dest='sample_dir', required=True, help='Where the output folder from bcbio is, the root '
                                                                       'for the sample results')
    plots = [PLOT_ALL, PLOT_CALLERS, PLOT_VIOLIN, PLOT_VAF1D, PLOT_VAF2D]
    parser.add_argument('-P', dest='toplot', action='append', choices=plots)
    parser.add_argument('-t', dest='time_points', action='append')
    parser.add_argument('-T', dest='time_point_normal', default='CR')
    options = parser.parse_args()
    print(options.toplot)
    task = ExonMutationPlotter(os.path.abspath(options.sample_dir))
    task.run(os.path.abspath(options.sample_dir), options.patientid, options.toplot, options.time_points, options.time_point_normal)

if __name__ == "__main__":  # detects if called from system
    cmdline()
