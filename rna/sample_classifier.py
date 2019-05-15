import argparse
import glob
import logging
import os
from os.path import join
from sklearn import preprocessing
from allpy.SimpleTask import SampleTask
from allpy.config import pipeline_config
from allpy.config.config_constants import NOTIFY_EMAIL


from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE


import numpy as np
from sklearn.calibration import CalibratedClassifierCV

import seaborn as sns

from allpy.notify.notify import send_email

CR = "CR"

logger = logging.getLogger()


class DataProcesser(object):
    def __init__(self):
        self.pca = None


    def rank_scale_reduce(self, df):
        df = self.rank(df)
        df = self.scale(df)
        df = self.reduce(df)
        return df

    def scale(self, df):
        return pd.DataFrame(preprocessing.minmax_scale(df),
                           index=df.index,
                           columns=df.columns)
    def rank(self, df):
        return df.T.rank().T

    def reduce(self, df):
        i = df.index

        if self.pca is None:
            logger.info("Fitting training data to 100 PCA components")
            self.pca = PCA(n_components=100,random_state=44)
            self.pca.fit(df)

        df = self.pca.transform(df)
        return pd.DataFrame(df, index=i)

    def joint_tsne(self, df):
        perplexities = [10, 30]
        learning_rates = [50,100,200]

        tsne_outputs = []

        logging.info("calculating TSNEs")

        for perplexity in perplexities:
            for learning_rate in learning_rates:
                tsne = TSNE(n_components=2, random_state=40,perplexity=perplexity,learning_rate=learning_rate)
                y = tsne.fit_transform(df)
                tmp_tsne = pd.DataFrame(y,index=df.index,columns=['x','y'])\
                                .assign(perplexity=perplexity,learning_rate=learning_rate)\
                                .assign(tsne_name="p{}l{}".format(int(perplexity),int(learning_rate)))
                tsne_outputs.append(tmp_tsne)


        tsne_df = pd.concat(tsne_outputs)

        tsne_pivot = tsne_df.reset_index()\
                .pivot_table(index='index',columns=['learning_rate','perplexity'],values=['x','y'])
        return tsne_df, tsne_pivot



class RNASampleClassifier(SampleTask):
    def __init__(self):
        super().__init__("SampleClassifier")
        self.data_processor = DataProcesser()

        lightgrey = '#D3D3D3'
        self.color_group_map = {
             'DUX4': '#299bd8',
             'ETV6--PAX5': '#fa9e1c',
             'ETV6--RUNX1': '#0a9f87',
             #'BCL2f': '#ff00e9',
             'BCL2f': '#a7f0ff',
             'KMT2Af': '#00DBFF',
             'MEF2Df': '#0650de',
             'NH-HeH': '#f08a8a',
             'LH-NT': '#c45923',
             'PAX5-plus': '#c7ff00',
             'PBX1--TCF3': '#7f390f',
             'Ph-like': '#d6e208',
             'Ph-like;Ph-pos': '#d6e208',
             'Ph-pos': '#2d9302',
             'CEBPEf': '#00008B',
             'ZNF384f': '#6e3f5f',
             'HOXA':'#ffa500',
             'TAL/LMO': '#fdd9ba',
             'TLX1': '#ffbeec',
             'immature': '#b3acff',
             #'TLX3':'#a7f0ff',
             'TLX3':'#ff00e9',
             'unknown': lightgrey,
             'unassigned-B': lightgrey,
             'unassigned-T': lightgrey
        }

    def derive_sample_name(self, task_input, cmd_sample_name):
        if cmd_sample_name is not None:
            patient = cmd_sample_name
        else:
            patient = os.path.basename(os.path.abspath(task_input))
        return patient

    def load_configuration(self):
        config = pipeline_config.Config()
        conf_dict = config.get_config_dict([NOTIFY_EMAIL])

        return conf_dict

    def _load_training_data(self):

        self.subgroup_map = {"Ph-pos": 'Ph-like;Ph-pos', "Ph-like": 'Ph-like;Ph-pos', 'unassigned':'unassigned-B'}
        samples_groupings = pd.read_table('/home/mpschr/mounts/all_data/sample_info/sample_groupings.txt',
                                               index_col=0)\
                                    .query('entity==["BCP-ALL","T-ALL"] and sample_name != "{}"'.format(self.sample))\
                                    .replace(self.subgroup_map)

        tall_subgroups = pd.read_table('/home/mpschr/Documents/projects/martin/tall-classifier/tall-classes.martin.txt',index_col=0)\
                            .replace({'other': 'unassigned-T'})
        tall_subgroups.index = tall_subgroups.index.rename('sample_name')
        tall_subgroups.head()
        samples_groupings.subgroup.update(tall_subgroups.Homminga.replace({'other':'unassigned-T'}))

        self.all_samples = samples_groupings
        self.X = self._load_dkh_rna_data(samples_groupings)
        dkh_blacklist = self._get_dkh_blacklist(samples_groupings)
        train_and_pred_samples = self.X.index.difference(dkh_blacklist)
        self.training_groups = samples_groupings.loc[train_and_pred_samples].subgroup.value_counts()
        self.X, self.prot_coding = self._preprocess_dkh_rna_data(samples_groupings, train_and_pred_samples)
        self.y = samples_groupings.ix[train_and_pred_samples].subgroup.dropna()


        self.ball_or_tall = samples_groupings[['subgroup','entity']].drop_duplicates().dropna().set_index('subgroup').entity.to_dict()
        
        self.selected_genes = self._load_classifier_genes()

        self.X = self.X.ix[self.y.index, self.selected_genes]

        logger.info("DKH Training samples: {} / Groups:\n{}".format(self.y.shape[0], self.y.fillna('-').value_counts()))
        logger.info("DKH Training set dimensions (samples/genes): {}".format(self.X.shape))

        return

    def _load_classifier_genes(self):
        subgroup_anovas_file = '/home/mpschr/Documents/projects/lorenz/rna_seq_cluster_medians/subgroup_anovas.20180724.tsv'
        subgroup_anovas = pd.read_table(subgroup_anovas_file,index_col=['gene_name', 'gene_id'])

        logger.info("Subgroup-Anovas data set {}".format(subgroup_anovas_file.split('.')[-2]))

        subgroup_anovas = subgroup_anovas.reset_index().ix[
            subgroup_anovas.reset_index().gene_name.isin(self.prot_coding.index)]
        subgroup_genes_150_selection = subgroup_anovas.loc[lambda df: df.gene_name.str.contains('-') == False] \
            .query('n_group > 3 and median_group > -1 and anova_q < 0.05') \
            .sort_values(['group', 'log2_fold_change'], ascending=[True, False]).groupby('group').head(50)
        subgroup_genes_150_selection = subgroup_genes_150_selection.append(
            subgroup_anovas.query("group == 'NH-HeH' and log2_fold_change < -4  and median_others > 1"))
        subgroup_genes_150_selection = subgroup_genes_150_selection.append(
            subgroup_anovas.query("group == 'Ph-like;Ph-pos' and log2_fold_change < -4  and median_others > 1"))
        subgroup_genes_150 = subgroup_genes_150_selection.gene_name.drop_duplicates()
        logger.info("Selected genes for classifier: {} / Subgroup specificities:\n{}".format(
            subgroup_genes_150.shape[0],
            subgroup_genes_150_selection.group.value_counts()
        ))
        return subgroup_genes_150

    def _get_dkh_blacklist(self, ball_samples_groupings):
        blacklist_bad_quality = ball_samples_groupings.query('data_rna != "complete"').index.tolist()
        blacklist_REL = ball_samples_groupings.query('status == "REL"').index.tolist()
        blacklist_Xeno = ['Berlin2265', 'Berlin2495', 'BerlinML453', 'Berlin2493', 'Berlin2494', 'Berlin2474']
        dkh_blacklist = blacklist_bad_quality + blacklist_REL + blacklist_Xeno
        return dkh_blacklist

    def _preprocess_dkh_rna_data(self, ball_samples_groupings, train_and_pred_samples):
        non_dup = self.X.columns.value_counts().loc[lambda x: x == 1].index
        non_null = (pd.isnull(self.X.ix[train_and_pred_samples, non_dup]).sum() < self.X.shape[0] * 0.15).loc[
            lambda x: x == True]
        protein_coding_genes = pd.read_table('/home/mpschr/Data/ensembl75.genesym-genebiotype.txt') \
            .rename(columns=lambda x: x.replace(' ', '_')).query("Gene_type == 'protein_coding'")
        prot_coding = non_null.loc[protein_coding_genes.Associated_Gene_Name].dropna()
        absolute_min = round(self.X.min().min()) - 1

        return self.X[prot_coding.index].fillna(absolute_min).loc[train_and_pred_samples], prot_coding

    def _load_dkh_rna_data(self, ball_samples_groupings):
        rna_data_file = glob.glob('/home/mpschr/Documents/projects/rna/output/cohort/tpm_log2.matrix.ALL.*.tsv')[-1]
        rna_data_file_version = rna_data_file.split('.')[-2]
        logger.info("AG Baldus RNA data set: {}".format(rna_data_file_version))
        combined_tpmlog2_matrix = pd.read_table(rna_data_file, index_col=0)
        combined_tpmlog2_matrix.head()
        return combined_tpmlog2_matrix.set_index('gene_name').transpose()



    def run(self, sample_dir, exec_dir='.', cmd_sample_name=None, patient_id=None):
        conf_dict = self.load_configuration()
        super().run(sample_dir,exec_dir,cmd_sample_name,patient_id)

        self.output_dir = join(sample_dir,'sample-classification')
        os.makedirs(self.output_dir,exist_ok=True)


        fh = logging.FileHandler(join(self.output_dir,'classification.log'))
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

        #self.output_dir = '/home/mpschr/Documents/projects/rna/subtype-classifer/devtests'

        self._load_training_data()
        self._sample_topredict = self._load_prediction_sample(sample_dir)

        #X_reduced_training = self.data_processor.rank_scale_reduce(self.X)
        #X_reduced_predict = self.data_processor.rank_scale_reduce(self._sample_topredict[self.X.columns])
        #X_reduced = X_reduced_training.append(X_reduced_predict)

        X_reduced = self.data_processor.rank_scale_reduce(self.X.append(self._sample_topredict[self.X.columns]))

        tsne_plotdata, X_tsne = self.data_processor.joint_tsne(X_reduced)
        self._plot_tsne(tsne_plotdata)

        self.X_training = X_reduced.join(X_tsne, how='left').ix[self.y.index]
        self.X_predict = X_reduced.join(X_tsne, how='left').ix[[self.sample]]


        logger.info("training data shape: {}".format(self.X_training.shape))
        logger.info("prediction data shape: {}".format(self.X_predict.shape))

        #self._train()

        prediction = self._train_predict()

        prediction.T.to_csv(join(self.output_dir,'{}.prediction.txt'.format(self.sample)),sep='\t')
        X_reduced.join(X_tsne,how='left').to_csv(join(self.output_dir,'{}.training-data.txt'.format(self.sample)),sep='\t')

        self._plot_voting(prediction)
        self._plot_heatmap(X_tsne, 'TSNE')
        self._plot_heatmap(X_reduced, 'PCA', 1)
        #self._plot_heatmap(self.X.append(self._sample_topredict[self.X.columns]), 'genes', zscore=1)

        elapsed_time = self.timer.get_elapsed_time()
        logging.info(elapsed_time)

        return


        # if max(returncodes.values()) == 0 and min(returncodes.values()) == 0:
        #     send_email('michael.p.schroeder@gmail.com', "{} RNA-SEQ-C done".format(reportid), elapsed_time)
        # else:
        #     send_email('michael.p.schroeder@gmail.com', "FAILURE of {}  RNA-SEQ-C".format(reportid),
        #                "{}\n\n{}".format(elapsed_time, returncodes))

    def _plot_heatmap(self, data, dataname,standard_scale=None, zscore=None, cbar_kws={}):

        logger.info('Plotting clustered heatmap')

        data = data.T

        colormap=self.y.replace(self.color_group_map).to_dict()
        colormap.update({self.sample: 'black'})


        colors = [colormap[x] if x in colormap else 'white' for x in data.columns]

        data.columns = [self.sample if x == self.sample else '' for x in data.columns]

        p = sns.clustermap(data,
                           standard_scale=standard_scale,
                           z_score=zscore,
                           col_colors=colors,
                           cbar_kws=cbar_kws,


        )
        #sns.heatmap()
        p.fig.savefig(
            join(self.output_dir, '{}.heatmap.{}.png'.format(self.sample,dataname))
        )

    def _plot_tsne(self, tsne_data):

        logger.info('Plotting TSNE')

        plotdf = tsne_data.query("learning_rate == 100 and perplexity == 30")
        plotdf = plotdf.join(self.all_samples[['subgroup']],how='left')
        plotdf['subgroup'] = plotdf.subgroup.fillna(self.sample)

        colormap=self.color_group_map.copy()
        colormap.update({self.sample: 'black'})

        grouporder = plotdf.subgroup.fillna('unassigned').sort_values().drop_duplicates().tolist()
        
        sns.set_style("whitegrid")
        fig, ax = plt.subplots(figsize=(12,10))
        ax.margins(0.05)
        fig.tight_layout()

        markers = {}
        for g in grouporder:
            if g == self.sample:
                markers[g] = "*" 
            elif self.ball_or_tall[g] == 'BCP-ALL':
                  markers[g] =  'o'
            elif self.ball_or_tall[g] == 'T-ALL':
                  markers[g] = 's'

        for group in grouporder:
            marker = '*' if  group == self.sample else 'o'
            color = colormap[group]
            ax.scatter(plotdf.query('subgroup == "{}"'.format(group)).x, plotdf.query('subgroup == "{}"'.format(group)).y,
                       marker=markers[group],
                       label=group,
                       color=color,
                       s=50 if group != self.sample else 100,
                       alpha=0.8)
        ax.legend(bbox_to_anchor=(1.15, 1.))

        fig.savefig(
            join(self.output_dir, '{}.tsne.svg'.format(self.sample)),
            bbox_inches='tight'
        )
        fig.savefig(
            join(self.output_dir, '{}.tsne.png'.format(self.sample)),
            bbox_inches='tight'
        )

    def _plot_voting(self, prediction):

        logger.info('Creating classifier voting plot')

        dkh_validation_res = pd.read_table(
            '/home/mpschr/Documents/projects/rna/subtype-classifer/dkh_classification_results.tsv', index_col=0)
        dkh_validation_res = dkh_validation_res.append(prediction.assign(subgroup='topredict').\
                                                       rename(columns={'highest_proba':'pred'}))
        dkh_validation_res.index = dkh_validation_res.index.rename('sample_name')

        plot_data = pd.melt(dkh_validation_res.drop(['classification'],axis=1).reset_index(), id_vars=['sample_name', 'subgroup', 'OK', 'pred'],
                         var_name='pred_clf', value_name='proba') \
                .assign(From_group=lambda df:
                            df.apply(lambda r: str(r.pred_clf == r.subgroup) if r.subgroup != 'topredict' else self.sample,
                            axis=1))

        p = sns.factorplot(
            data=plot_data.sort_values('From_group'),
            x='proba', y='pred_clf', hue='From_group', alpha=0.7,
            kind='strip', jitter=True, hue_order=['False','True',self.sample],
            palette=['blue','green','black'],
            size=10
        )
        p.fig.savefig(
            join(self.output_dir, '{}.subgroup-votes.svg'.format(self.sample))
        )

    def _train_predict(self):
        rf_params = {'criterion': 'entropy', 'min_samples_split': 2, 'n_estimators': 1000}
        logger.info("training random forest classifier with {}".format(rf_params))
        rf_clf = RandomForestClassifier(random_state=44, **rf_params)
        rf_clf = CalibratedClassifierCV(rf_clf, cv=None)

        rf_clf = rf_clf.fit(self.X_training.ix[self.y.index], self.y)

        res_rf = (
            pd.Series(rf_clf.predict(self.X_predict), index=self.X_predict.index).rename('highest_proba').to_frame()
                .join(
                pd.DataFrame(rf_clf.predict_proba(self.X_predict), index=self.X_predict.index, columns=rf_clf.classes_)
                )
        )
        res_rf.insert(0,'classification',
                      res_rf.apply(lambda r: r.highest_proba if r.ix[rf_clf.classes_].max() > 0.35 else 'unassigned',axis=1))

        res_rf.index = res_rf.index.rename('sample_name')
        pred = res_rf.classification.ix[self.sample]
        logging.info("CLASSIFICATION of {}: {}".format(self.sample, pred))
        return res_rf

    def _load_prediction_sample(self, sample_dir):

        sample_expression_file = join(sample_dir, 'stringtie', '{}.tsv'.format(self.sample))
        logger.info("loading expression from file {}".format(sample_expression_file))
        sample_df = pd.read_table(sample_expression_file).set_index('Gene Name')[['TPM']]
        sample_df = sample_df.ix[self.selected_genes].sort_values('TPM',ascending=False)\
                        .reset_index().drop_duplicates(subset='Gene Name').set_index('Gene Name')\
                        .rename(columns={"TPM":self.sample}).T
        absolute_min = round(sample_df.dropna().min().min()) - 1

        sample_df = sample_df.fillna(absolute_min)

        return  sample_df
        #sample_df_reduced = self.data_processor.rank_scale_reduce(sample_df)

        #return sample_df_reduced




def absolute_path(path):
    if path is None:
        return None

    return os.path.abspath(path)


def cmdline_rna():
    cmdline(seqtype='rna')


def cmdline(seqtype='exon'):
    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-S', dest='sample_dir', default=None,
                       help="Sample folder, containing the 'stringtie' folder with the .tsv expression file")

    parser.add_argument('-p', dest='patientid', required=False)

    options = parser.parse_args()
    print(options)

    task = RNASampleClassifier()
    task.run(absolute_path(options.sample_dir), options.patientid)


if __name__ == "__main__":  # detects if called from system
    cmdline()
