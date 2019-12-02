#!/usr/bin/env python3

"""
figS6_full_TCGA_SRA_heatmap.py

Python 3 code for clustering TCGA junctions by prevalence across all TCGA
cancer types and subtypes and all non-cancer SRA samples.

"""

import argparse
from datetime import datetime
import glob
import logging
from math import floor
from matplotlib import use; use('pdf')
import numpy as np
import os
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import sys
try:
    from utilities.utilities import _FULL_HEATMAP_COLORS, collect_mutual_jxs
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _FULL_HEATMAP_COLORS, collect_mutual_jxs
from utilities.utilities import _ALL_ABBR, masked_double_heatmap, _SRA_ABBR
from utilities.utilities import collect_metasra_count, snaptron_results_to_jxs
from utilities.utilities import jx_dir_to_df


def cluster_jxs_together(db_jx_df, expt_jx_dict, out_path, count_dict,
                         sra_threshold, now):
    """

    :param db_jx_df:
    :param expt_jx_dict:
    :param mode:
    :param bin_sizes:
    :param bin_max:
    :param cancer:
    :return:
    """
    tcga_cols = db_jx_df.drop(['jx'], axis=1).columns.values.tolist()
    renames = {}
    xlabels_init = []
    for name in tcga_cols:
        can = name.split('_Sample_Percents')[0]
        abbr = _ALL_ABBR[can]
        xlabels_init.append(abbr)
        renames[name] = abbr
    db_jx_df.rename(renames, axis=1, inplace=True)
    short_cancers = ['DT', 'CASC', 'MPNT', 'SYNS']
    db_jx_df.drop(short_cancers, axis=1, inplace=True)
    db_jx_df = db_jx_df[(db_jx_df!=0).any(axis=1)]

    all_sra_jxs = set()
    for key in expt_jx_dict:
        expt_jxs = expt_jx_dict[key]
        jxs = collect_mutual_jxs(db_jx_df, expt_jxs, key, ret_type='jxs')
        all_sra_jxs = all_sra_jxs.union(jxs)

    full_df = db_jx_df.fillna(0)
    full_df = full_df[full_df['jx'].isin(all_sra_jxs)]
    full_df = full_df.reset_index(drop=True)
    tcga_master_cols = full_df.drop('jx', axis=1).columns.values.tolist()

    sra_cols = list(expt_jx_dict.keys())
    sra_abbr_cols = []
    for item in sra_cols:
        try:
            sra_abbr_cols.append(_SRA_ABBR[item])
        except KeyError:
            continue

    sra_master = pd.DataFrame(full_df['jx'])
    for key in expt_jx_dict:
        try:
            name = _SRA_ABBR[key]
        except KeyError:
            continue
        sra_master[name] = sra_master.jx.apply(
            lambda x:
            expt_jx_dict[key].get(x, 0) / count_dict[key]
        )

    null_thresh = floor(sra_threshold * len(sra_master))

    sra_master = sra_master.replace(to_replace=0, value=np.nan)
    sra_master = sra_master.dropna(thresh=null_thresh, axis=1).fillna(0)

    metric = 'seuclidean'
    method = 'ward'

    xlabelset = [
        'THCA', 'MESO', 'LUAD', 'mesenchymal stem cells',
        'COAD', 'READ', 'PAAD', 'pancreatic islet primary cell',
        'UCEC', 'ECAD', 'BRCA', 'ACC', 'placenta tissue',
        'embryo tissue', 'UCS', 'LIHC', 'hepatocyte primary cell',
        'epithelial primary cell', 'CHOL', 'THYM',
        'thymus primary cell', 'thymus tissue',
        'hematopoietic cell line', 'leukocyte cell line',
        'lymphocyte cell line', 'lymphocyte primary cell', 'DLBC',
        'PRAD', 'KIRC', 'KIRP', 'KICH',
        'somatic stem cells', 'TGCT', 'iPS cell line',
        'pluripotent stem cell line', 'embryo stem cells',
        'zygote primary cell', 'embryo primary cell', 'LAML',
        'UVM', 'SKCM', 'melanocyte cell line', 'neonate cellline',
        'SARC', 'LMS', 'UPLS', 'MFS', 'PCPG', 'PGG', 'PCHC', 'GBM',
        'LGG', 'OAC', 'ODG', 'AC', 'CESC', 'CSC', 'HNSC', 'BLCA',
        'LUSC', 'STAD', 'ESAD', 'epithelial cell line', 'ESCA',
        'ESSC', 'embryo cell line', 'OV'
    ]
    full_cluster = sns.clustermap(
        full_df.drop('jx', axis=1), method=method, metric=metric,
        standard_scale=1, col_cluster=False
    )

    row_order = full_cluster.dendrogram_row.reordered_ind
    tcga_cluster = full_df.reindex(row_order, axis='rows')
    tcga_cluster = tcga_cluster.reset_index(drop=True)

    full_merge = pd.merge(tcga_cluster, sra_master, on=['jx'], how='inner')
    full_merge = full_merge.drop('jx', axis=1)
    merge_df = full_merge[xlabelset]
    logging.info('total junction top_x is: {}'.format(len(merge_df)))
    fig_name = ('figS2A_full_TCGA_SRA_heatmap_{}.pdf'.format(now))
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at: {}'.format(fig_file))

    label_ranges = []
    df_cols = merge_df.columns.values.tolist()
    for i in range(len(df_cols)):
        label_ranges.append([i, i+1])
    label_ranges = np.array(label_ranges)
    colors = [
        'xkcd:{}'.format(_FULL_HEATMAP_COLORS[x]) for x in df_cols
    ]
    texts = None
    fontcolors = None

    colorbar_dict = {
        'ranges': label_ranges, 'colors': colors, 'labels': texts,
        'fontcolors': fontcolors, 'height_percent': '1.75%'
    }
    masked_double_heatmap(
        merge_df, tcga_master_cols, fig_file, colorbar=colorbar_dict,
        other_cbar_label='SRA cell type prevalence', size=(5, 4),
        bottom_pad=-3.5, xlabel_fontsize=4,
        cbar_labels=['0.1%', '1%', '10%', '100%']
    )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cluster TCGA junctions by prevalence across all TCGA '
                    'cancer types and subtypes and all non-cancer SRA samples.'
    )
    parser.add_argument(
        '--snaptron-results',
        help='directory containing multiple .txt files containing results '
             'from a previous snaptron search.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output figure and log file'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--database-junction_directory', '-d',
        help='specify a directory containing .csv files with junctions '
             'extracted via a jx_indexer query, each containing prevalence '
             'values for one cancer type.'
    )
    parser.add_argument(
        '--exptlist-directory', '-e',
        help='If snaptron junctions have already been collected, provide the '
             'directory where the lists of recount/snaptron-available SRA '
             'experiments for each cell type are stored.'
    )

    args = parser.parse_args()
    snap_dir = args.snaptron_results
    out_path = args.output_path
    log_mode = args.log_level
    jx_dir = args.database_junction_directory
    exptlist_dir = args.exptlist_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'double_heatmap_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    # Load TCGA junction dataframe
    top_x = 200
    min_per = 0.01
    max_per = 1.0
    chunk_it = True
    jx_df = jx_dir_to_df(
        jx_dir, min_per, max_per, chunk_it, top_x=top_x, percentile=1
    )

    # Load SRA junctions
    meta_sra_expt_counts = {}
    snaptron_jxs = {}
    sra_threshold = 0.01
    snap_glob = os.path.join(snap_dir, '*rawresults*.txt')
    snap_files = glob.glob(snap_glob)
    for snap_file in snap_files:
        name_tag = os.path.basename(snap_file).split('.')[0]
        name_tag = name_tag.split('_rawresults')[0]
        try:
            name_tag = name_tag.split('metaSRA-runs_')[1]
        except IndexError:
            pass

        meta_sra_expt_counts[name_tag] = collect_metasra_count(
            name_tag, exptlist_dir
        )
        if meta_sra_expt_counts[name_tag] < 20:
            continue

        logging.info('{}: {}'.format(name_tag, meta_sra_expt_counts[name_tag]))
        with open(snap_file) as lines:
            snaptron_jxs[name_tag] = snaptron_results_to_jxs(lines)

    # Cluster SRA and TCGA junctions and plot heatmap
    cluster_jxs_together(
        jx_df, snaptron_jxs, out_path, meta_sra_expt_counts, sra_threshold,
        now
    )
