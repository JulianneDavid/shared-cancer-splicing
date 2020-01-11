#!/usr/bin/env python3

"""
fig2A_TCGA_heatmap.py

Python 3 code for clustering TCGA-found junctions by prevalence in all TCGA
cancer types.

Note "true" prevalence files are not required since these use the top 200
shared junctions only, so all have at least 2 reads per junction.

"""

import argparse
from datetime import datetime
import logging
from matplotlib import use; use('pdf')
import numpy as np
import os
import seaborn as sns; sns.set(color_codes=True)
import sys
try:
    from utilities.utilities import  masked_double_heatmap, jx_dir_to_df
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import masked_double_heatmap, jx_dir_to_df

from utilities.utilities import _CANCER_COLORS,_ALL_ABBR, _TCGA_CANCER_TYPES
from utilities.utilities import mm2inch

def cluster_tcga_jxs(db_jx_df, out_path, now):
    """

    :param db_jx_df:
    :param out_path:
    :param now:
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
    db_jx_df = db_jx_df[(db_jx_df!=0).any(axis=1)]
    full_df = db_jx_df.fillna(0)

    metric = 'seuclidean'
    method = 'ward'
    xlabelset = [
        'ESCA', 'STAD', 'OV', 'BLCA', 'LUSC', 'HNSC', 'CESC', 'THYM',
        'DLBC', 'ACC', 'PCPG', 'PAAD', 'LUAD', 'THCA', 'MESO', 'SARC',
        'UCS', 'UCEC', 'BRCA', 'TGCT', 'PRAD', 'LIHC', 'CHOL', 'KIRC',
        'KIRP', 'KICH', 'LAML', 'GBM', 'LGG', 'READ', 'COAD', 'UVM',
        'SKCM'
    ]

    full_cluster = sns.clustermap(
        full_df.drop('jx', axis=1), method=method, metric=metric,
        standard_scale=1, col_cluster=False
    )

    row_order = full_cluster.dendrogram_row.reordered_ind
    tcga_cluster = full_df.reindex(row_order, axis='rows')
    tcga_cluster = tcga_cluster.reset_index(drop=True)

    full_merge = tcga_cluster.drop('jx', axis=1)
    merge_df = full_merge[xlabelset]

    logging.info('total junction top_x is: {}'.format(len(merge_df)))

    fig_name = ('fig2A_TCGA_heatmap_{}.pdf'.format(now))
    fig_file = os.path.join(out_path, fig_name)
    label_ranges = []
    df_cols = merge_df.columns.values.tolist()
    for i in range(len(df_cols)):
        label_ranges.append([i, i+1])
    label_ranges = np.array(label_ranges)
    colors = [
        'xkcd:{}'.format(_CANCER_COLORS[x]) for x in df_cols
    ]
    texts = None
    fontcolors = None

    colorbar_dict = {
        'ranges': label_ranges, 'colors': colors, 'labels': texts,
        'fontcolors': fontcolors, 'height_percent': '2.0%'
    }

    fig_size = mm2inch(89, 84)
    masked_double_heatmap(
        merge_df, merge_df.columns.values.tolist(), fig_file, size=fig_size,
        colorbar=colorbar_dict, cbarticks = [0.001, .01, .1, .8],
        cbar_font_adder=2, xlabel_fontsize=4, bottom_pad=-3.5
    )
    logging.info('saving figure at {}'.format(fig_file))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cluster TCGA junctions by prevalence across cancer types.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output figure and log file.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--database-junction_directory', '-d',
        help='specify a directory containing .csv files with junctions '
             'extracted via a jx_indexer query, each containing prevalence '
             'values for one cancer type. Please run with direct output from '
             'jx_indexer query, not results of set membership analysis.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    log_mode = args.log_level
    jx_dir = args.database_junction_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'fig2A_TCGA_heatmap_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    count = 200
    min_per = 0.01
    max_per = 1.0
    chunk_it = True
    order = _TCGA_CANCER_TYPES
    tcga_jx_df = jx_dir_to_df(
        jx_dir, min_per, max_per, chunk_it, top_x=count, percentile=1,
        can_order=order
    )
    cluster_tcga_jxs(tcga_jx_df, out_path, now)
