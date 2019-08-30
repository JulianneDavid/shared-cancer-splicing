#!/usr/bin/env python3

"""
fig2B_TCGA_subtype_heatmaps.py
Python 3 code for clustering TCGA-found junctions by prevalence in some TCGA
cancer types and their histological subtypes.

"""

import argparse
from datetime import datetime
import logging
from matplotlib import use; use('pdf')
import matplotlib.cm as cm
import numpy as np
import os
import seaborn as sns; sns.set(color_codes=True)
import sys
try:
    from utilities.utilities import jx_dir_to_df, masked_double_heatmap
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import jx_dir_to_df, masked_double_heatmap
from utilities.utilities import _ALL_ABBR, _CANCER_COLORS, _FONT_COLORS


def cluster_subtype_jxs(db_jx_df, out_path, now):
    """

    :param db_jx_df:
    :param expt_jx_dict:
    :param mode:
    :param bin_sizes:
    :param bin_max:
    :param cancer:
    :return:
    """
    tcga_master_cols = ['CESC', 'PCPG', 'SARC', 'LGG', 'ESCA']
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

    metric = 'cityblock'
    method = 'ward'
    xlabelset = [
        'CESC', 'CSC', 'ECAD', 'CASC',
        'PCPG', 'PCHC', 'PGG',
        'SARC', 'LMS', 'UPLS', 'DT', 'SYNS', 'MFS', 'MPNT',
        'LGG', 'AC', 'OAC', 'ODG',
        'ESCA', 'ESSC', 'ESAD'
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

    fig_name = ('fig2B_TCGA_subtype_heatmap_{}.pdf'.format(now))
    fig_file = os.path.join(out_path, fig_name)
    label_ranges = np.array(
        [[0, 4], [4, 7], [7, 14], [14, 18], [18, 21]]
    )
    texts = ['CESC', 'PCPG', 'SARC', 'LGG', 'ESCA',]
    colors = [
        'xkcd:{}'.format(_CANCER_COLORS[x]) for x in texts
    ]
    fontcolors = {
        text: _FONT_COLORS[text] for text in texts
    }
    colorbar_dict = {
        'ranges': label_ranges, 'colors': colors, 'labels': texts,
        'fontcolors': fontcolors
    }
    vline_pos = [4, 7, 14, 18]
    masked_double_heatmap(
        merge_df, tcga_master_cols, fig_file, masked_cmap=cm.Greys,
        other_cmap=cm.Blues, colorbar=colorbar_dict,
        vline_pos=vline_pos,
        other_cbar_label='cancer subtype prevalence'
    )
    logging.info('saving figure at: {}'.format(fig_file))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cluster TCGA junctions by prevalence across cancer types '
                    'and their subtypes.'
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
             'values for one cancer type.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    log_mode = args.log_level
    jx_dir = args.database_junction_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'fig2B_TCGA_subtype_heatmap_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    count = 200
    min_per = 0.01
    max_per = 1.0
    chunk_it = True

    subtype_order = [
        'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma',
        'Cervical_Adenosquamous', 'Cervical_Squamous_Cell_Carcinoma',
        'Endocervical_Adenocarcinoma',
        'Brain_Lower_Grade_Glioma',
        'Astrocytoma', 'Oligoastrocytoma', 'Oligodendroglioma',
        'Pheochromocytoma_and_Paraganglioma',
        'Paraganglioma', 'Pheochromocytoma',
        'Sarcoma',
        'Desmoid_Tumor', 'Leiomyosarcoma',
        'Malignant_Peripheral_Nerve_Sheath_Tumors', 'Myxofibrosarcoma',
        'Synovial_Sarcoma', 'Undifferentiated_Pleomorphic_Sarcoma',
        'Esophageal_Carcinoma',
        'Esophagus_Adenocarcinoma', 'Esophagus_Squamous_Cell_Carcinoma'
    ]
    subtype_jx_df = jx_dir_to_df(
        jx_dir, min_per, max_per, chunk_it, top_x=count, percentile=1,
        can_order=subtype_order
    )
    cluster_subtype_jxs(subtype_jx_df, out_path, now)
