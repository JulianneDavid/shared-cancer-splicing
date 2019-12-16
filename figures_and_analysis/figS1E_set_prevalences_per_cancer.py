#!/usr/bin/env python3

"""
fig2_set_prevalences_per_cancer.py

Python 3 code for plotting prevalences per cancer type of tissue-matched,
core normal, and other junctions from TCGA samples.

"""

import argparse
from datetime import datetime
import glob
import logging
from matplotlib import use; use('pdf')
import os
import pandas as pd
import sys
try:
    from utilities.utilities import _TCGA_ABBR, _PER
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _TCGA_ABBR, _PER
from utilities.utilities import  grouped_boxplots_with_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot cancer-type prevalence boxplots per cancer type for '
                    'overall sets.'
    )
    parser.add_argument(
        '--set-memberships', '-s',
        help='Give the directory containing set membership files for all TCGA '
             'cancer types.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output figure and log file'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )

    args = parser.parse_args()
    set_dir = args.set_memberships
    out_path = args.output_path
    log_mode = args.log_level

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'figS2_all_cancer_boxplots_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    data_label_dict = {}
    grouped_data_dict = {}
    for cancer, abbr in _TCGA_ABBR.items():
        per_col = cancer + _PER
        logging.info('starting {}'.format(cancer))
        jx_path = os.path.join(
            set_dir, '{}_piechart_annotation*.csv'.format(cancer)
        )
        jx_file = glob.glob(jx_path)[0]
        max_chunk = 1e6
        jx_df = pd.DataFrame()
        chunks = pd.read_table(jx_file, sep=',', chunksize=max_chunk)
        for chunk in chunks:
            chunk = chunk.fillna(0)
            chunk = chunk[chunk[per_col] >= 0.01]
            jx_df = pd.concat([jx_df, chunk])

        in_paired = jx_df[jx_df['paired'] == 1][per_col].tolist()
        in_gtex = jx_df[
            (jx_df['gtex'] == 1) & (jx_df['paired'] == 0)
        ][per_col].tolist()
        non_gtex = jx_df[jx_df['gtex'] == 0][per_col].tolist()

        grouped_data_dict[abbr] = {}
        grouped_data_dict[abbr]['data'] = [in_paired, in_gtex, non_gtex]
        grouped_data_dict[abbr]['table_data'] = [
            len(in_paired), len(in_gtex), len(non_gtex)
        ]

    plot_info_dict = {}
    plot_info_dict['light colors'] = [
        'xkcd:bright sky blue', 'xkcd:kermit green', 'xkcd:saffron'
    ]
    plot_info_dict['dark colors'] = plot_info_dict['light colors']
    plot_info_dict['legend'] = [
        'tissue-matched\nnormal junctions', 'core normal\njunctions',
        'other\njunctions'
    ]
    plot_info_dict['row font color'] = ['black', 'black', 'black']
    plot_info_dict['row labels'] = plot_info_dict['legend']
    plot_info_dict['row colors'] = plot_info_dict['light colors']

    fig_size = (13.0, 4.0)
    fig_name = 'figS2_all_cancer_boxplots_{}.pdf'.format(now)
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at {}'.format(fig_file))
    grouped_boxplots_with_table(
        grouped_data_dict, plot_info_dict, fig_file, fig_size
    )
