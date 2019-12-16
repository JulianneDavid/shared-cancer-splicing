#!/usr/bin/env python3

"""
fig1A_overall_set_barplots.py

Python 3 code for plotting relative proportions of cancer junctions found in
tissue-matched normal GTEx and TCGA samples, all GTEx and TCGA normal samples
('core normals'), or no GTEx or TCGA normal samples..

"""

import argparse
from datetime import datetime
import glob
import json
import logging
from matplotlib import use; use('pdf')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import os
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
try:
    from utilities.utilities import _TCGA_ABBR
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _TCGA_ABBR

from utilities.utilities import _PER, _CANCER_COLORS, _FONT_COLORS


def barplots_with_table(data_dict, plot_dict, out_path, now, flag=None):
    plt.rcParams.update({'figure.autolayout': True})
    if flag == 'all':
        plt.rcParams['figure.figsize'] = 7.0, 5.0
    elif flag == 'unexpl':
        plt.rcParams['figure.figsize'] = 4.68, 4.0
    else:
        plt.rcParams['figure.figsize'] = 13.0, 4.0

    sns.set_context("paper")
    sns.set_style("whitegrid")
    sns.set(style="ticks")

    bot_data = data_dict['paired']
    mid_data = data_dict['gtex']
    top_data = data_dict['other']
    bot_color = plot_dict['light colors'][0]
    mid_color = plot_dict['light colors'][1]
    top_color = plot_dict['light colors'][2]
    f, (ax, ax2) = plt.subplots(
        nrows=2, ncols=1, gridspec_kw={'height_ratios': [10000, 1]}
    )
    plt.sca(ax)

    num_groups = len(data_dict['paired'])

    ind_l = list(range(1, (3 * num_groups) + 1, 3))
    ind_l = [x - 0.15 for x in ind_l]
    ind_m = [x + 0.85 for x in ind_l]
    ind_r = [x + 0.85 for x in ind_m]

    barwidths = 0.75
    plt.bar(ind_l, bot_data, barwidths, color=bot_color, linewidth=0)
    plt.bar(ind_m, mid_data, barwidths, color=mid_color, linewidth=0)
    plt.bar(ind_r, top_data, barwidths, color=top_color, linewidth=0)

    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    edge_buffer = .7
    ax.set_xlim(left=ind_l[0] - edge_buffer, right=ind_r[-1] + edge_buffer)

    plt.yscale('log')
    plt.ylabel('relative abundance (%)', fontsize=10)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:g}'.format(y))
    )
    plt.setp(
        ax.yaxis.get_majorticklabels(), fontsize=8, color='black'
    )
    columns = data_dict['abbr']

    # Add Table
    rows = ['null']
    whitefont_cols = []
    col_cols = []
    for i, abbr in enumerate(columns):
        try:
            font_color = _FONT_COLORS[abbr]
            col_cols.append('xkcd:{}'.format(_CANCER_COLORS[abbr]))
        except KeyError:
            for can_abbr in _FONT_COLORS.keys():
                if can_abbr in abbr:
                    font_color = _FONT_COLORS[can_abbr]
            for can_abbr in _CANCER_COLORS.keys():
                if can_abbr in abbr:
                    col_cols.append('xkcd:{}'.format(_CANCER_COLORS[can_abbr]))
                    break
        if font_color == 'white':
            whitefont_cols.append(i)

    table_vals = [['' for _ in columns]]
    row_label_cols = ['white']
    the_table = ax.table(
        cellText=table_vals, rowLabels=rows, colLabels=columns, loc='bottom',
        cellLoc='center', colColours=col_cols, rowColours=row_label_cols,
        bbox=[0, -0.09, 1, .075]
    )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5)

    for (row, col), cell in the_table.get_celld().items():
        if (row == 0):
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=7.5)
            )
            if len(columns[0]) > 5:
                cell.set_height(cell.get_height() * 1.5)
            if col in whitefont_cols:
                cell._text.set_color('white')
        else:
            cell.set_height(0)
        if col == -1:
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=6.5)
            )
            cell._text.set_color('white')

    ax2.yaxis.grid(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.patch.set_alpha(0)

    plt.sca(ax2)
    plt.plot([])
    plt.xticks([], [])
    plt.yticks([], [])

    fig = plt.gcf()
    fig_name = 'fig1A_grouped_barplots_{}.pdf'.format(now)
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at {}'.format(fig_file))
    fig.savefig(fig_file)
    plt.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Assigns TCGA cancer junctions to overall groups (tissue-'
                    'matched normals, core normals, or other) and creates '
                    'bar plots comparing groups by cancer type.'
    )
    parser.add_argument(
        '--set-memberships', '-s',
        help='Give the directory containing set membership files for all TCGA '
             'cancer types.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output figure'
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
        out_path, 'fig1A_barplots_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    json_file = os.path.join(out_path, 'barcharts_json_data.txt')

    if os.path.exists(json_file):
        with open(json_file) as recovered_data:
            grouped_data_dict, plot_info_dict = json.load(recovered_data)
    else:
        data_label_dict = {}
        grouped_data_dict = {'paired': [], 'gtex': [], 'other': [], 'abbr': []}
        sort_dict = {}
        for cancer, abbr in _TCGA_ABBR.items():
            grouped_data_dict['abbr'].append(abbr)
            per_col = cancer + _PER
            logging.info('{}'.format(cancer))
            jx_path = os.path.join(
                set_dir, '{}_piechart_annotation*.csv'.format(cancer)
            )
            jx_file = glob.glob(jx_path)[0]
            max_chunk = 1e6
            jx_df = pd.DataFrame()
            chunks = pd.read_table(jx_file, sep=',', chunksize=max_chunk)
            for chunk in chunks:
                chunk = chunk.fillna(0)
                jx_df = pd.concat([jx_df, chunk])

            num_jxs = len(jx_df)
            in_paired = jx_df[jx_df['paired'] == 1][per_col].count()
            in_gtex = jx_df[
                (jx_df['gtex'] == 1) & (jx_df['paired'] == 0)
            ][per_col].count()
            non_gtex = jx_df[jx_df['gtex'] == 0][per_col].count()
            grouped_data_dict['paired'].append(100 * in_paired / num_jxs)
            grouped_data_dict['gtex'].append(100 * in_gtex / num_jxs)
            grouped_data_dict['other'].append(100 * non_gtex / num_jxs)

        grouped_df = pd.DataFrame(grouped_data_dict)
        grouped_df.sort_values(by=['other'], ascending=True, inplace=True)
        grouped_data_dict = grouped_df.to_dict(orient='list')

        plot_info_dict = {}
        plot_info_dict['light colors'] = [
            'xkcd:bright sky blue', 'xkcd:kermit green', 'xkcd:saffron'
        ]

        with open(json_file, 'w') as output:
            json.dump([grouped_data_dict, plot_info_dict], output)

    barplots_with_table(grouped_data_dict, plot_info_dict, out_path, now)
