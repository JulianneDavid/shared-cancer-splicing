#!/usr/bin/env python3

"""
figS1BC_S9BC_junction_sharedness.py

Python 3 code for determining inter- and intra-cancer sharedness.

"""

import argparse
from datetime import datetime
import logging
from math import floor
from matplotlib import use; use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import sys
try:
    from utilities.utilities import _CANCER_COLORS, _TCGA_ABBR, _PER
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _CANCER_COLORS, _TCGA_ABBR, _PER
from utilities.utilities import get_jx_prev_filename, jx_df_from_file


def count_shared_jxs(jx_dir):
    """

    :param jx_dir:
    :return:
    """
    jx_cancer_counts = {}
    intercan_dict = {'jx': [], 'top_x': []}
    sharedness_threshold = 0.1
    over_thresh_dict = {'jx': [], 'top_x': [], 'cancers': []}
    unshared_dict = {'jx': [], 'top_x': []}
    unshared_count = {}
    all_jxs = set()
    intracan_shared_jxs = set()
    intercan_over_thresh = {}
    unshared_jxs = set()
    avg_unshared_total = 0
    nonunique_all_count = 0
    over5_any = set()

    intercan_all = {}
    intercan_all_dict = {'jx': [], 'top_x': []}

    for cancer in _TCGA_ABBR:
        logging.info('starting cancer {}'.format(cancer))
        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            logging.info('no file, continuing')
            continue

        init_df = jx_df_from_file(
            jx_file, min_per=0, max_per=1.0, chunk_it=True,
            glob_form=prev_glob, sample=False, top_x=False
        )

        for jx in init_df['jx']:
            try:
                intercan_all[jx] += 1
            except KeyError:
                intercan_all[jx] = 1

        nonunique_all_count += len(init_df)
        all_jxs.update(init_df.jx)
        per_col = cancer + _PER
        min_sharedness = init_df[per_col].min()
        shared_df = init_df[init_df[per_col] > min_sharedness]
        unshared_df = init_df[init_df[per_col] == min_sharedness]
        intracan_shared_jxs.update(shared_df.jx)
        over5_any.update(shared_df[shared_df[per_col] > 0.05]['jx'])
        unshared_jxs.update(unshared_df.jx)
        avg_unshared_total += len(unshared_df)

        for jx in shared_df['jx']:
            try:
                jx_cancer_counts[jx] += 1
            except KeyError:
                jx_cancer_counts[jx] = 1

        for jx in shared_df[shared_df[per_col] > sharedness_threshold]['jx']:
            percent = (
                '{:,g}'.format(
                    100 * float(shared_df[shared_df['jx'] == jx][per_col])
                )
            )
            try:
                intercan_over_thresh[jx]['top_x'] += 1
                intercan_over_thresh[jx]['cancers'].append(
                    '{}: {}%'.format(cancer, percent)
                )
            except KeyError:
                intercan_over_thresh[jx] = {
                    'top_x': 1, 'cancers': ['{}: {}%'.format(cancer, percent)]
                }

        for jx in unshared_df['jx']:
            try:
                unshared_count[jx] += 1
            except KeyError:
                unshared_count[jx] = 1

    all_unique = len(all_jxs)

    for jx, count in intercan_all.items():
        intercan_all_dict['jx'].append(jx)
        intercan_all_dict['top_x'].append(count)

    intercan_all_df = pd.DataFrame(intercan_all_dict)
    intercan_all_jxs = len(intercan_all_df[intercan_all_df['top_x'] > 1])

    for jx, count in jx_cancer_counts.items():
        intercan_dict['jx'].append(jx)
        intercan_dict['top_x'].append(count)

    intercan_df = pd.DataFrame(intercan_dict)
    intercan_shared_jxs = intercan_df[intercan_df['top_x'] > 1]
    num_intercan_shared_jxs = len(intercan_shared_jxs)

    for jx, count in unshared_count.items():
        unshared_dict['jx'].append(jx)
        unshared_dict['top_x'].append(count)

    unshared_count_df = pd.DataFrame(unshared_dict)
    true_unshared_count = len(
        unshared_count_df[unshared_count_df['top_x'] == 1]
    )

    for jx, data in intercan_over_thresh.items():
        over_thresh_dict['jx'].append(jx)
        over_thresh_dict['top_x'].append(data['top_x'])
        over_thresh_dict['cancers'].append(';\n'.join(data['cancers']))

    over_thresh_df = pd.DataFrame(over_thresh_dict)
    over_thresh_shared_jxs = over_thresh_df[over_thresh_df['top_x'] > 1]

    over_thresh_name = (
        'over_{}%_intercancer_jxs_{}.csv'
        ''.format(sharedness_threshold * 100, flag)
    )
    over_thresh_file = os.path.join(out_path, over_thresh_name)
    with open(over_thresh_file, 'w') as output:
        over_thresh_shared_jxs.to_csv(output, index=False)

    all_shared = len(intracan_shared_jxs)
    over5_shared = len(over5_any)
    logging.info(
        'total number of non-normal tcga cancer junctions: {}'
        ''.format(all_unique)
    )
    logging.info(
        'total number of shared junctions: {}'
        ''.format(all_shared, all_shared / all_unique)
    )
    logging.info(
        'total number of junctions in 5%+ in 1+ cancers: {}, {}'
        ''.format(over5_shared, over5_shared / all_unique)
    )

    logging.info(
        'total number of ANY junctions shared in more than 1 cancer type: {}'
        ''.format(intercan_all_jxs, intercan_all_jxs / all_unique)
    )

    logging.info(
        'total number of inter- and intra-cancer shared junctions: {}'.format(
            num_intercan_shared_jxs, num_intercan_shared_jxs / all_shared
        )
    )
    logging.info(
        'total number of cancer-type specific shared jxs: {}, {}'.format(
            all_shared - len(intercan_shared_jxs),
            (all_shared - len(intercan_shared_jxs)) / all_shared
        )
    )
    logging.info(
        'total number of shared junctions with more than 5% in more than one '
        'cancer type: {}, {}'.format(
            len(over_thresh_shared_jxs),
            len(over_thresh_shared_jxs) / all_shared
        )
    )

    logging.info('unshared jxs:')
    logging.info(
        'total number of unshared non-unique: {}'.format(avg_unshared_total)
    )
    logging.info(
        'avg percent per cancer: {}'
        ''.format(avg_unshared_total / nonunique_all_count)
    )
    logging.info('total number unique: {}'.format(len(unshared_jxs)))
    logging.info(
        'total percent unique: {}'.format(len(unshared_jxs) / all_unique)
    )
    logging.info(
        'total number of junctions in only one sample: {}, {}'
        ''.format(true_unshared_count, true_unshared_count / all_unique)
    )
    return


def prepare_line_data(jx_prev_file, group_col, globform):
    """

    :param jx_prev_file:
    :param group_col:
    :param globform:
    :return:
    """
    init_df = jx_df_from_file(
        jx_prev_file, min_per=0, max_per=1.0, chunk_it=True,
        glob_form=globform, sample=False, top_x=False
    )
    tempgroup = init_df.groupby([group_col])['jx'].count()
    neojx_sharedness_counts = pd.DataFrame(
        {'prev_level': tempgroup.index,
         'neojx_count': tempgroup.values}
    )
    return neojx_sharedness_counts


def prepare_intercancer_data(jx_dir):
    """

    :param jx_dir:
    :return:
    """
    jx_counts = {}
    count_df_dict = {'jx': [], 'top_x': []}
    for cancer in _TCGA_ABBR:
        logging.info('starting cancer {}'.format(cancer))
        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            logging.info('no file, continuing')
            continue

        init_df = jx_df_from_file(
            jx_file, min_per=0, max_per=1.0, chunk_it=True,
            glob_form=prev_glob, sample=False, top_x=False
        )
        jx_list = init_df['jx'].tolist()

        for jx in jx_list:
            try:
                jx_counts[jx] += 1
            except KeyError:
                jx_counts[jx] = 1

    for jx, count in jx_counts.items():
        count_df_dict['jx'].append(jx)
        count_df_dict['top_x'].append(count)

    intercan_df = pd.DataFrame(count_df_dict)

    shared_jxs = intercan_df[intercan_df['top_x'] > 1].sort_values(
        by='top_x', ascending=False
    )
    shared_file = os.path.join(
        out_path, 'intracancer_shared_{}_jxs.csv'.format(flag)
    )
    with open(shared_file, 'w') as output:
        shared_jxs.to_csv(output, index=False)

    tempgroup = intercan_df.groupby(['top_x'])['jx'].count()

    intercan_counts = pd.DataFrame(
        {'num_cancers': tempgroup.index,
         'neojx_count': tempgroup.values}
    )
    return intercan_counts, flag


def intra_cancer_scatterplots(jx_dir):
    """

    :param jx_dir:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 7.0, 2.9
    sns.set_context("paper")
    sns.set_style("whitegrid")

    full_df = pd.DataFrame()
    cancer_colors = []
    for cancer in _TCGA_ABBR:
        logging.info('starting cancer {}'.format(cancer))
        abbr = _TCGA_ABBR[cancer]
        can_color = _CANCER_COLORS[_TCGA_ABBR[cancer]]
        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            continue

        prev_col = cancer + _PER
        can_sharedness = prepare_line_data(jx_file, prev_col, prev_glob)
        can_sharedness['color'] = can_color
        can_sharedness['abbr'] = abbr
        full_df = full_df.append(can_sharedness, sort=True)
        cancer_colors.append(can_color)

    sns.scatterplot(
        x='prev_level', y='neojx_count', hue='abbr', data=full_df,
        legend=False,
        palette=sns.xkcd_palette(cancer_colors),
        edgecolor='None', s=4, alpha=1

    )
    plt.yscale('log')
    return


def inter_cancer_histogram(jx_dir, out_path, now):
    """

    :param jx_dir:
    :param out_path:
    :param now:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 7.0, 2.9
    sns.set_context("paper")
    sns.set_style("whitegrid")

    jx_cantype_counts, flag = prepare_intercancer_data(jx_dir)

    plt.bar(
        jx_cantype_counts['num_cancers'], jx_cantype_counts['neojx_count'],
        width=1,
        edgecolor='black', facecolor='xkcd:saffron'
    )

    ax = plt.gca()
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.yaxis.set_minor_locator(ticker.LogLocator())
    ax.yaxis.set_minor_formatter(ticker.LogFormatter())
    plt.setp(
        ax.yaxis.get_majorticklabels(), fontsize=7, color='black'
    )
    ax.tick_params(which='minor', length=4, color='k')
    plt.setp(
        ax.xaxis.get_majorticklabels(), fontsize=7, color='black'
    )
    plt.xticks(range(1, jx_cantype_counts['num_cancers'].max() + 1))

    plt.xlabel('TCGA cancer types (#)')
    if flag == 'all':
        fig_name = 'figS1C_intercancer_histogram_{}.pdf'.format(now)
        ylabel_text = (
            'junction top_x (#)'
        )
    else:
        fig_name = 'figS9B_intercancer_histogram_{}.pdf'.format(now)
        ylabel_text = (
            'junction top_x (#)\nunexplained stage 3+ junctions'
        )

    plt.ylabel(ylabel_text, fontsize=8)
    plt.yscale('log')
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
    )
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    fig.savefig(fig_file)
    plt.close()
    return


def moving_average(data_array, window):
    moving_array = np.cumsum(data_array)
    moving_array[window:] = moving_array[window:] - moving_array[:-window]
    return moving_array[window - 1:] / window


def intra_cancer_percancer_lineplots(jx_dir, out_path, now):
    """

    :param jx_dir:
    :param out_path:
    :param now:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 7.0, 2.9
    sns.set_context("paper")
    sns.set_style("whitegrid")

    full_df = pd.DataFrame()
    flag = ''
    for cancer in _TCGA_ABBR:
        logging.info('starting cancer {}'.format(cancer))
        abbr = _TCGA_ABBR[cancer]
        can_color = _CANCER_COLORS[_TCGA_ABBR[cancer]]

        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            continue

        prev_col = cancer + _PER
        can_sharedness = prepare_line_data(jx_file, prev_col, prev_glob)
        can_sharedness['color'] = can_color
        can_sharedness['abbr'] = abbr
        num_pts = len(can_sharedness)
        full_df = full_df.append(can_sharedness, sort=True)
        if flag == 'all':
            even_split = np.array_split(can_sharedness, 5)
            mid_df = pd.concat(even_split[1:-1])
            end_split = max(5, floor(len(even_split[-1]) * 0.25))
            end_df = np.array_split(even_split[-1], end_split)
            frontend_df = pd.concat(end_df[:-1])
            split_df = [even_split[0], mid_df, frontend_df, end_df[-1]]
            xvals = []
            yvals = []
            window_lengths = (
                1, max(1, int(floor(len(mid_df) / 4))), 5, 1
            )

        else:
            maj_window_len = max(1, int(floor(len(can_sharedness) / 10)))
            if num_pts > 1:
                num_splits = max(1, floor(num_pts * 0.20))
                even_split = np.array_split(can_sharedness, num_splits)
                try:
                    mid_df = pd.concat(even_split[1:-1])
                    split_df = [even_split[0], mid_df, even_split[-1]]

                except ValueError:
                    split_df = even_split
                xvals = []
                yvals = []
                window_lengths = (1, maj_window_len, 1)
            else:
                split_df = [can_sharedness]
                window_lengths = [1]

        for i, df in enumerate(split_df):
            new_y = moving_average(
                df['neojx_count'].tolist(), window_lengths[i]
            )
            yvals.extend(new_y)
            new_x = moving_average(
                df['prev_level'].tolist(), window_lengths[i]
            )
            xvals.extend(new_x)
        opacity = 0.9
        plt.plot(
            xvals, yvals, '-', color='xkcd:{}'.format(can_color), alpha=opacity
        )

    ax = plt.gca()
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:g}'.format(y))
    )
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x, _: '{}%'.format('{:,g}'.format(100 * x))
        )
    )
    plt.xlabel('TCGA cancer-type jx prevalence', fontsize=8)
    plt.yscale('log')
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
    )
    if flag == 'all':
        ylabel_text = 'junction top_x (#)'
        fig_name = (
            'figS1B_per_cancer_line_and_scatter_plot_{}.pdf'.format(now)
        )
    elif flag == 'unexpl':
        ylabel_text = 'junction top_x (#)\nunexplained level 3+ junctions'
        fig_name = (
            'figS9C_per_cancer_line_and_scatter_plot_{}.pdf'.format(now)
        )
    else:
        ylabel_text = 'junction top_x (#)\ncheck filter'
        fig_name = (
            'per_cancer_line_and_scatter_plot_{}.pdf'.format(now)
        )

    plt.ylabel(ylabel_text, fontsize=8)
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving fig at {}'.format(fig_file))
    fig.savefig(fig_file)
    plt.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate and plot inter- and intra-cancer type junction '
                    'recurrence.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output files and figures..'
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
        out_path, 'junction_sharedness_analysis_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('command line: {}'.format(' '.join(sys.argv)))

    inter_cancer_histogram(jx_dir, out_path, now)

    intra_cancer_scatterplots(jx_dir)
    intra_cancer_percancer_lineplots(jx_dir, out_path, now)

    count_shared_jxs(jx_dir)
