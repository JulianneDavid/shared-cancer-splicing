#!/usr/bin/env python3

"""
fig1B_S1B_ncn_jx_counts_per_sample.py

Python 3 code for creating sorted strip plots of non-core-normal cancer
junction counts per cancer sample per cancer type.

"""

import argparse
import sys
import csv; csv.field_size_limit(sys.maxsize)
from datetime import datetime
import glob
from matplotlib import use; use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
import seaborn as sns; sns.set()
try:
    from utilities.utilities import _TCGA_ABBR, _CANCER_COLORS
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _TCGA_ABBR, _CANCER_COLORS, mm2inch


_SORT_ORDER = [
    'Brain_Lower_Grade_Glioma',
    'Thyroid_Carcinoma',
    'Glioblastoma_Multiforme',
    'Head_and_Neck_Squamous_Cell_Carcinoma',
    'Prostate_Adenocarcinoma',
    'Mesothelioma',
    'Esophageal_Carcinoma',
    'Kidney_Renal_Clear_Cell_Carcinoma',
    'Skin_Cutaneous_Melanoma',
    'Kidney_Chromophobe',
    'Kidney_Renal_Papillary_Cell_Carcinoma',
    'Liver_Hepatocellular_Carcinoma',
    'Breast_Invasive_Carcinoma',
    'Lung_Adenocarcinoma',
    'Stomach_Adenocarcinoma',
    'Pancreatic_Adenocarcinoma',
    'Lung_Squamous_Cell_Carcinoma',
    'Pheochromocytoma_and_Paraganglioma',
    'Sarcoma',
    'Testicular_Germ_Cell_Tumors',
    'Cholangiocarcinoma',
    'Acute_Myeloid_Leukemia',
    'Bladder_Urothelial_Carcinoma',
    'Uterine_Carcinosarcoma',
    'Colon_Adenocarcinoma',
    'Ovarian_Serous_Cystadenocarcinoma',
    'Rectum_Adenocarcinoma',
    'Uterine_Corpus_Endometrial_Carcinoma',
    'Thymoma',
    'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma',
    'Adrenocortical_Carcinoma',
    'Uveal_Melanoma',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma'
]


def double_sorted_stripplot(top_df, bottom_df, out_path, grayed, now):
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 7.0, 4.5
    sns.set_context("paper")
    sns.set_style("whitegrid")
    f, (ax1, ax2) = plt.subplots(
        nrows=2, ncols=1, gridspec_kw = {'hspace':0.06}
    )
    cancer_abbrs = top_df.color.unique()
    color_set = []
    midpoints = {}
    x_locs = []
    x_labels = cancer_abbrs.tolist()
    gray_ranges1 = []
    for abbr in cancer_abbrs:
        color_set.append(_CANCER_COLORS[abbr])
        can_xrange = top_df[top_df.color == abbr]['x_val']
        sides = 75
        if abbr in grayed:
            gray_ranges1.append(
                [can_xrange.min() - sides, can_xrange.max() + sides]
            )

    plt.axes(ax1)
    sns.scatterplot(
        x='x_val', y='neojx_count', hue='color', data=top_df,
        legend=False, palette=sns.xkcd_palette(color_set),
        edgecolor='None', s=4, ax=ax1
    )
    if grayed:
        for gray_range in gray_ranges1:
            plt.axvspan(gray_range[0], gray_range[1], color='gray', alpha=0.1)
    plt.xlabel('')
    plt.xticks([], [])
    plt.setp(ax2.get_xticklabels(), rotation=90, fontsize=6)
    plt.yscale('log')
    y_toplabel = (
        'junctions/sample (#), excluding:\n'
        'tissue-matched GTEx and TCGA normal'
    )
    plt.ylabel(y_toplabel, fontsize=7)

    plt.axes(ax2)
    gray_ranges2 = []
    for abbr in cancer_abbrs:
        can_xrange = bottom_df[bottom_df.color == abbr]['x_val']
        midpoints[abbr] = can_xrange.mean()
        x_locs.append(midpoints[abbr])
        sides = 75
        if abbr in grayed:
            gray_ranges2.append(
                [can_xrange.min() - sides, can_xrange.max() + sides]
            )
    sns.scatterplot(
        x='x_val', y='neojx_count', hue='color', data=bottom_df,
        legend=False, palette=sns.xkcd_palette(color_set),
        edgecolor='None', s=4, ax=ax2
    )

    if grayed:
        for gray_range in gray_ranges2:
            plt.axvspan(gray_range[0], gray_range[1], color='gray', alpha=0.1)

    fig_name = 'figS1B_ncnjx_count_double_sorted_stripplot_{}.pdf'.format(now)

    plt.xlabel('')
    plt.xticks(x_locs, x_labels)
    plt.setp(ax2.get_xticklabels(), rotation=90, fontsize=6)
    plt.yscale('log')

    ax1.xaxis.grid(False)
    ax2.xaxis.grid(False)

    y_botlabel = (
        'junctions/sample (#), excludng:\nall GTEx and TCGA normal'
    )
    plt.ylabel(y_botlabel, fontsize=7)
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    fig.savefig(fig_file)
    plt.close()
    return


def sorted_stripplot(prepped_df, out_path, gray, now):
    """Plots waterfall plot/sorted stripplot.

    Input:
    prepped_df: df prepared by function prep_sorted_stripplot (pandas df)
    out_path: path for storing output files (string)

    Returns none
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = mm2inch((178, 69.8))
    sns.set_context("paper")
    sns.set_style("whitegrid")

    cancer_abbrs = prepped_df.color.unique()
    color_set = []
    midpoints = {}
    x_locs = []
    x_labels = cancer_abbrs.tolist()
    gray_ranges = []
    for abbr in cancer_abbrs:
        color_set.append(_CANCER_COLORS[abbr])
        can_xrange = prepped_df[prepped_df.color == abbr]['x_val']
        midpoints[abbr] = can_xrange.mean()
        x_locs.append(midpoints[abbr])
        sides = 0
        if abbr in gray:
            gray_ranges.append(
                [can_xrange.min() - sides, can_xrange.max() + sides]
            )

    ax = sns.scatterplot(
        x='x_val', y='neojx_count', hue='color', data=prepped_df,
        legend=False, palette=sns.xkcd_palette(color_set),
        edgecolor='None', s=4
    )
    if gray:
        for gray_range in gray_ranges:
            plt.axvspan(gray_range[0], gray_range[1], color='gray', alpha=0.1)

    fig_name = 'fig1B_ncnjx_count_sorted_stripplot_{}.pdf'.format(now)

    plt.xlabel('')
    plt.xticks(x_locs, x_labels)
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=6)
    plt.yscale('log')
    ax.xaxis.grid(False)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
    )
    ylabel_text = ('junctions/sample (#)')
    plt.ylabel(ylabel_text, fontsize=7)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=5, color='black')
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    fig.savefig(fig_file, dpi=300)
    plt.close()
    return


def prep_sorted_stripplot_df(jx_dir, pre_sort):
    """Prepares dataframe for creating waterfall/sorted stripplot.

    Input:
    jx_dir: directory containing jx_indexer query results, .csv files
        containing  sample ids (recount or tcga ids) with neojx top_x for each.
        The file name format is [cancer type]*neojx_counts_per_sample*.csv
        (string)
    k_sort: whether to order cancer types to match the Kahles paper (boolean)
    our_sort: whether to order cancer types by non-paired-normal (useful for
        non-GTEx junction directory) only (boolean)

    If no sort order is passed, first sorts entries by mean neojunction top_x.
    Assigns x value for each sample (per cancer type) by increasing order of
    neojx top_x for each sample.  Adds 150 blank spaces in between cancer
    cancer types.

    Returns prepared df that creates sorted stripplot.
    """
    full_df = pd.DataFrame()
    first_count = 1

    if pre_sort:
        sort_order = _SORT_ORDER
    else:
        means = {}
        for cancer in _TCGA_ABBR:
            jx_glob = os.path.join(
                jx_dir, '{}*neojx_counts_per_sample*'.format(cancer)
            )
            try:
                jx_file = glob.glob(jx_glob)[0]
            except IndexError:
                continue
            curr_df = pd.read_table(jx_file, sep=',')
            mean = curr_df['neojx_count'].mean()
            means[mean] = cancer

        sort_order = []
        for key in sorted(means.keys()):
            sort_order.append(means[key])

    # for i, cancer in enumerate(_SORT_ORDER_KAHLES, 1):
    for i, cancer in enumerate(sort_order, 1):
        abbr = _TCGA_ABBR[cancer]
        inter_space = 150 * i
        jx_glob = os.path.join(
            jx_dir, '{}*neojx_counts_per_sample*'.format(cancer)
        )
        try:
            jx_file = glob.glob(jx_glob)[0]
        except IndexError:
            continue

        curr_df = pd.read_table(jx_file, sep=',').fillna(0)
        curr_df['color'] = abbr
        curr_df = curr_df.sort_values(by=['neojx_count'], ascending=True)
        curr_df = curr_df.reset_index(drop=True)
        try:
            curr_df = curr_df.drop(['recount_id'], axis=1)
        except KeyError:
            pass
        try:
            curr_df = curr_df.drop(['tcga_id'], axis=1)
        except KeyError:
            pass
        curr_df['x_val'] = curr_df.index + first_count
        full_df = pd.concat([full_df, curr_df])
        first_count = len(full_df) + inter_space
    return full_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='generated sorted strip plots of non-core normal '
                    'junctions per patient per cancer type.'
    )
    parser.add_argument(
        '--non-gtex-junction-directory', '-j',
        help='Directory containing .csv files generated by jx_indexer.py '
             'query: should contain neojx (not in gtex) counts for one TCGA '
             'cancer type per file.'
    )
    parser.add_argument(
        '--non-paired-junction-directory', '-p',
        help='Directory containing .csv files generated by jx_indexer.py '
             'query: should contain neojx (not in tissue) counts for one TCGA '
             'cancer type per file.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give path for output figures'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--prepared-sort-order', action='store_true',
        help='Choose this option to plot the cancer types in the same order '
             'as our non-paired-normal jxs.'
    )
    parser.add_argument(
        '--gray', '-g', nargs='*', default=[],
        choices=[
            'LAML', 'ACC', 'BLCA', 'LGG', 'BRCA', 'CESC', 'CHOL', 'COAD',
            'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD',
            'LUSC', 'DLBC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
            'SARC', 'SKCM', 'STAD', 'TGCT', 'THYM', 'THCA', 'UCS', 'UCEC',
            'UVM'
        ],
        help='Choose this option to add gray bar plotting for specified '
             'cancer types.'
    )
    parser.add_argument(
        '--double-plot', '-d', action='store_true',
        help='Choose this option to plot both non-gtex and non-paired normal '
             'juntion counts per sample. If this option is selected, '
             'directories for both types of top_x files must be provided.'
    )

    args = parser.parse_args()
    non_gtex_dir = args.non_gtex_junction_directory
    non_paired_dir = args.non_paired_junction_directory
    out_path = args.output_path
    log_mode = args.log_level
    our_sort = args.prepared_sort_order
    grayed = args.gray
    double_plot = args.double_plot

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    if double_plot:
        if not (non_gtex_dir and non_paired_dir):
            print(
                'both non-paired and non-GTEx directories required; exiting. '
                'Please run again with both arguments included.'
            )
            exit()
        non_paired_df = prep_sorted_stripplot_df(non_paired_dir, our_sort)
        non_gtex_df = prep_sorted_stripplot_df(non_gtex_dir, our_sort)
        double_sorted_stripplot(
            non_paired_df, non_gtex_df, out_path, grayed, now
        )

    else:
        if (non_gtex_dir and non_paired_dir):
            print(
                'only one of non-paired and non-GTEx directories required; '
                'non-GTEx (non-core normal) data will be used to generate '
                'this plot.'
            )
            exit()
        if non_paired_dir:
            jx_dir = non_paired_dir
        else:
            jx_dir = non_gtex_dir

        plot_df = prep_sorted_stripplot_df(jx_dir, our_sort)
        sorted_stripplot(plot_df, out_path, grayed, now)
