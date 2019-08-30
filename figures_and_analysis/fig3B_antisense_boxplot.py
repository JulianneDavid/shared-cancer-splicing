#!/usr/bin/env python3

"""
Fig3B_antisense_boxplots_and_analyze_setmemberships.py

Python 3 code for parsing the results of the set membership assignment files
and creating boxplots for antisense proportions across sets.

"""

import argparse
from datetime import datetime
import logging
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import sys
try:
    from utilities.utilities import get_jx_prev_filename, jx_df_from_file
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import get_jx_prev_filename, jx_df_from_file
from utilities.utilities import _TCGA_ABBR, grouped_boxplots_with_table


def antisense_boxplot(jx_dir, out_path, now):
    """

    :param jx_dir:
    :param all_files:
    :return:
    """
    all_cancers = list(_TCGA_ABBR.keys())
    main_header = 'pan cancer\njunction counts per group (#)'
    data_dict = {main_header: {'data': [[], [], [], [], []]}}
    gtex_total = 0
    adult_total = 0
    sc_total = 0
    dev_total = 0
    un_total = 0
    counts = {'gtex': [], 'adult': [], 'sc': [], 'dev': [], 'un': []}

    can_count = 0
    for cancer in all_cancers:
        logging.info('starting cancer {}'.format(cancer))
        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            continue

        can_count += 1
        jx_df = jx_df_from_file(
            jx_file, 0, 1, True, glob_form=prev_glob, sample=False,
            top_x=False, drop_ann=False
        )
        gtex = jx_df[(jx_df.gtex == 1)]
        adult = jx_df[(jx_df.gtex == 0) & (jx_df.sra_adult == 1)]
        un = jx_df[
            (jx_df.gtex == 0) &
            (jx_df.sra_stemcells == 0) &
            (jx_df.sra_adult == 0) &
            (jx_df.sra_developmental == 0)
        ]

        # Prep boxplots
        dev = jx_df[
            (~jx_df.jx.isin(un.jx)) & (~jx_df.jx.isin(adult.jx)) &
            (~jx_df.jx.isin(gtex.jx)) & (jx_df.sra_developmental == 1)
        ]
        sc = jx_df[
            (~jx_df.jx.isin(un.jx)) & (~jx_df.jx.isin(adult.jx)) &
            (~jx_df.jx.isin(gtex.jx)) &(jx_df.sra_stemcells == 1)
        ]

        gtex_antisense = gtex[gtex.antisense == 1]
        gtex_antiratio = len(gtex_antisense) / len(gtex)
        adult_antisense = adult[adult.antisense == 1]
        adult_anti_ratio = len(adult_antisense) / len(adult)
        dev_antisense = dev[dev.antisense == 1]
        dev_anti_ratio = len(dev_antisense) / len(dev)
        sc_antisense = sc[sc.antisense == 1]
        sc_ratio = len(sc_antisense) / len(sc)
        un_antisense = un[un.antisense == 1]
        un_ratio = len(un_antisense) / len(un)

        gtex_total += len(gtex)
        adult_total += len(adult)
        dev_total += len(dev)
        sc_total += len(sc)
        un_total += len(un)

        counts['gtex'].append(len(gtex))
        counts['adult'].append(len(adult))
        counts['dev'].append(len(dev))
        counts['sc'].append(len(sc))
        counts['un'].append(len(un))

        data_dict[main_header]['data'][0].extend([gtex_antiratio])
        data_dict[main_header]['data'][1].extend([adult_anti_ratio])
        data_dict[main_header]['data'][2].extend([dev_anti_ratio])
        data_dict[main_header]['data'][3].extend([sc_ratio])
        data_dict[main_header]['data'][4].extend([un_ratio])

    print_df = pd.DataFrame({
        'Cancer type': all_cancers,
        'Core normals jx top_x': counts['gtex'],
        'Core normals anti %': data_dict[main_header]['data'][0],
        'Other adult non-cancer jx top_x': counts['adult'],
        'Other adult non-cancer anti %': data_dict[main_header]['data'][1],
        'Developmental jx top_x': counts['dev'],
        'Developmental anti %': data_dict[main_header]['data'][2],
        'Stem cell jc': counts['sc'],
        'Stem cell ar': data_dict[main_header]['data'][3],
        'Unexplained jc': counts['un'],
        'Unexplained ar': data_dict[main_header]['data'][4]
    })
    supp_table = os.path.join(out_path, 'antisense_table_S5.csv')
    with open(supp_table, 'w') as output:
        print_df.to_csv(output, index=False)

    # Assemble boxplot
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 7.0, 2.9
    sns.set_context("paper")
    sns.set_style("whitegrid")

    counts_df = pd.DataFrame(counts)

    table_data = []
    for column in counts_df.columns.values.tolist():
        table_data.append(
            'median: {:,}; IQR: {:,}-{:,}'.format(
                round(counts_df[column].median()),
                int(counts_df[column].quantile([0.25])[0.25]),
                int(counts_df[column].quantile([0.75])[0.75])
            )
        )

    data_dict[main_header]['table_data'] = table_data

    plot_info_dict = {}
    plot_info_dict['light colors'] = [
        'xkcd:kermit green',
        'xkcd:light teal', 'xkcd:pale purple', 'xkcd:apricot',
        'xkcd:light red'
    ]

    plot_info_dict['dark colors'] = plot_info_dict['light colors']
    plot_info_dict['legend'] = [
        'Core normals',
        'Other adult non-cancer',
        'Developmental',
        'Stem cell',
        'Unexplained'
    ]

    plot_info_dict['row colors'] = plot_info_dict['light colors']
    plot_info_dict['row font color'] = [
        'black', 'black', 'black', 'black', 'black'
    ]
    plot_info_dict['row labels'] = plot_info_dict['legend']
    fig_name = 'fig3B_antisense_ratios_boxplot_{}.pdf'.format(now)
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at {}'.format(fig_file))
    grouped_boxplots_with_table(
        data_dict, plot_info_dict, fig_file, logscale=False,
        y_label='antisense junctions (%)', percent=False,
        right_lim_shift=3
    )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Analyzes set memberships and generates antisense figure.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store log file and antisense figure output.'
    )
    parser.add_argument(
        '--full-set-membership-directory', '-s',
        help='Provide the directory containing set membership files generated '
             'by the set membership annotation script.'
    )

    args = parser.parse_args()
    log_mode = args.log_level
    jx_dir = args.full_set_membership_directory
    out_path = args.output_path

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'fig3B_antisense_boxplot_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('command line: {}'.format(' '.join(sys.argv)))

    antisense_boxplot(jx_dir, out_path, now)
