#!/usr/bin/env python3

"""
figS7_SKCM_jx_similarity.py

Python 3 code for comparing sample-wise junction sets across groups: TCGA
cancer (SKCM), tissue-matched samples from GTEx and TCGA, and SRA tissue of
origin (melanocytes).

"""

import argparse
from datetime import datetime
import glob
import json
import logging
from matplotlib import use; use('pdf')
import os
import pandas as pd
import sys
try:
    from utilities.utilities import grouped_boxplots_with_table, _SRA_ABBR
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import grouped_boxplots_with_table, _SRA_ABBR
from utilities.utilities import _MATCHED_NORMALS, _ABBR_TO_CAN


def snaptron_results_to_sets(snaptron_result_lines):
    """

    :param snaptron_result_lines:
    :return:
    """
    jx_sets = {}
    for line in snaptron_result_lines:
        items = line.split('\t')
        if not items[0].endswith('I'):
            continue

        chrom, left, right, strand = items[2], items[3], items[4], items[6]
        left = str(int(left) - 1)
        right = str(int(right) - 1)
        jx = ';'.join([chrom, left, right, strand])

        samples = items[12].split(',')
        samples.remove('')

        for sample in samples:
            samp_id, cov = sample.split(':')
            try:
                jx_sets[samp_id].append(jx)
            except KeyError:
                jx_sets[samp_id] = [jx]

    for sample, junctions in jx_sets.items():
        jx_sets[sample] = set(junctions)

    return jx_sets


def collect_db_sample_jx_sets(name_tag, tumor_normal, db_directory):
    """

    :param name_tag:
    :param tumor_normal:
    :param db_directory:
    :return:
    """
    max_chunk = 1e6
    glob_form = '{}_all_jxs_{}_*.txt'.format(name_tag, tumor_normal)
    file_path = os.path.join(db_directory, glob_form)
    jx_file = glob.glob(file_path)[0]
    jx_df = pd.DataFrame()
    chunks = pd.read_table(jx_file, sep=',', chunksize=max_chunk)
    for chunk in chunks:
        chunk = chunk.fillna(0)
        jx_df = pd.concat([jx_df, chunk])

    samples = jx_df.groupby('recount_id')
    set_dict = dict(samples.apply(lambda x: set(x['jx'])))

    return set_dict


def compare_jx_sets(jx_set_1, jx_set_2):
    overlap_vals = []
    for s1 in jx_set_1.values():
        for s2 in jx_set_2.values():
            min_len = min(len(s1), len(s2))
            overlap_vals.append(len(s1.intersection(s2)) / min_len)

    return overlap_vals


def collect_sample_overlaps(snap_dir, out_path, abbr):
    """
    
    :param snap_dir: 
    :param out_path: 
    :param abbr: 
    :return: 
    """
    json_file = os.path.join(out_path, '{}_figS7_json.txt'.format(abbr))
    if os.path.exists(json_file):
        logging.info('json file exists - starting recovery.')
        with open(json_file) as recover:
            overlaps, sample_counts = json.load(recover)
            for set in overlaps:
                try:
                    while 0.0 in set:
                        set.remove(0.0)
                except ValueError:
                    continue
        logging.info('json file data recovery complete.')

    else:
        logging.info('loading SRA junctions per sample')
        if not snap_dir.endswith('.txt'):
            txt_path = os.path.join(snap_dir, '*rawresults*.txt')
            snap_files = glob.glob(txt_path)
        else:
            snap_files = [snap_dir]

        sra_norms = {}
        for jx_file in snap_files:
            name_tag = os.path.basename(jx_file).split('.')[0]
            name_tag = name_tag.split('_rawresults')[0]
            try:
                name_tag = name_tag.split('metaSRA-runs_')[1]
            except IndexError:
                pass
            if _SRA_ABBR[name_tag] not in sra_abbrs:
                continue
            with open(jx_file) as lines:
                sra_norms.update(
                    snaptron_results_to_sets(lines)
                )
            logging.info(sra_norms.keys())

        logging.info('finished loading SRA sample junctions.')
        logging.info('loading GTEx sample junctions')
        gtex_norms = collect_db_sample_jx_sets(gtex_match, 'normal', db_dir)
        logging.info('finished loading GTEx sample junctions')
        logging.info('starting loading TCGA cancer sample junctions')
        tcga_cans = collect_db_sample_jx_sets(cancer, 'tumor', db_dir)
        logging.info('finished loading TCGA cancer sample junctions')
        logging.info('loading TCGA normal sample junctions')
        tcga_norms = collect_db_sample_jx_sets(cancer, 'normal', db_dir)
        logging.info('finished loading TCGA normal sample junctions')

        logging.info('starting cancer-sra overlaps:')
        can_sra_overlaps = compare_jx_sets(tcga_cans, sra_norms)
        logging.info('starting cancer-gtex overlaps:')
        can_gtex_overlaps = compare_jx_sets(tcga_cans, gtex_norms)
        logging.info('starting cancer-matched normal overlaps:')
        can_tnorm_overlaps = compare_jx_sets(tcga_cans, tcga_norms)
        logging.info('starting cancer-cancer overlaps')
        can_can_overlaps = compare_jx_sets(tcga_cans, tcga_cans)

        logging.info('all overlaps calculated!')
        overlaps = [
            can_can_overlaps, can_tnorm_overlaps, can_sra_overlaps,
            can_gtex_overlaps
        ]
        sample_counts = [
            len(tcga_cans), len(tcga_norms), len(sra_norms), len(gtex_norms)
        ]
        logging.info('starting json dump.')
        with open(json_file, 'w') as dump_file:
            json.dump([overlaps, sample_counts], dump_file)
        logging.info('json dump finished.')
        
    return overlaps, sample_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare samplewise all junctions per sample across '
                    'groups: TCGA cancer samples with TCGA SKCM matched '
                    'normal samples, SRA melanocyte samples, and GTEx bulk '
                    'skin samples.'
    )
    parser.add_argument(
        '--snaptron-results',
        help='.txt file containing results from a previous snaptron search, '
             'or directory containing multiple of these.'
    )
    parser.add_argument(
        '--db-jx-dir', '-d',
        help='Give the directory containing files with all junctions for each '
             'sample for appropriate cancer.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output log file and figure.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )

    args = parser.parse_args()
    snap_dir = args.snaptron_results
    db_dir = args.db_jx_dir
    out_path = args.output_path
    log_mode = args.log_level

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'SKCM_jx_similarities_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    grouped_data_dict = {}
    abbr = 'SKCM'
    cancer = _ABBR_TO_CAN[abbr]
    gtex_match = _MATCHED_NORMALS[cancer][0]
    sra_abbrs = {'mel_cl', 'mel_pc'}
    sra_name = 'melanocyte'
    gtex_name = gtex_match.lower()

    overlaps, sample_counts = collect_sample_overlaps(snap_dir, out_path, abbr)

    logging.info('starting data dict assembly.')
    col_header = (
        '{}\nGTEx match: {}\nSRA match: {}'
        ''.format(abbr, gtex_name, sra_name)
    )
    grouped_data_dict[col_header] = {
        'data': overlaps, 'table_data': sample_counts
    }

    logging.info('creating plot info dict')
    plot_info_dict = {}
    plot_info_dict['light colors'] = [
        'xkcd:dark tan', 'xkcd:apple green', 'xkcd:rose pink',
        'xkcd:dusty blue'
    ]
    plot_info_dict['dark colors'] = [
        'xkcd:earth', 'xkcd:grass', 'xkcd:deep pink', 'xkcd:denim'
    ]
    plot_info_dict['legend'] = [
        '# tumor samples', '# normal samples', '# SRA matched samples',
        '# GTEx matched samples'
    ]
    plot_info_dict['row colors'] = plot_info_dict['light colors']
    plot_info_dict['row font color'] = ['black', 'black', 'black', 'black']
    plot_info_dict['row labels'] = plot_info_dict['legend']
    fig_name = 'figS7_SKCM_jx_overlap_boxplot_{}.pdf'.format(now)
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at {}'.format(fig_file))
    grouped_boxplots_with_table(
        grouped_data_dict, plot_info_dict, fig_file, logscale=False,
        y_label='proportion of junctions shared\nbetween individual samples',
        right_lim_shift=2.6
    )

