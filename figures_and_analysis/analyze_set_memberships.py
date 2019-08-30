#!/usr/bin/env python3

"""
analyze_set_memberships.py

Python 3 code for parsing the results of the set membership assignment files
and creating prevalence files for TCGA cancer types that do not include
junctions with only a single TCGA read.

"""

import argparse
from datetime import datetime
import logging
import os
import sys
try:
    from utilities.utilities import get_jx_prev_filename, jx_df_from_file
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import get_jx_prev_filename, jx_df_from_file
from utilities.utilities import _TCGA_ABBR, _PER


def analyze_sets(jx_dir, out_path, now):
    """Parses set membership files and returns various overall statistics

    Input:
    jx_dir: (string) path to directory containing set membership

    """
    all_files = list(_TCGA_ABBR.keys())
    shared_threshold = 0.05
    can_count = 0
    all_set = set()
    gtex_set = set()
    non_gtex_set = set()
    single_count_dict = {}
    percent_unexplained_sum = []
    percent_developmental_sum = []
    percent_stemcell_sum = []
    percent_otheradult_sum = []
    percent_antisense_sum = []
    percent_non_tis_in_gtex_sum = []
    percent_shared_sum = []
    percent_non_gtex_sum = []
    percent_in_gtex_sum = []
    percent_cancer_loci_sum = []
    all_shared_set = set()
    distinguishing_shared = set()
    over_threshold = set()
    multi_cancer_over_thresh = set()
    single_sample_set = set()
    multi_cancer_any = set()
    for cancer in all_files:
        logging.info('starting cancer {}'.format(cancer))
        jx_file, flag, prev_glob = get_jx_prev_filename(jx_dir, cancer)
        if not jx_file:
            continue

        can_count += 1
        jx_df = jx_df_from_file(
            jx_file, 0, 1, True, glob_form=prev_glob, sample=False,
            top_x=False, drop_ann=False
        )
        per_col = cancer + _PER
        min_sharedness = jx_df[per_col].min()

        all_set.update(jx_df.jx)
        all_gtex = jx_df[jx_df.gtex == 1]
        gtex_set.update(all_gtex.jx)
        neojxs = jx_df[jx_df.paired == 0]
        nongtex = jx_df[jx_df.gtex == 0]

        nongtex_to_write = nongtex[['jx', 'annotation', per_col]]
        new_nongtex_file = os.path.join(
            out_path, '{}_not-in-GTEx_jxs_{}.csv'.format(cancer, now)
        )
        with open(new_nongtex_file, 'w') as output:
            nongtex_to_write.to_csv(output, index=False)

        curr_non_gtex_jxs = nongtex['jx']
        multi_cancer_any.update(non_gtex_set.intersection(curr_non_gtex_jxs))
        non_gtex_set.update(curr_non_gtex_jxs)
        unexpl = nongtex[
            (nongtex.sra_adult == 0) &
            (nongtex.sra_developmental == 0) &
            (nongtex.sra_stemcells == 0)
        ]
        non_gtex_count = len(nongtex)

        percent_unexplained_sum.append(len(unexpl) / non_gtex_count)

        non_gtex_count = len(nongtex)
        dev_df = nongtex[nongtex.sra_developmental == 1]
        percent_developmental_sum.append(len(dev_df) / non_gtex_count)
        sc_df = nongtex[nongtex.sra_stemcells == 1]
        percent_stemcell_sum.append(len(sc_df) / non_gtex_count)
        oth_adult_df = nongtex[nongtex.sra_adult == 1]
        percent_otheradult_sum.append(len(oth_adult_df) / non_gtex_count)
        anti_df = unexpl[unexpl.antisense == 1]
        percent_antisense_sum.append(len(anti_df) / len(unexpl))
        num_gtex_nonmatch = len(neojxs[neojxs.gtex == 1])
        percent_non_tis_in_gtex_sum.append(num_gtex_nonmatch / len(neojxs))
        percent_non_gtex_sum.append(len(nongtex) / len(jx_df))
        percent_in_gtex_sum.append(len(all_gtex) / len(jx_df))
        can_locus_df = unexpl[unexpl.cancer_locus == 1]
        percent_cancer_loci_sum.append(len(can_locus_df) / len(unexpl))

        shared_nongtex = nongtex[nongtex[per_col] > min_sharedness]
        percent_shared_sum.append(len(shared_nongtex) / len(nongtex))
        all_shared_set.update(shared_nongtex['jx'])
        distinguishing_shared.symmetric_difference_update(shared_nongtex['jx'])

        high_shared = nongtex[nongtex[per_col] > shared_threshold]
        curr_jxs = set(high_shared['jx'])
        multi_cancer_over_thresh.update(over_threshold.intersection(curr_jxs))
        over_threshold.update(curr_jxs)

        single_sample = nongtex[nongtex[per_col] == min_sharedness]
        for jx in single_sample['jx']:
            try:
                single_count_dict[jx] += 1
            except KeyError:
                single_count_dict[jx] = 1

        single_sample_set.symmetric_difference_update(single_sample['jx'])

    total_single_sample_jxs = 0
    for jx, count in single_count_dict.items():
        if count == 1:
            total_single_sample_jxs += 1

    logging.info('summary:')
    logging.info(
        '{}% of non-tissue matched junctions, on average per cancer type, are '
        'present in other GTEx and TCGA normal samples.'
        ''.format(sum(percent_non_tis_in_gtex_sum) / can_count)
    )
    logging.info(
        '{}% of non-core-normal junctions are shared by multiple samples '
        'within a cancer type; {}% of these are unique to one cancer type'
        ''.format(
            len(all_shared_set) / len(non_gtex_set),
            len(distinguishing_shared) / len(all_shared_set)
        )
    )
    logging.info(
        '{}% of non-core-normal junctions on average per cancer type are in '
        'developmental samples'
        ''.format(sum(percent_developmental_sum) / can_count)
    )
    logging.info('result section 1:')
    logging.info(
        '{}% of non-tissue matched junctions, on average per cancer type, are '
        'present in other GTEx and TCGA normal samples.'
        ''.format(sum(percent_non_tis_in_gtex_sum) / can_count)
    )
    logging.info(
        '{}% of cancer junctions are in core normals'
        ''.format(sum(percent_in_gtex_sum) / can_count)
    )
    logging.info(
        '{}% of non-core normal junctions are from individual samples only'
        '{}% are shared by >{}% of samples of at least one cancer type'
        ''.format(
            total_single_sample_jxs / len(non_gtex_set),
            len(over_threshold) / len(non_gtex_set),
            shared_threshold * 100
        )
    )
    logging.info(
        '{}% of non-core normal junctions are shared between multiple cancer '
        'types; {} junctions are in >{}% of samples of more than 1 TCGA cancer '
        'type'.format(
            len(multi_cancer_any) / len(non_gtex_set),
            len(multi_cancer_over_thresh),
            shared_threshold * 100
        )
    )
    logging.info('results section 3:')
    logging.info('On average per cancer type:')
    logging.info(
        '{}% of non-core-normal junctions are in SRA developmental samples'
        ''.format(sum(percent_developmental_sum) / can_count)
    )
    logging.info(
        '{}% of non-core-normal junctions are in SRA non-cancer adult samples'
        ''.format(sum(percent_otheradult_sum) / can_count)
    )
    logging.info(
        '{}% of non-core-normal junctions are in SRA stem cell samples'
        ''.format(sum(percent_stemcell_sum) / can_count)
    )
    logging.info(
        '{}% of junctions on average per cancer are unexplained (not found in '
        'any of the target tissue or cell types)'
        ''.format(sum(percent_unexplained_sum) / can_count)
    )
    logging.info(
        'Of the unexplained junctions: {} on average per cancer type are '
        'associated with cancer genes'
        ''.format(sum(percent_cancer_loci_sum) / can_count)
    )
    logging.info(
        'Of unexplained junctions: {}% on average per cancer type are in '
        'antisense transcripts'
        ''.format(sum(percent_antisense_sum) / can_count)
    )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Analyzes set memberships.'
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
        out_path, 'set_membership_analysis_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('command line: {}'.format(' '.join(sys.argv)))

    prev_dir = os.path.join(jx_dir, 'true_TCGA_prevalence_files')
    os.makedirs(prev_dir, exist_ok=True)

    analyze_sets(jx_dir, prev_dir, now)
