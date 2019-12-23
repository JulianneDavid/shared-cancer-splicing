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
from numpy import std
import os
from statistics import mean
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
    unexplained_list = []
    developmental_list = []
    in_stemcell_list = []
    in_otheradult_list = []
    antisense_list = []
    non_tis_in_gtex_list = []
    percent_shared_sum = []
    percent_non_gtex_sum = []
    in_gtex_list = []
    in_cancerloci_list = []
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
            top_x=False, drop_ann=False, cancer=cancer
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

        unexplained_list.append(len(unexpl) / non_gtex_count)

        non_gtex_count = len(nongtex)
        dev_df = nongtex[nongtex.sra_developmental == 1]
        developmental_list.append(len(dev_df) / non_gtex_count)
        sc_df = nongtex[nongtex.sra_stemcells == 1]
        in_stemcell_list.append(len(sc_df) / non_gtex_count)
        oth_adult_df = nongtex[nongtex.sra_adult == 1]
        in_otheradult_list.append(len(oth_adult_df) / non_gtex_count)
        anti_df = unexpl[unexpl.antisense == 1]
        antisense_list.append(len(anti_df) / len(unexpl))
        num_gtex_nonmatch = len(neojxs[neojxs.gtex == 1])
        non_tis_in_gtex_list.append(num_gtex_nonmatch / len(neojxs))
        percent_non_gtex_sum.append(len(nongtex) / len(jx_df))
        in_gtex_list.append(len(all_gtex) / len(jx_df))
        can_locus_df = unexpl[unexpl.cancer_locus == 1]
        in_cancerloci_list.append(len(can_locus_df) / len(unexpl))

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

    logging.info('abstract:')
    logging.info(
        'averaging across cancer types, {:.1%} of exon-exon junctions '
        'thought to be cancer-specific based on comparison with '
        'tissue-matched samples (σ = {:.1%}) are in fact present in other '
        'adult non-cancer tissues throughout the body'
        ''.format(mean(non_tis_in_gtex_list), std(non_tis_in_gtex_list))
    )
    logging.info(
        '{:.1%} of junctions not present in any GTEx or TCGA normal tissues '
        'are shared by multiple samples within at least one cancer type '
        'cohort, and {:.1%} of these distinguish between different cancer '
        'types'.format(
            len(all_shared_set) / len(non_gtex_set),
            len(distinguishing_shared) / len(all_shared_set)
        )
    )
    logging.info(
        'many of these junctions not found in GTEx or TCGA normal tissues '
        '({:.1%} on average, σ = {:.1%}) are also found in embryological '
        'and other developmentally associated cells'
        ''.format(mean(developmental_list), std(developmental_list))
    )
    logging.info('result section 1:')
    logging.info(
        'We found that on average, across cancer types, {:.1%} of junctions '
        'potentially thought to be cancer-specific based on comparison only '
        'with tissue-matched samples (σ = {:.1%}) are in fact present in '
        'other adult non-cancer tissues and cell types throughout the body.'
        ''.format(mean(non_tis_in_gtex_list), std(non_tis_in_gtex_list))
    )
    logging.info(
        'Across cancer types, an average of {:.1%} of all junctions found in '
        'cancer samples (σ = {:.1%}) are also present in one or more adult '
        'normal samples from GTEx or TCGA [“core normals”]'
        ''.format(mean(in_gtex_list), std(in_gtex_list))
    )
    logging.info(
        'We observed that over half ({:.1%}) of these junctions are confined '
        'to individual samples, although a small but significant subset '
        '({:.2%}) is shared across at least 5% of samples in at least one '
        'cancer-type cohort'.format(
            total_single_sample_jxs / len(non_gtex_set),
            len(over_threshold) / len(non_gtex_set), shared_threshold * 100
        )
    )
    logging.info(
        'We also noted that {:.1%} of novel junctions are shared between '
        'multiple cancer types, with a total of {} junctions present in at '
        'least {:,} of samples each across two or more TCGA cancer cohorts'
        ''.format(
            len(multi_cancer_any) / len(non_gtex_set),
            len(multi_cancer_over_thresh), shared_threshold * 100
        )
    )
    logging.info('results section 3:')
    logging.info('On average per cancer type:')
    logging.info(
        '{:.1%} of these cancer junctions (σ = {:.1%}) occur in SRA '
        'developmental cell or tissue samples.'
        ''.format(mean(developmental_list), std(developmental_list))
    )
    logging.info(
        'We also considered samples from SRA normal stem cell samples and '
        'from selected SRA normal adult tissues and cell types: {:.1%} '
        '(σ = {:.1%}) and {:.1%} (σ = {:.1%}) of cancer junctions not '
        'found in core normals occur in stem cell and selected adult tissues'
        ''.format(
            mean(in_stemcell_list), std(in_stemcell_list),
            mean(in_otheradult_list), std(in_otheradult_list)
        )
    )
    logging.info(
        'The remaining significant majority of these cancer junctions not '
        'found in core normals were also not present found in any non-cancer '
        'SRA tissue or cell type studied ({:.1%} on average per cancer type '
        'cohort (σ = {:.1%})'
        ''.format(mean(unexplained_list), std(unexplained_list))
    )
    logging.info(
        'Among all otherwise unexplained junctions, an average of {:.2%} '
        '(σ = {:.2%}) across cancer types are associated with known '
        'cancer-predisposing or cancer-relevant loci.'
        ''.format(mean(in_cancerloci_list), std(in_cancerloci_list))
    )
    logging.info(
        'Further, an elevated proportion of otherwise unexplained junctions '
        '(on average, {:.1%}, σ = {:.1%}) occur in likely antisense '
        'transcripts'.format(mean(antisense_list), std(antisense_list))
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
