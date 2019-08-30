#!/usr/bin/env python3

"""
figS5_S9A_SRA_cancer_shared_prevalence.py

Python 3 code for comparing TCGA cancer-type cohort prevalences for junctions
not found in core normals and either (a) found in cancer type-matched SRA
samples or (b) not found in these samples.

all junctions:
time python ../jxtapose_experiments/figS5_S9A_sra_cancer_sharedness_boxplots.py
--ontology-df sra_cancer_comp/r -o paper1/results/
-d paper1/not_gencode_filtered_data/prevalence_files/
-e /sra_cancer_comp/exptlists/


unexplained junctions:
paper1 davidju$ time python
../../jxtapose_experiments/figS5_S9A_sra_cancer_sharedness_boxplots.py
--ontology-df ../sra_cancer_comp/raw_results/ -o results/
-d not_gencode_filtered_data/set_memberships_final/unexplained/
-e ../sra_cancer_comp/exptlists/ --jx-type 'unexplained'

"""

import argparse
from datetime import datetime
import glob
import logging
import os
from scipy import stats
import sys
try:
    from utilities.utilities import _ABBR_TO_CAN, _PER, snaptron_results_to_jxs
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _ABBR_TO_CAN, _PER, snaptron_results_to_jxs
from utilities.utilities import get_jx_prev_filename, collect_metasra_count
from utilities.utilities import  jx_df_from_file, grouped_boxplots_with_table


_TCGA_TO_SRA = {
    'THCA': 'sra_thca',
    'THYM': 'sra_thym',
    'UVM':'sra_uvm',
    'LIHC': 'sra_lihc',
    'CHOL': 'sra_chol',
    'GBM': 'sra_gbm',
    'LAML': 'sra_laml',
    'ACC': 'sra_acc',
    'BRCA': 'sra_brca_ductal',
    'CESC': 'sra_cesc',
    'COAD': 'sra_coad',
    'HNSC': 'sra_hnsc',
    'LUAD': 'sra_luad',
    'LUSC': 'sra_lusc',
    'DLBC': 'sra_dlbc',
    'PRAD': 'sra_prad',
    'KIRC': 'sra_renalcarc',
    'SARC': 'sra_sarc',
    'SKCM': 'sra_skcm',
    'TGCT': 'sra_tgct'
}

_SRA_TO_TCGA = {value: key for key, value in _TCGA_TO_SRA.items()}

_DEC_ORDER_UNEXPL = [
    'KIRC', 'CHOL', 'UVM', 'THYM', 'GBM', 'THCA', 'DLBC', 'PRAD', 'LAML',
    'LIHC', 'SARC', 'LUAD', 'TGCT', 'CESC',
]
_REJECT_UNEXPLAINED = ['SARC', 'TGCT', 'CESC', 'LUAD']
_REJECT_ALL = ['LUAD', 'TGCT', 'CESC']
_DEC_ORDER_ALL = [
    'THYM', 'CHOL','UVM', 'GBM', 'KIRC', 'THCA', 'DLBC', 'PRAD', 'LAML',
    'LIHC', 'SARC', 'LUAD', 'TGCT', 'CESC',
]

_SRA_CAN_ABBR = {
    'Acute_Myeloid_Leukemia': 'sra_laml',
    'Adrenocortical_Carcinoma': 'sra_acc',
    'Breast_Invasive_Carcinoma_ductal': 'sra_brca_ductal',
    'Cervical_Carcinoma': 'sra_cesc',
    'Cholangiocarcinoma': 'sra_chol',
    'Colon_Adenocarcinoma': 'sra_coad',
    'Glioblastoma_Multiforme': 'sra_gbm',
    'Head_and_Neck_Squamous_Cell_Carcinoma': 'sra_hnsc',
    'Liver_Hepatocellular_Carcinoma': 'sra_lihc',
    'Lung_Adenocarcinoma': 'sra_luad',
    'Lung_Squamous_Cell_Carcinoma': 'sra_lusc',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma': 'sra_dlbc',
    'no_name': 'no_name',
    'Prostate_Adenocarcinoma': 'sra_prad',
    'Renal_Cell_Carcinoma': 'sra_renalcarc',
    'Sarcoma': 'sra_sarc',
    'Skin_Cutaneous_Melanoma': 'sra_skcm',
    'Testicular_Germ_Cell_Tumors': 'sra_tgct',
    'Thymoma': 'sra_thym',
    'Thyroid_Carcinoma': 'sra_thca',
    'Uveal_Melanoma': 'sra_uvm'
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare TCGA cancer-type cohort prevalences for '
                    'junctions not found in core normals and either (a) found '
                    'in cancer type-matched SRA samples or (b) not found in '
                    'these samples.'
    )
    parser.add_argument(
        '--snaptron-results',
        help='directory containing multiple .txt files containing results of '
             'a previous snaptron search for SRA cancer samples.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output figure and log file, including statistics.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--database-junction-directory', '-d',
        help='FOR FIGURE S5: directory containing .csv files with junctions '
             'extracted via a jx_indexer query, each containing prevalence '
             'values for one cancer type, to analyze all non-core normal '
             'junctions. FOR FIGURE S9A: this should be the "unexplained" '
             'directory output from running set_membership_annotation.py.'
    )
    parser.add_argument(
        '--exptlist-directory', '-e',
        help='If snaptron junctions have already been collected, provide the '
             'directory where the lists of recount/snaptron-available SRA '
             'experiments for each cell type are stored.'
    )
    parser.add_argument(
        '--unexplained-junctions', action='store_true',
        help='Select this option for only unexplained junctions also not '
             'found in any selected SRA non-cancer samples) for figure S9A.'
    )

    args = parser.parse_args()
    snap_dir = args.snaptron_results
    out_path = args.output_path
    log_mode = args.log_level
    jx_dir = args.database_junction_directory
    exptlist_dir = args.exptlist_directory
    unexplained_jxs = args.unexplained_junctions

    if unexplained_jxs:
        dec_order = _DEC_ORDER_UNEXPL
        fig_flag = 'figS9A_unexplained_jxs'
        fig_size = 4.68, 4.0
    else:
        dec_order = _DEC_ORDER_ALL
        fig_flag = 'figS5_all_jxs'
        fig_size = 7.0, 4.0

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, '{}_SRA-cancer-shared-prevs_log_{}.txt'.format(fig_flag, now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    sra_abbr_to_nametag = {}
    for nametag, abbr in _SRA_CAN_ABBR.items():
        sra_abbr_to_nametag[abbr] = nametag

    grouped_data_dict = {}
    pvals = []
    H_stats = []
    for cancer_abbr in dec_order:
        cancer = _ABBR_TO_CAN[cancer_abbr]
        sra_abbr = _TCGA_TO_SRA[cancer_abbr]
        logging.info('starting {}'.format(cancer))

        can_file, flag, all_jxs_name = get_jx_prev_filename(jx_dir, cancer)

        if flag == 'all':
            if cancer_abbr in _REJECT_ALL:
                continue
        elif flag == 'unexpl':
            if cancer_abbr in _REJECT_UNEXPLAINED:
                continue

        name_tag = sra_abbr_to_nametag[sra_abbr]
        snap_glob = os.path.join(
            snap_dir, '{}_rawresults*.txt'.format(name_tag)
        )
        try:
            snap_file = glob.glob(snap_glob)[0]
        except IndexError:
            logging.info('snaptron results globbed file not found:')
            logging.info(snap_glob)

        sra_expt_count = collect_metasra_count(name_tag, exptlist_dir)
        with open(snap_file) as lines:
            sra_can_jxs = set(
                snaptron_results_to_jxs(lines, min_sample_count=1)
            )

        jx_df = jx_df_from_file(
            can_file, 0.01, 1.0, chunk_it=True, glob_form=all_jxs_name,
            sample=False, top_x=False, drop_ann=True, cancer=cancer
        )
        per_col = cancer + _PER

        jx_df['sra_can_type'] = jx_df.jx.apply(
            lambda x: x in sra_can_jxs
        ).astype(int)

        in_sra = jx_df[jx_df['sra_can_type'] == 1][per_col].tolist()
        non_sra = jx_df[jx_df['sra_can_type'] == 0][per_col].tolist()
        stat, pval = stats.kruskal(in_sra, non_sra)
        logging.info("kruskal statistic:".format(stat))
        logging.info("kruskal t-test p-value:".format(pval))
        pvals.append(pval)
        H_stats.append(stat)
        if pval < 0.00001:
            pval_for_table = '<0.00001'
        else:
            pval_for_table = '{:.4g}'.format(pval)

        grouped_data_dict[cancer_abbr] = {}
        grouped_data_dict[cancer_abbr]['data'] = [in_sra, non_sra]
        grouped_data_dict[cancer_abbr]['table_data'] = [
            sra_expt_count, len(in_sra), len(non_sra), pval_for_table
        ]

    logging.info('ranges for all cancer types calculated:')
    logging.info('H-stat range: {}-{}'.format(min(H_stats), max(H_stats)))
    logging.info('p-value range: {}-{}'.format(min(pvals), max(pvals)))

    plot_info_dict = {}
    plot_info_dict['light colors'] = ['xkcd:tangerine', 'xkcd:cerulean']
    plot_info_dict['dark colors'] = ['xkcd:pumpkin', 'xkcd:ocean blue']
    plot_info_dict['row colors'] = [
        'white', 'xkcd:pumpkin', 'xkcd:ocean blue', 'white'
    ]
    plot_info_dict['row font color'] = ['black', 'white', 'white', 'black']
    plot_info_dict['row labels'] = [
        'SRA expts', 'SRA-shared\nneojunctions', 'TCGA-only\nneojunctions',
        'p value'
    ]

    fig_name = '{}_SRA-cancer-shared_boxplots_{}.pdf'.format(fig_flag, now)
    fig_file = os.path.join(out_path, fig_name)
    logging.info('saving figure at {}'.format(fig_file))

    grouped_boxplots_with_table(
        grouped_data_dict, plot_info_dict, fig_file, fig_size
    )
