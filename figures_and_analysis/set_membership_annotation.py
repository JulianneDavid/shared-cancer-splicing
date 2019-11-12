#!/usr/bin/env python3

"""
set membership_annotation.py

Python 3 code for annotating a junction list with presence in annotation,
gene boundaries, antisense transcripts, and presence in various adult, stem
cell, and developmental samples.

"""

import argparse
from datetime import datetime
import glob
import json
import logging
import os
import pandas as pd
import sqlite3 as sql
import sys
try:
    import utilities.utilities as jx_util
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    import utilities.utilities as jx_util
from utilities.utilities import _TCGA_CANCER_TYPES, _SRA_ADULT, _SRA_DEV
from utilities.utilities import _SRA_STEMCELLS
from utilities.utilities import snaptron_results_to_jxs, create_jx_id_map
from utilities.utilities import  make_id_name_dict, gtf_to_cds, cds_to_tree
from utilities.utilities import cds_to_antisense,  extract_splice_sites
from utilities.utilities import jx_df_from_file, jx_gene_overlap


def collect_cancer_loci(cancer_gene_file, symbol_column, separator):
    can_gene_df = pd.read_table(
        cancer_gene_file, sep=separator, usecols=[symbol_column]
    )
    cancer_loci = set(can_gene_df[symbol_column].tolist())
    return cancer_loci


def check_cancer_loci(coding_region_entry, cancer_locus_set):
    cancer_locus_presence = 0
    regions = coding_region_entry.split(';')
    all_genes = []
    for gene_set in regions:
        for gene in gene_set.split(','):
            all_genes.append(gene)
    for gene in all_genes:
        if gene in cancer_locus_set:
            cancer_locus_presence = 1

    return cancer_locus_presence


def load_non_cancer_SRA_jxs(snap_in, min_samps, min_reads, overall_set_count):
    """

    :param snap_in:
    :return:
    """
    sra_adult_jxs = set()
    sra_dev_jxs = set()
    sra_sc_jxs = set()
    sra_emb_all = set()
    sra_fet = set()
    sra_zyg = set()
    sra_oo = set()
    sra_plc = set()
    sra_emb_ect = set()
    sra_emb_emb = set()
    sra_emb_late = set()
    sra_emb_mes = set()
    sra_emb_myo = set()
    dev_sets = [
        sra_dev_jxs, sra_emb_all, sra_fet, sra_zyg, sra_oo, sra_plc,
        sra_emb_ect, sra_emb_emb, sra_emb_late, sra_emb_mes, sra_emb_myo
    ]

    if not snap_in.endswith('.txt'):
        txt_path = os.path.join(snap_in, '*rawresults*.txt')
        ont_files = glob.glob(txt_path)
    else:
        ont_files = [snap_in]

    adult_samp_counts = {}
    sc_samp_counts = {}
    dev_samp_counts = {}
    for ont in ont_files:
        name_tag = os.path.basename(ont).split('.')[0]
        name_tag = name_tag.split('_rawresults')[0]
        try:
            name_tag = name_tag.split('metaSRA-runs_')[1]
        except IndexError:
            pass
        with open(ont) as lines:
            jxs = snaptron_results_to_jxs(
                lines, min_sample_count=min_samps, min_read_count=min_reads
            )
            if name_tag in _SRA_ADULT:
                sra_adult_jxs.update(jxs)
                for jx, samp_count in jxs.items():
                    try:
                        adult_samp_counts[jx] += samp_count
                    except KeyError:
                        adult_samp_counts[jx] = samp_count

            if name_tag in _SRA_DEV:
                sra_dev_jxs.update(jxs)
                for jx, samp_count in jxs.items():
                    try:
                        dev_samp_counts[jx] += samp_count
                    except KeyError:
                        dev_samp_counts[jx] = samp_count

            if name_tag in _SRA_STEMCELLS:
                sra_sc_jxs.update(jxs)
                for jx, samp_count in jxs.items():
                    try:
                        sc_samp_counts[jx] += samp_count
                    except KeyError:
                        sc_samp_counts[jx] = samp_count
            if name_tag in jx_util._SRA_EMB_ALL:
                sra_emb_all.update(jxs)
            if name_tag in jx_util._SRA_FET:
                sra_fet.update(jxs)
            if name_tag in jx_util._SRA_ZYG:
                sra_zyg.update(jxs)
            if name_tag in jx_util._SRA_OO:
                sra_oo.update(jxs)
            if name_tag in jx_util._SRA_PLC:
                sra_plc.update(jxs)
            if name_tag in jx_util._SRA_EMB_ECT:
                sra_emb_ect.update(jxs)
            if name_tag in jx_util._SRA_EMB_EMB:
                sra_emb_emb.update(jxs)
            if name_tag in jx_util._SRA_EMB_LATE:
                sra_emb_late.update(jxs)
            if name_tag in jx_util._SRA_EMB_MES:
                sra_emb_mes.update(jxs)
            if name_tag in jx_util._SRA_EMB_MYO:
                sra_emb_myo.update(jxs)

    uncategorized_jxs = set()
    if overall_set_count > 1:
        for jx, samp_count in adult_samp_counts.items():
            if samp_count < overall_set_count:
                sra_adult_jxs.remove(jx)
                uncategorized_jxs.update([jx])
        for jx, samp_count in sc_samp_counts.items():
            if samp_count < overall_set_count:
                sra_sc_jxs.remove(jx)
                uncategorized_jxs.update([jx])

        dev_removals = set()
        for jx, samp_count in dev_samp_counts.items():
            if samp_count < overall_set_count:
                dev_removals.update([jx])
                uncategorized_jxs.update([jx])
        for dev_set in dev_sets:
            dev_set.difference_update(dev_removals)

    jx_sets = [
        sra_adult_jxs, sra_dev_jxs, sra_sc_jxs, sra_emb_all, sra_fet, sra_zyg,
        sra_oo, sra_plc, sra_emb_ect, sra_emb_emb, sra_emb_late, sra_emb_mes,
        sra_emb_myo, uncategorized_jxs
    ]
    return jx_sets


def load_cancer_SRA_jxs(sra_can):
    txt_path = os.path.join(sra_can, '*rawresults*.txt')
    sra_can_files = glob.glob(txt_path)
    sra_can_jxs = set()
    for file in sra_can_files:
        with open(file) as lines:
            jxs = set(
                snaptron_results_to_jxs(
                    lines, min_sample_count=min_samps, min_read_count=min_reads
                )
            )
        sra_can_jxs = sra_can_jxs.union(jxs)

    return sra_can_jxs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Check TCGA junctions for developmental tissue evidence.'
    )
    parser.add_argument(
        '--snaptron-results',
        help='.txt file containing junction results from a previous snaptron '
             'search, or directory containing multiple of these.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output set membership files. subdirectories for '
             'developmental and unexplained junctions will be created within '
             'this directory.'
    )
    parser.add_argument(
        '--db-path', default='./',
        help='give the path for storing the created sql database.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--database-junction-directory', '-d',
        help='Specify a directory containing .csv files with junctions '
             'extracted via a jx_indexer query, each containing prevalence '
             'values for one cancer type.'
    )
    parser.add_argument(
        '--nongtex-junctions-directory', '-g',
        help='Specify a directory containing .csv files of cancer type '
             'specific junctions that are not present in GTEx.'
    )
    parser.add_argument(
        '--nonpaired-junctions-directory', '-p',
        help='Specify a directory containing .csv files of cancer type '
             'specific junctions that are not present in the paired normal '
             'tissue type(s) on GTEx.'
    )
    parser.add_argument(
        '--gtf-file',
        help='gtf file containing GENCODE annotation.'
    )
    parser.add_argument(
        '--min-SRA-sample-top_x', default=1, type=int,
        help='provide the minimum number of samples a junction must occur in '
             'in a given SRA sample type cohort to top_x as "in" the cohort.'
    )
    parser.add_argument(
        '--min-SRA-read-top_x', default=1, type=int,
        help='provide the minimum read coverage sum a junction must have in a '
             'given SRA sample type cohort to top_x as "in" the cohort.'
    )
    parser.add_argument(
        '--single-read-jx-json', required=True,
        help='Give the json file containing the list of single-read TCGA jxs, '
             'created by select_1-read_jxs.py.'
    )
    parser.add_argument(
        '--cancer-sra-directory',
        help='Specify the directory containing snaptron results for junctions '
             'from SRA cancer samples.'
    )
    parser.add_argument(
        '--cancer-gene-census',
        help='Provide the file with COSMIC cancer gene census data, '
             'cancer_gene_census.tsv'
    )
    parser.add_argument(
        '--oncokb-cancer-genes',
        help='Provide the file containing the onco-KB cancer gene list.'
    )
    parser.add_argument(
        '--min-overall-set-top_x', default=1, type=int,
        help='provide the minimum number of samples a junction must occur in '
             'in a broad SRA category (e.g. adult, stem cell, devlopmental) '
             'for membership in that category.'
    )

    args = parser.parse_args()
    snap_results = args.snaptron_results
    out_path = args.output_path
    db_path = args.db_path
    log_mode = args.log_level
    jx_dir = args.database_junction_directory
    nongtex_dir = args.nongtex_junctions_directory
    nonpair_dir = args.nonpaired_junctions_directory
    gtf_path = args.gtf_file
    min_samps = args.min_SRA_sample_count
    min_reads = args.min_SRA_read_count
    single_read_file = args.single_read_jx_json
    sra_can = args.cancer_sra_directory
    cancer_census = args.cancer_gene_census
    oncokb = args.oncokb_cancer_genes
    overall_set_count = args.min_overall_set_count

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'piechart_annotation_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    dev_dir = os.path.join(out_path, 'developmental')
    os.makedirs(dev_dir, exist_ok=True)

    unexpl_dir = os.path.join(out_path, 'unexplained')
    os.makedirs(unexpl_dir, exist_ok=True)

    # load junctions with only 1 read in TCGA
    with open(single_read_file) as recovered_data:
        single_read_jxs = json.load(recovered_data)

    single_read_jxs = set(single_read_jxs)

    # create connection to the database and pull junction-id info
    try:
        db_name = os.path.join(db_path, 'new_jx_index.db')
    except sql.OperationalError:
        print('If OperationalError is "unable to open database file":')
        print('make sure -d gives the PATH to the database directory,')
        print('not the database itself.')
        raise sql.OperationalError
    conn = sql.connect(db_name)
    index_db = conn.cursor()
    jx_id_map = create_jx_id_map(conn)

    # prepare annotation dictionaries from .gtf file
    coding_regions = gtf_to_cds(gtf_path)
    logging.info('coding regions discovered')
    CDS_interval_tree = cds_to_tree(coding_regions)
    logging.info('CDS tree created')
    antisense_interval_tree = cds_to_antisense(coding_regions)
    logging.info('antisense region tree created')
    id_name_dict = make_id_name_dict(gtf_path)
    jx_annotations = extract_splice_sites(gtf_path)

    # collect cancer genes
    cancer_loci = set()
    cancer_loci.update(collect_cancer_loci(
            cancer_census, symbol_column='Gene Symbol', separator=','
    ))
    cancer_loci.update(collect_cancer_loci(
        oncokb, symbol_column='Hugo Symbol', separator='\t'
    ))

    # Load SRA junctions
    jx_sets = load_non_cancer_SRA_jxs(
        snap_results, min_samps, min_reads, overall_set_count
    )
    sra_adult_jxs, sra_dev_jxs, sra_sc_jxs, sra_emb_all, sra_fet, sra_zyg, \
    sra_oo, sra_plc, sra_emb_ect, sra_emb_emb, sra_emb_late, sra_emb_mes, \
    sra_emb_myo, uncategorized_jxs = jx_sets

    sra_can_jxs = load_cancer_SRA_jxs(sra_can)

    for cancer in _TCGA_CANCER_TYPES:
        annotated_file = os.path.join(
            out_path,
            '{}_piechart_annotation_min_{}_samples.csv'
            ''.format(cancer, overall_set_count)
        )
        if os.path.exists(annotated_file):
            continue

        logging.info('starting {}'.format(cancer))
        all_jxs_name = '{}_all_jxs*.csv'.format(cancer)
        file = glob.glob(os.path.join(jx_dir, all_jxs_name))[0]

        jx_df = jx_df_from_file(
            file, 0.0, 1.0, chunk_it=True, glob_form=all_jxs_name,
            sample=False, top_x=False, drop_ann=False
        )
        jx_df['gencode'] = jx_df.annotation.apply(lambda x: x == 3).astype(int)

        # Load and annotate non-GTEx jxs
        nongtex_name = '{}_all_neojxs*non_GTEx*.txt'.format(cancer)
        nongtex_file = glob.glob(os.path.join(nongtex_dir, nongtex_name))[0]
        nongtex_df = pd.read_table(nongtex_file, sep=',').fillna(0)
        if 'jx' not in nongtex_df.columns.values:
            nongtex_df['jx'] = nongtex_df.jx_id.apply(lambda x: jx_id_map[x])
        nongtex_jxs = set(nongtex_df['jx'].tolist())
        jx_df['gtex'] = jx_df.jx.apply(
            lambda x: x not in nongtex_jxs
        ).astype(int)

        # Load and annotate non-tissue-matched jxs
        nonpair_name = '{}_all_neojxs*non_paired_normal*.txt'.format(cancer)
        nonpair_file = glob.glob(os.path.join(nonpair_dir, nonpair_name))[0]
        nonpair_df = pd.read_table(nonpair_file, sep=',').fillna(0)
        if 'jx' not in nonpair_df.columns.values:
            nonpair_df['jx'] = nonpair_df.jx_id.apply(lambda x: jx_id_map[x])
        nonpair_jxs = set(nonpair_df['jx'].tolist())
        jx_df['paired'] = jx_df.jx.apply(
            lambda x: x not in nonpair_jxs
        ).astype(int)

        # Annotate all SRA cell types
        jx_df['sra_stemcells'] = jx_df.jx.apply(
            lambda x: x in sra_sc_jxs
        ).astype(int)
        jx_df['sra_developmental'] = jx_df.jx.apply(
            lambda x: x in sra_dev_jxs
        ).astype(int)
        jx_df['sra_adult'] = jx_df.jx.apply(
            lambda x: x in sra_adult_jxs
        ).astype(int)
        jx_df['sra_embryo_all'] = jx_df.jx.apply(
            lambda x: x in sra_emb_all
        ).astype(int)
        jx_df['sra_embryo_ectoderm'] = jx_df.jx.apply(
            lambda x: x in sra_emb_ect
        ).astype(int)
        jx_df['sra_embryo_embryo'] = jx_df.jx.apply(
            lambda x: x in sra_emb_emb
        ).astype(int)
        jx_df['sra_embryo_lateembryo'] = jx_df.jx.apply(
            lambda x: x in sra_emb_late
        ).astype(int)
        jx_df['sra_embryo_mesenchyme'] = jx_df.jx.apply(
            lambda x: x in sra_emb_mes
        ).astype(int)
        jx_df['sra_embryo_myoblast'] = jx_df.jx.apply(
            lambda x: x in sra_emb_myo
        ).astype(int)
        jx_df['sra_neonate_fetal'] = jx_df.jx.apply(
            lambda x: x in sra_fet
        ).astype(int)
        jx_df['sra_zygote'] = jx_df.jx.apply(
            lambda x: x in sra_zyg
        ).astype(int)
        jx_df['sra_oocyte'] = jx_df.jx.apply(
            lambda x: x in sra_oo
        ).astype(int)
        jx_df['sra_placenta'] = jx_df.jx.apply(
            lambda x: x in sra_plc
        ).astype(int)

        # Remove single-read TCGA jxs
        poss_unexpl = jx_df[
            (jx_df.gtex == 0) &
            (jx_df.sra_stemcells == 0) &
            (jx_df.sra_adult == 0) &
            (jx_df.sra_developmental == 0)
        ]
        logging.info('init unexplained: {}'.format(len(poss_unexpl)))
        init_unexpl_count = len(poss_unexpl)
        poss_unexpl = poss_unexpl[~poss_unexpl.jx.isin(sra_can_jxs)]
        logging.info('removing sra cancer jxs: {}'.format(len(poss_unexpl)))
        poss_unexpl = poss_unexpl[~poss_unexpl.jx.isin(uncategorized_jxs)]
        logging.info(
            'removing uncategorized sra jxs: {}'.format(len(poss_unexpl))
        )
        jxs_to_kill = set(
            poss_unexpl[poss_unexpl.jx.isin(single_read_jxs)]['jx'].tolist()
        )
        logging.info(
            '{} of these have 1 TCGA read and will be killed:'
            ''.format(len(jxs_to_kill))
        )
        jx_df = jx_df[~jx_df.jx.isin(jxs_to_kill)]

        # Add other annotations
        jx_df['coding_regions'] = jx_df.jx.apply(
            lambda x: jx_gene_overlap(x, CDS_interval_tree, id_name_dict)
        )
        jx_df['antisense_regions'] = jx_df.jx.apply(
            lambda x: jx_gene_overlap(x, antisense_interval_tree, id_name_dict)
        )
        jx_df['sense'] = (jx_df.coding_regions != ';').astype(int)
        jx_df['antisense'] = (jx_df.antisense_regions != ';').astype(int)
        jx_df['cancer_locus'] = jx_df.coding_regions.apply(
            lambda x: check_cancer_loci(x, cancer_loci)
        )

        with open(annotated_file, 'w') as output:
            jx_df.to_csv(output, index=False)

        # select and save developmental jxs
        developmental_df = jx_df[
            (jx_df.gtex == 0) &
            (jx_df.sra_adult == 0) &
            (jx_df.sra_stemcells == 0) &
            (jx_df.sra_developmental == 1)

        ]
        dev_file = os.path.join(
            dev_dir,
            '{}_piechart_annotation_developmental_min_{}_samples.csv'
            ''.format(cancer, overall_set_count)
        )
        with open(dev_file, 'w') as output:
            developmental_df.to_csv(output, index=False)

        # select and save unexplained jxs
        unexplained_df = jx_df[
            (jx_df.gtex == 0) &
            (jx_df.sra_stemcells == 0) &
            (jx_df.sra_adult == 0) &
            (jx_df.sra_developmental == 0)
        ]
        logging.info('{} final unexplained jxs'.format(len(unexplained_df)))
        logging.info(
            'this is {}% of the originals'
            ''.format(len(unexplained_df) / init_unexpl_count)
        )
        unexpl_file = os.path.join(
            unexpl_dir,
            '{}_piechart_annotation_unexplained_min_{}_samples.csv'
            ''.format(cancer, overall_set_count)
        )
        with open(unexpl_file, 'w') as output:
            unexplained_df.to_csv(output, index=False)
