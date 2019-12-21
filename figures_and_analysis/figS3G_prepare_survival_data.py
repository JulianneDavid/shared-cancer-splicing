#!/usr/bin/env python3

"""
paper1_reviewerresponse_assemble_survival_data.py


"""

import argparse
import os
import pandas as pd
import sqlite3 as sql
import sys
try:
    from utilities.utilities import _JX_ANN_TABLE, _JX_SAMP_TABLE, _PHEN_TABLE
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _JX_ANN_TABLE, _JX_SAMP_TABLE, _PHEN_TABLE


def access_time(followup, death_time):
    try:
        return int(death_time)
    except ValueError:
        try:
            return int(followup)
        except ValueError:
            return 0


def check_msln_status(case_id, target_list):
    if case_id in target_list:
        return 1
    else:
        return 0


def check_death(death_time):
    try:
        int(death_time)
        return 1
    except ValueError:
        return 0


def query_db(db_path):
    try:
        db_name = os.path.join(db_path, 'new_jx_index.db')
        db_conn = sql.connect(db_name)
    except sql.OperationalError:
        print('If OperationalError is "unable to open database file": ')
        print('make sure -d gives the PATH to the database directory, ')
        print('not the database itself.')
        raise sql.OperationalError

    samp_pull_query = (
        'SELECT * FROM '
        '   (SELECT recount_id recid, coverage FROM {js} INNER JOIN {ja} '
        '    ON {ja}.jx_id == {js}.jx_id WHERE jx == "chr16;766903;768491;-") '
        'INNER JOIN {sp} ON recid == {sp}.recount_id '
        'WHERE primary_type == "Ovarian_Serous_Cystadenocarcinoma";'
        ''.format(
            ja=_JX_ANN_TABLE, js=_JX_SAMP_TABLE, sp=_PHEN_TABLE
        )
    )
    query_response = pd.read_sql_query(samp_pull_query, db_conn)
    return set(query_response['case_id'].tolist())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates file for survival response curve plot.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store survival data output.'
    )
    parser.add_argument(
        '--db-path', '-d', default='./',
        help='give the path for the created sql database.'
    )
    parser.add_argument(
        '--TCGA-phenotype-file', '-P',
        help='TCGA.tsv file from recount2 with sample phenotype information.'
    )

    args = parser.parse_args()
    phen_table = args.TCGA_phenotype_file
    out_path = args.output_path
    db_path = args.db_path

    MSLN_sampset = query_db(db_path)

    tcga_phen = pd.read_table(
        phen_table,
        usecols=['xml_days_to_last_followup',
                 'gdc_cases.diagnoses.days_to_death',
                 'gdc_cases.case_id',
                 'gdc_cases.project.name'],
        dtype=str
    )
    new_names = {
        'xml_days_to_last_followup': 'followup',
        'gdc_cases.diagnoses.days_to_death': 'death_time',
        'gdc_cases.case_id': 'case_id',
        'gdc_cases.project.name': 'cancer'
    }
    tcga_phen.rename(new_names, axis='columns', inplace=True)
    tcga_phen = tcga_phen[
        tcga_phen['cancer'] == 'Ovarian Serous Cystadenocarcinoma'
    ]
    tcga_phen.drop(['cancer'], axis=1, inplace=True)
    tcga_phen['with_MSLN_jx'] = tcga_phen.case_id.apply(
        lambda x: check_msln_status(x, MSLN_sampset)
    )
    tcga_phen['death'] = tcga_phen.death_time.apply(
        lambda x: check_death(x)
    )
    tcga_phen['time'] = tcga_phen.apply(
        lambda x: access_time(x['followup'], x['death_time']),
        axis=1
    )
    tcga_phen = tcga_phen[tcga_phen['time'] > 0]
    tcga_phen.drop(['followup', 'death_time', 'case_id'], axis=1, inplace=True)
    outfile = os.path.join(out_path, 'OV_mslnjx_survival_data.csv')
    with open(outfile, 'w') as output:
        tcga_phen.to_csv(outfile, index=False)
