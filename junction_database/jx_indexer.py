#!/usr/bin/env python3

"""
jx_indexer.py
Python 3.6 code for indexing and scoring RNA junctions.
This software requires junction coverage and BED files, phenotype information
    for each sample, and an index to reference sample IDs with their project
    and run accession numbers.

These can be found at recount2: https://jhubiostatistics.shinyapps.io/recount/

GTEx information was downloaded from the GTEx tab: jx_bed and jx_cov files
    under the junctions column, and the link under "phenotypes".

TCGA information was obtained in the same way, from the TCGA tab.

The sample ID-accession number index can be found under Documentation; the link
    is given in the description of a Junction Raw Coverage File.

"""
from datetime import datetime
import os
import sqlite3 as sql


if __name__ == '__main__':
    from jx_parser import parser
    from importlib import import_module

    args = parser.parse_args()
    db_path = args.db_path
    log_mode = args.log_level
    testing = args.testing

    db_name = os.path.join(db_path, 'new_jx_index.db')
    if args.subparser_name == 'index' and os.path.exists(db_name):
        options = ['y', 'n']
        user_check = input(
            '\nThe database file already exists. Would you like to erase the '
            'current database file and begin generating the database from '
            'scratch? [y/n]\n(Note: the full database is 249GB and takes '
            'significant resources to generate. Deleting the database file is '
            'only recommended if a previous indexing run did not complete '
            'successfully.\nWould you like to delete the database file? '
            '[y/n]: '
        )
        while (user_check not in options):
            user_check = input(
                "Your response isn't recognized. Please enter 'y' if you "
                "would like to delete the existing database file and generate "
                "the database from scratch and 'n' otherwise: "
            )
        if user_check == 'n':
            print('Database will not be generated.  Exiting.')
            exit()

        else:
            print('Deleting existing database file.')
            os.remove(db_name)

    try:
        conn = sql.connect(db_name)
        index_db = conn.cursor()
    except sql.OperationalError:
        print('If OperationalError is "unable to open database file": ')
        print('make sure -d gives the PATH to the database directory, ')
        print('not the database itself.')
        raise sql.OperationalError

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    submain = import_module(args.subparser_name)
    submain.main(args, now, conn, index_db)
