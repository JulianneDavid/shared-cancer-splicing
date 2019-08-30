#!/usr/bin/env python3

"""
count_unique_SRA_expts.py

Python 3 code to determine total number of SRA experiments queried.

"""

import argparse
from datetime import datetime
import glob
import logging
import os
import sys
try:
    from utilities.utilities import _SRA_ADULT, _SRA_DEV, _SRA_STEMCELLS
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _SRA_ADULT, _SRA_DEV, _SRA_STEMCELLS


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count SRA experiments used.'
    )
    parser.add_argument(
        '--non-cancer-expt-lists', '-n',
        help='Give the directory with experiment list files for non-cancer '
             'SRA samples.'
    )
    parser.add_argument(
        '--cancer-expt-lists', '-c',
        help='Give the directory with experiment list files for cancer SRA '
             'samples.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output log file.'
    )

    args = parser.parse_args()
    noncan_exptlist_dir = args.non_cancer_expt_lists
    can_exptlist_dir = args.cancer_expt_lists
    log_mode = args.log_level
    out_path = args.output_path

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'SRA_sample_count_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    can_ids = set()
    non_can_ids = set()
    adult_ids = set()
    dev_ids = set()
    stem_cell_ids = set()
    all_sra_ids = set()

    noncan_path = os.path.join(noncan_exptlist_dir, '*_exptlist.txt')
    noncan_exptlists = glob.glob(noncan_path)

    can_path = os.path.join(can_exptlist_dir, '*_exptlist.txt')
    can_exptlists = glob.glob(can_path)

    for samp_type in noncan_exptlists:
        name_tag = os.path.basename(samp_type).split('.')[0]
        name_tag = name_tag.split('_recount_exptlist')[0]
        try:
            name_tag = name_tag.split('metaSRA-runs_')[1]
        except IndexError:
            pass

        with open(samp_type) as expt_list:
            for line in expt_list:
                if line != '\n':
                    all_sra_ids.update([line])
                    non_can_ids.update([line])
                    if name_tag in _SRA_ADULT:
                        adult_ids.update([line])
                    if name_tag in _SRA_DEV:
                        dev_ids.update([line])
                    if name_tag in _SRA_STEMCELLS:
                        stem_cell_ids.update([line])

    for samp_type in can_exptlists:
        with open(samp_type) as expt_list:
            for line in expt_list:
                if line != '\n':
                    all_sra_ids.update([line])
                    can_ids.update([line])

    logging.info('total number of SRA samples: {}'.format(len(all_sra_ids)))
    logging.info('total cancer samples: {}'.format(len(can_ids)))
    logging.info('total non-cancer samples: {}'.format(len(non_can_ids)))
    logging.info('adult non-cancer: {}'.format(len(adult_ids)))
    logging.info('developmental: {}'.format(len(dev_ids)))
    logging.info('stem cells: {}'.format(len(stem_cell_ids)))
