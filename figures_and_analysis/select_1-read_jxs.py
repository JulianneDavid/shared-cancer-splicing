#!/usr/bin/env python3

"""
select_1-read_jxs.py

Python3 script for collecting junctions with only a single read in TCGA.

"""

import argparse
import csv
import json
import os
import time

class JXIndexError(Exception):
    pass


def collect_single_read_jxs(cov_file, bed_file):
    """Parses coverage and .bed files and returns list of single-read junctions

    :param cov_file:
    :param bed_file:

    return: list of single-read junctions.
    """
    single_read_jxs = []
    fill_int = time.time()
    start_time = fill_int
    with open(cov_file) as cov_file, open(bed_file) as bed:
        print('starting tcga junctions')
        jx_cov = csv.reader(cov_file, delimiter='\t')
        jx_bed = csv.reader(bed, delimiter='\t')
        for i, (line, jx_info) in enumerate(zip(jx_cov, jx_bed)):
            jx_id, ids, covs = (line[0], line[1], line[2])
            chrom, left, right, strand = (
                jx_info[0], jx_info[1], jx_info[2], jx_info[5]
            )
            if jx_id != jx_info[3].split('|')[0]:
                raise JXIndexError(
                    "jx ids at row {row} don't match\n"
                    "jx from cov file: {cj}\njx from bed file: {bj}\n"
                    "".format(row=i, cj=jx_id, bj=jx_info[3].split('|')[0])
                )

            try:
                if int(covs) == 1:
                    jx = ';'.join([chrom, left, right, strand])
                    single_read_jxs.append(jx)
            except ValueError:
                continue

            if (i % 1000000) == 0:
                fill_2000 = time.time()
                print('\n{}th entry'.format(i))
                print('intermediate fill time is', fill_2000 - fill_int)
                print('total fill time is', fill_2000 - start_time)
                print('number of single-read jxs:', len(single_read_jxs))
                fill_int = time.time()
    return single_read_jxs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='collect junctions with only 1 read in TCGA.'
    )
    parser.add_argument(
        '--tcga-jx-cov', '-C', required=True,
        help='File with junction coverages for TCGA.'
    )
    parser.add_argument(
        '--tcga-jx-bed', '-B', required=True,
        help='BED file for TCGA junction location information.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store output json file.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    tcga_cov = args.tcga_jx_cov
    tcga_bed = args.tcga_jx_bed

    single_read_jxs = collect_single_read_jxs(tcga_cov, tcga_bed)

    out_file = os.path.join(out_path, 'single-read-tcga-jxs.json')
    with open(out_file, 'w') as output:
        json.dump(single_read_jxs, output)
