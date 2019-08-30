#!/usr/bin/env python3

"""
create_snaptron_queries.py

Python 3 code for creating snaptron batch query files and a slurm script for
running the downloads serially.

"""
import argparse
import csv
from datetime import datetime
import glob
import logging
import os
import sys
try:
    from utilities.utilities import _CHR_REGIONS, _SRA_ABBR
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _CHR_REGIONS, _SRA_ABBR


def collect_expts(sra_file, id_file, name_tag, out_path):
    """

    :param sra_file:
    :param id_file:
    :return:
    """
    expt_list = []
    not_in_recount = []
    sample_ids = {}
    expt_file = os.path.join(
        out_path, 'metaSRA-runs_{}_recount_exptlist.txt'.format(name_tag)
    )
    with open(id_file) as ids:
        csv_reader = csv.reader(ids, delimiter='\t')
        for line in csv_reader:
            sample_ids[line[2]] = line[0]
    with open(sra_file) as sra_array, open(expt_file, 'w') as output:
        next(sra_array)
        csv_reader = csv.reader(sra_array, delimiter=',')
        for i, experiment in enumerate(csv_reader):
            try:
                samp_id = str(sample_ids[experiment[5]])
                expt_list.append(samp_id)
                output.write('{}\n'.format(samp_id))
            except KeyError:
                not_in_recount.append(experiment[5])
    logging.info(
        '{} SRA expts collected:'.format(len(expt_list))
    )
    logging.info(
        '{} SRA expts without recount ids skipped'.format(len(not_in_recount))
    )
    return expt_list


def create_single_chrom_batch_query(expt_list, name_tag, util_path, proj_path,
                                    query_dir, dl_directory='snaptron_dl',
                                    slurm=False):
    """

    :param expt_list:
    :param outpath:
    :param name_tag:
    :param now:
    :return:
    """
    expts = ','.join(expt_list)
    script_file = os.path.join(
        proj_path, query_dir,
        '{}_indiv_chrom_batch_query_script.sh'.format(name_tag)
    )
    final_results = '{}_rawresults.txt'.format(name_tag)
    with open(script_file, 'w') as output:
        output.write('#!/usr/bin/bash\n')
        if slurm:
            output.write('#SBATCH --cpus-per-task=1\n')
            output.write('#SBATCH --mem-per-cpu=6G\n')
            output.write('#SBATCH --partition=PARTITION\n')
            output.write('#SBATCH -t 35:59:59\n\n')
        output.write('utilitypath={}\n'.format(util_path))
        output.write('projectpath={}\n\n'.format(proj_path))

    for region in _CHR_REGIONS:
        chr = region.split(':')[0]
        query_name = '{}_{}_snaptron_batch_query.tsv'.format(name_tag, chr)
        temp_snaptron_results = (
            'temp_snaptron_results_{}_{}.txt'.format(name_tag, chr)
        )
        query_file = os.path.join(proj_path, query_dir, query_name)
        with open(query_file, 'w') as output:
            output.write('region\teither\tsamples\n')
            output.write('{}\t2\t{}\n'.format(region, expts))

        with open(script_file, 'a') as output:
            output.write('echo "{}"\n'.format(query_name))
            output.write('n=0\n')
            output.write('until [ $n -ge 5 ]\ndo\n')
            output.write(
                'python ${up}/snaptron-experiments/client/query_snaptron.py '
                '--bulk-query-file ${pp}/{qd}/{qn} --bulk-query-stdout '
                '> ${pp}/{qd}/{tr} && break\n'.format(
                    up='{utilitypath}', pp='{projectpath}', qn=query_name,
                    dd=dl_directory, qd=query_dir, tr=temp_snaptron_results
                )
            )
            output.write('n=$[$n+1]\n')
            output.write('echo "Failure: sleeping 15 seconds"\n')
            output.write('sleep 15\n')
            output.write('done\n')
            if chr == 'chr1':
                output.write(
                    'head -1 ${pp}/{qd}/{tr} > ${pp}/{dd}/{fr}\n'
                    ''.format(
                        pp='{projectpath}', dd=dl_directory, qd=query_dir,
                        tr=temp_snaptron_results, fr=final_results
                     )
                )
            output.write(
                'tail -n +2 -q ${pp}/{qd}/{tr} >> ${pp}/{dd}/{fr}\n\n'
                ''.format(
                    pp='{projectpath}', dd=dl_directory, qd=query_dir,
                    tr=temp_snaptron_results, fr=final_results
                )
            )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Check TCGA junctions for developmental tissue evidence.'
    )
    parser.add_argument(
        '--MetaSRA-directory', '-s',
        help='directory containing multiple .csv files with SRA accession '
             'numbers from meta-SRA searches. Scripts will be created to '
             'download junctions from snaptron for the samples in each file, '
             'if possible.'
    )
    parser.add_argument(
        '--sample-id-file', '-i', required=True,
        help='File for cross referencing sample IDs with project and run '
             'accession numbers, e.g. the sample_ids.tsv from recount2.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output files: sampled spots, aligned junction '
             'reads, and SRA numbers with their p-values.'
    )
    parser.add_argument(
        '--log-level', '-l', default='INFO', choices=['INFO'],
        help='choose what logging mode to run (only INFO currently supported)'
    )
    parser.add_argument(
        '--project-path', '-p', default='./',
        help='If other than the current working directory, give the path '
             'to the main directory where the snaptron downloaded results '
             'should be based.'
    )
    parser.add_argument(
        '--utilities-path', '-u', default='./',
        help='If other than the current working directory, give the path '
             'to the directory where the query_snaptron tool can be accessed.'
    )

    args = parser.parse_args()
    sra_in = args.MetaSRA_directory
    id_file = args.sample_id_file
    out_path = args.output_path
    log_mode = args.log_level
    proj_path = args.project_path
    util_path = args.utilities_path

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'snaptron_batch_query_prep_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    proj_path = os.path.abspath(proj_path)
    util_path = os.path.abspath(util_path)

    can_result_dir = 'SRA_cancer_rawresults'
    os.makedirs(os.path.join(proj_path, can_result_dir), exist_ok=True)
    can_exptlist_path = os.path.join(proj_path, 'SRA_cancer_exptlists')
    os.makedirs(can_exptlist_path, exist_ok=True)

    noncan_result_dir = 'SRA_noncancer_rawresults'
    os.makedirs(os.path.join(proj_path, noncan_result_dir), exist_ok=True)
    noncan_exptlist_path = os.path.join(proj_path, 'SRA_noncancer_exptlists')
    os.makedirs(noncan_exptlist_path, exist_ok=True)

    query_dir = 'snaptron_batch_queries'
    os.makedirs(os.path.join(proj_path, query_dir), exist_ok=True)

    sra_files_path = os.path.join(sra_in, 'metaSRA-runs*.csv')
    sra_files = glob.glob(sra_files_path)
    for sra in sra_files:
        logging.info('starting snaptron prep for samples in {}'.format(sra))
        name_tag = os.path.basename(sra).split('.')[0]
        try:
            name_tag = name_tag.split("metaSRA-runs_")[1]
        except IndexError:
            name_tag = 'no_name'

        logging.info('base name is: {}'.format(name_tag))

        if name_tag in _SRA_ABBR:
            expt_path = noncan_exptlist_path
            result_dir = noncan_result_dir
        else:
            expt_path = can_exptlist_path
            result_dir = can_result_dir

        expt_list = collect_expts(sra, id_file, name_tag, expt_path)
        if not expt_list:
            continue

        create_single_chrom_batch_query(
            expt_list, name_tag, util_path, proj_path, query_dir,
            dl_directory=result_dir
        )
