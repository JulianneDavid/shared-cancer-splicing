import argparse


_help_intro = """jx_indexer.py has two modes: index and experiment.

In index mode, we create a database of RNA splicing junctions in cancer and
normal tissue, along with the samples in which they are found, listed under
their relevant tissue or cancer type, and the coverage for each sample.

In experiment mode, cancer junction files are generated for data processing.

"""


def help_formatter(prog):
    """ So formatter_class's max_help_position can be changed. """
    return argparse.HelpFormatter(prog, max_help_position=40)


parser = argparse.ArgumentParser(
    description=_help_intro,
    formatter_class=help_formatter
)
subparser = parser.add_subparsers(
    help='',
    dest='subparser_name'
)
index_parser = subparser.add_parser(
    'index',
    help='Creates index of normal and cancer tissue junctions.'
)
experiment_parser = subparser.add_parser(
    'experiment',
    help='run jx_indexer experiments: counts and collects non-paired-'
         'normal and non-GTEx neojxs.'
)

parser.add_argument(
    '--db-path', '-d', default='./',
    help='give the path for storing the created sql database.'
)
parser.add_argument(
    '--log-level', '-l',
    choices=['DEBUG', 'INFO', 'ERROR', 'WARNING', 'CRITICAL'],
    default='INFO',
    help='Choose what logging mode to run.'
)
parser.add_argument(
    '--testing', action='store_true',
    help='select this option for running test on mini mock database.'
)

index_parser.add_argument(
    '--tcga-jx-cov', '-C', required=True,
    help='File with junction coverages for TCGA.'
)
index_parser.add_argument(
    '--tcga-jx-bed', '-B', required=True,
    help='BED file for TCGA junction location information.'
)
index_parser.add_argument(
    '--tcga-phenotype', '-P', required=True,
    help='File with TCGA phenotype information.'
)
index_parser.add_argument(
    '--gtex-jx-cov', '-c', required=True,
    help='File with junction coverages for GTEx.'
)
index_parser.add_argument(
    '--gtex-jx-bed', '-b', required=True,
    help='BED for for GTEx junction location information.'
)
index_parser.add_argument(
    '--gtex-phenotype', '-p', required=True,
    help='File with GTEx phenotype information'
)
index_parser.add_argument(
    '--sample-id-file', '-s', required=True,
    help='File for cross referencing sample IDs with project and run '
         'accession numbers.'
)
index_parser.add_argument(
    '--gtf-file', '-g', required=True,
    help='gtf file containing GENCODE annotation.'
)
experiment_parser.add_argument(
    '--output-path', '-o', default='./',
    help='Give the path to store neojx experiment output.'
)
experiment_parser.add_argument(
    '--batch-number', type=int,
    help='For running queries for multiple cancer types in parallel via '
         'slurm: batch numbers correspond to cancer types per batch dict '
         'above. This should start at 1, not 0!'
)
