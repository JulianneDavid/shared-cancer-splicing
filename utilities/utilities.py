from collections import defaultdict
import csv
import glob
import gzip
from intervaltree import IntervalTree
import json
import logging
from math import ceil
from math import floor
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import mmh3
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import re
import requests
from scipy import stats
import seaborn as sns; sns.set(color_codes=True)
import subprocess as sp
import time

_ID = '_IDs'
_COV = '_Covs'
_PER = '_Sample_Percents'
_MED_COV = '_Median_Coverage'
_STG = '_Stages'

_JX_SAMP_TABLE = 'jx_sample_map'
_JX_ANN_TABLE = 'jx_annotation_map'
_PHEN_TABLE = 'sample_phenotype_map'

_ID_TYPE_MAP = {
    'TCGA': 'TCGA_ID', 'file': 'TCGA_file_ID', 'recount': 'recount_ID',
    'case': 'TCGA_case_ID'
}

_ONT_PREFIXES = ('EFO', 'DOID', 'CL', 'UBERON', 'CVCL')

_TRANSLATION_DICT = {
    "GCT": "A", "GCG": "A", "GCA": "A", "GCC": "A", "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G", "ATT": "I", "ATC": "I", "ATA": "I", "CTT": "L",
    "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "TTT": "F", "TTC": "F", "TGG": "W", "TAT": "Y", "TAC": "Y",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "CAT": "H", "CAC": "H",
    "AAA": "K", "AAG": "K", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGT": "C", "TGC": "C", "ATG": "M", "AAT": "N", "AAC": "N", "CAA": "Q",
    "CAG": "Q", "TAA": "*", "TAG": "*", "TGA": "*"
    }

_SUPPORTED_NUCLEOTIDES = ['G', 'C', 'A', 'T']

_START_CODON = 'ATG'

_CHR_REGIONS = [
    'chr1:1-248956422', 'chr2:1-242193529', 'chr3:1-198295559',
    'chr4:1-190214555', 'chr5:1-181538259', 'chr6:1-170805979',
    'chr7:1-159345973', 'chr8:1-145138636', 'chr9:1-138394717',
    'chr10:1-133797422', 'chr11:1-135086622', 'chr12:1-133275309',
    'chr13:1-114364328', 'chr14:1-107043718', 'chr15:1-101991189',
    'chr16:1-90338345', 'chr17:1-83257441', 'chr18:1-80373285',
    'chr19:1-58617616', 'chr20:1-64444167', 'chr21:1-46709983',
    'chr22:1-50818468', 'chrX:1-156040895', 'chrY:1-57227415', 'chrM:1-16569'
]

_TCGA_ABBR = {
    'Acute_Myeloid_Leukemia': 'LAML',
    'Adrenocortical_Carcinoma': 'ACC',
    'Bladder_Urothelial_Carcinoma': 'BLCA',
    'Brain_Lower_Grade_Glioma': 'LGG',
    'Breast_Invasive_Carcinoma': 'BRCA',
    'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma': 'CESC',
    'Cholangiocarcinoma': 'CHOL',
    'Colon_Adenocarcinoma': 'COAD',
    'Esophageal_Carcinoma': 'ESCA',
    'Glioblastoma_Multiforme': 'GBM',
    'Head_and_Neck_Squamous_Cell_Carcinoma': 'HNSC',
    'Kidney_Chromophobe': 'KICH',
    'Kidney_Renal_Clear_Cell_Carcinoma': 'KIRC',
    'Kidney_Renal_Papillary_Cell_Carcinoma': 'KIRP',
    'Liver_Hepatocellular_Carcinoma': 'LIHC',
    'Lung_Adenocarcinoma': 'LUAD',
    'Lung_Squamous_Cell_Carcinoma': 'LUSC',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma': 'DLBC',
    'Mesothelioma': 'MESO',
    'Ovarian_Serous_Cystadenocarcinoma': 'OV',
    'Pancreatic_Adenocarcinoma': 'PAAD',
    'Pheochromocytoma_and_Paraganglioma': 'PCPG',
    'Prostate_Adenocarcinoma': 'PRAD',
    'Rectum_Adenocarcinoma': 'READ',
    'Sarcoma': 'SARC',
    'Skin_Cutaneous_Melanoma': 'SKCM',
    'Stomach_Adenocarcinoma': 'STAD',
    'Testicular_Germ_Cell_Tumors': 'TGCT',
    'Thymoma': 'THYM',
    'Thyroid_Carcinoma': 'THCA',
    'Uterine_Carcinosarcoma': 'UCS',
    'Uterine_Corpus_Endometrial_Carcinoma': 'UCEC',
    'Uveal_Melanoma': 'UVM'
}

_PRIMARY_ABBR = {
    'Acute_Myeloid_Leukemia': 'LAML', 'Adrenocortical_Carcinoma': 'ACC',
    'Astrocytoma': 'AC', 'Bladder_Urothelial_Carcinoma': 'BLCA',
    'Breast_Invasive_Carcinoma': 'BRCA', 'Cervical_Adenosquamous': 'CASC',
    'Cervical_Squamous_Cell_Carcinoma': 'CSC', 'Cholangiocarcinoma': 'CHOL',
    'Colon_Adenocarcinoma': 'COAD', 'Dedifferentiated_Liposarcoma': 'DDLS',
    'Desmoid_Tumor': 'DT', 'Endocervical_Adenocarcinoma': 'ECAD',
    'Esophagus_Adenocarcinoma': 'ESAD',
    'Esophagus_Squamous_Cell_Carcinoma': 'ESSC',
    'Glioblastoma_Multiforme': 'GBM',
    'Head_and_Neck_Squamous_Cell_Carcinoma': 'HNSC',
    'Kidney_Chromophobe': 'KICH', 'Kidney_Renal_Clear_Cell_Carcinoma': 'KIRC',
    'Kidney_Renal_Papillary_Cell_Carcinoma': 'KIRP', 'Leiomyosarcoma': 'LMS',
    'Liver_Hepatocellular_Carcinoma': 'LIHC', 'Lung_Adenocarcinoma': 'LUAD',
    'Lung_Squamous_Cell_Carcinoma': 'LUSC',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma': 'DLBC',
    'Malignant_Peripheral_Nerve_Sheath_Tumors': 'MPNT',
    'Mesothelioma': 'MESO', 'Myxofibrosarcoma': 'MFS',
    'Oligoastrocytoma': 'OAC', 'Oligodendroglioma': 'ODG',
    'Ovarian_Serous_Cystadenocarcinoma': 'OV',
    'Pancreatic_Adenocarcinoma': 'PAAD', 'Paraganglioma': 'PGG',
    'Pheochromocytoma': 'PCHC', 'Prostate_Adenocarcinoma': 'PRAD',
    'Rectum_Adenocarcinoma': 'READ', 'Skin_Cutaneous_Melanoma': 'SKCM',
    'Stomach_Adenocarcinoma': 'STAD', 'Synovial_Sarcoma': 'SYNS',
    'Testicular_Germ_Cell_Tumors': 'TGCT', 'Thymoma': 'THYM',
    'Thyroid_Carcinoma': 'THCA',
    'Undifferentiated_Pleomorphic_Sarcoma': 'UPLS',
    'Uterine_Carcinosarcoma': 'UCS',
    'Uterine_Corpus_Endometrial_Carcinoma': 'UCEC', 'Uveal_Melanoma': 'UVM'
}

_TCGA_CANCER_TYPES = _TCGA_ABBR.keys()
_CANCER_TYPES_PRIMARY = _PRIMARY_ABBR.keys()
_CONSTANT_CANCER_TYPES = list(
    set(_CANCER_TYPES_PRIMARY).intersection(_TCGA_CANCER_TYPES)
)
_ALL_CANCERS = list(
    set(_CANCER_TYPES_PRIMARY).union(set(_TCGA_CANCER_TYPES))
)
_ALL_ABBR = _TCGA_ABBR.copy()
_ALL_ABBR.update(_PRIMARY_ABBR)

_ABBR_TO_CAN = {value: key for key, value in _ALL_ABBR.items()}

_FONT_COLORS = {
    'THCA': 'black', 'GBM': 'white', 'LGG': 'black',
    'PAAD': 'white', 'MESO': 'white', 'READ': 'black',
    'SARC': 'white', 'COAD': 'black', 'BRCA': 'white',
    'ACC': 'white', 'LUAD': 'black', 'CHOL': 'white',
    'SKCM': 'black', 'CESC': 'black', 'PCPG': 'black',
    'HNSC': 'black', 'UCEC': 'black', 'PRAD': 'white',
    'KIRP': 'black', 'KIRC': 'black', 'LUSC': 'white',
    'BLCA': 'black', 'DLBC': 'white', 'UCS': 'black',
    'THYM': 'black', 'UVM': 'white', 'KICH': 'black',
    'TGCT': 'white', 'STAD': 'white', 'ESCA': 'white',
    'LIHC': 'black', 'OV': 'white', 'LAML': 'black',
    'pan cancer': 'black',
    'SKCM-to-normal\nsamplewise shared junctions': 'black',
    'MESO-to-normal\nsamplewise shared junctions': 'white',
    'pan cancer\njunction counts per group (#)': 'black',
}

_CANCER_COLORS = {
    'THCA': 'dandelion', 'GBM': 'medium purple', 'LGG': 'light mauve',
    'PAAD': 'blue grey', 'MESO': 'royal purple', 'READ': 'very light blue',
    'SARC': 'teal', 'COAD': 'light blue', 'BRCA': 'cerise',
    'ACC': 'yellow brown', 'LUAD': 'pale lilac', 'CHOL': 'marine blue',
    'SKCM': 'yellowish green', 'CESC': 'pale orange', 'PCPG': 'gold',
    'HNSC': 'sage green', 'UCEC': 'light peach', 'PRAD': 'dark red',
    'KIRP': 'rosy pink', 'KIRC': 'light rose', 'LUSC': 'deep lilac',
    'BLCA': 'light pink', 'DLBC': 'cobalt', 'UCS': 'tangerine',
    'THYM': 'taupe', 'UVM': 'emerald green', 'KICH': 'light red',
    'TGCT': 'scarlet', 'STAD': 'azure', 'ESCA': 'nice blue',
    'LIHC': 'light blue grey', 'OV': 'dull orange', 'LAML': 'melon',
    'pan cancer': 'off white',
    'SKCM-to-normal\nsamplewise shared junctions': 'yellowish green',
    'MESO-to-normal\nsamplewise shared junctions': 'royal purple',
    'pan cancer\njunction counts per group (#)': 'off white',

}

_TISSUE_TYPES = [
    'Adipose_Tissue', 'Adrenal_Gland', 'Bladder', 'Blood', 'Blood_Vessel',
    'Brain', 'Breast', 'Cervix_Uteri', 'Colon', 'Esophagus', 'Fallopian_Tube',
    'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pancreas',
    'Pituitary', 'Prostate', 'Salivary_Gland', 'Skin', 'Small_Intestine',
    'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina'
    ]

_ALL_TISSUE_CHOICES = _TISSUE_TYPES + _ALL_CANCERS

_ACCEPTABLE_CHROMS = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY'
}

_MATCHED_NORMALS = {
    'Acute_Myeloid_Leukemia': ['Blood'],
    'Adrenocortical_Carcinoma': ['Adrenal_Gland'],
    'Bladder_Urothelial_Carcinoma': ['Bladder'],
    'Brain_Lower_Grade_Glioma': ['Brain'],
    'Breast_Invasive_Carcinoma': ['Breast'],
    'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma':
        ['Cervix_Uteri'],
    'Colon_Adenocarcinoma': ['Colon'],
    'Esophageal_Carcinoma': ['Esophagus'],
    'Glioblastoma_Multiforme': ['Brain'],
    'Kidney_Chromophobe': ['Kidney'],
    'Kidney_Renal_Clear_Cell_Carcinoma': ['Kidney'],
    'Kidney_Renal_Papillary_Cell_Carcinoma': ['Kidney'],
    'Liver_Hepatocellular_Carcinoma': ['Liver'],
    'Lung_Adenocarcinoma': ['Lung'],
    'Lung_Squamous_Cell_Carcinoma': ['Lung'],
    'Mesothelioma': ['Lung'],
    'Ovarian_Serous_Cystadenocarcinoma': ['Ovary', 'Fallopian_Tube'],
    'Pancreatic_Adenocarcinoma': ['Pancreas'],
    'Pheochromocytoma_and_Paraganglioma': ['Adrenal_Gland', 'Nerve'],
    'Prostate_Adenocarcinoma': ['Prostate'],
    'Rectum_Adenocarcinoma': ['Colon'],
    'Sarcoma': ['Adipose_Tissue', 'Muscle'],
    'Skin_Cutaneous_Melanoma': ['Skin'],
    'Stomach_Adenocarcinoma': ['Stomach'],
    'Testicular_Germ_Cell_Tumors': ['Testis'],
    'Thyroid_Carcinoma': ['Thyroid'],
    'Uterine_Carcinosarcoma': ['Uterus'],
    'Uterine_Corpus_Endometrial_Carcinoma': ['Uterus'],
    'Head_and_Neck_Squamous_Cell_Carcinoma': ['Skin'],
    'Cholangiocarcinoma': ['Liver'],
    'Thymoma': [],
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma': [],
    'Uveal_Melanoma': []
}

_SRA_FET = {
    'neonate_cellline': 'nn_cl', 'neonate_primarycell': 'nn_pc',
    'neonate_tissue': 'nn_tis',
}

_SRA_ZYG = {
    'zygote_primarycell': 'zyg_pc',
}

_SRA_OO = {
    'oocyte_primarycell': 'ooc_pc'
}

_SRA_PLC = {
    'placenta_cellline': 'plc_cl', 'placenta_primarycell': 'plc_pc',
    'placenta_tissue': 'plc_tis'
}

_SRA_EMB_ECT = {
    'ectoderm_cellline': 'ect_cl'
}

_SRA_EMB_EMB = {
    'embryo_cellline': 'emb_cl', 'embryo_primarycell': 'emb_pc',
    'embryo_stemcells': 'emb_sc', 'embryo_tissue': 'emb_tis'
}

_SRA_EMB_LATE = {
    'lateembryo_cellline': 'le_cl', 'lateembryo_primarycell': 'le_pc',
    'lateembryo_tissue': 'le_tis'
}

_SRA_EMB_MES = {
    'mesenchyme_primarycells': 'mes_pc'
}

_SRA_EMB_MYO = {
    'myoblast_cellline': 'myb_cl', 'myoblast_primarycell': 'myb_pc'
}

_SRA_EMB_ALL = _SRA_EMB_EMB.copy()
_SRA_EMB_ALL.update(_SRA_EMB_ECT)
_SRA_EMB_ALL.update(_SRA_EMB_LATE)
_SRA_EMB_ALL.update(_SRA_EMB_MES)
_SRA_EMB_ALL.update(_SRA_EMB_MYO)

_SRA_DEV = _SRA_EMB_ALL.copy()
_SRA_DEV.update(_SRA_FET)
_SRA_DEV.update(_SRA_OO)
_SRA_DEV.update(_SRA_PLC)
_SRA_DEV.update(_SRA_ZYG)

_SRA_STEMCELLS = {
    'inducedpluripotentstemcell_cellline': 'ips_cl',
    'pluripotentstemcell_cellline': 'pps_cl',
    'pluripotentstemcell_stemcells': 'pps_sc',
    'somaticstemcell_stemcells': 'ssc_sc',
    'mesenchymalstemcell_stemcells': 'msc_sc',
    'hematopoietic_stemcells': 'hpc_sc',
    'glialcell_stemcells': 'gc_sc',
    'astrocyte_stemcells': 'ast_sc',
    'epithelialcell_stemcells': 'ec_sc',
    'fibroblast_stemcells': 'fb_sc',
    'mesenchyme_stemcells': 'mes_sc'

}

_SRA_ADULT = {
    'melanocyte_cellline': 'mel_cl', 'aorta_cellline': 'aor_cl',
    'aorta_tissue': 'aor_tis', 'astrocyte_cellline': 'ast_cl',
    'astrocyte_primarycell': 'ast_pc',
    'biliarytree_all': 'bt_all', 'bone_cellline': 'bone_cl',
    'bone_primarycell': 'bone_pc', 'bone_tissue': 'bone_tis',
    'epithelialcell_cellline': 'ec_cl',
    'epithelialcell_primarycell': 'ec_pc',
    'eye_cellline': 'eye_cl', 'eye_primarycell': 'eye_pc',
    'eye_tissue': 'eye_tis', 'fallopiantube_cellline': 'ft_cl',
    'fallopiantube_tissue': 'ft_tis', 'fibroblast_cellline': 'fb_cl',
    'fibroblast_primarycell': 'fb_pc',
    'gallbladder_primarycell': 'gb_pc', 'glialcell_cellline': 'gc_cl',
    'glialcell_primarycell': 'gc_pc', 'hematopoietic_cellline': 'hpc_cl',
    'hematopoietic_primarycell': 'hpc_pc', 'hepatocyte_cellline': 'hep_cl',
    'hepatocyte_primarycell': 'hep_pc', 'isletoflangerhans_cellline': 'ilh_cl',
    'isletoflangerhans_primarycell': 'ilh_pc', 'leukocyte_cellline': 'lk_cl',
    'leukocyte_primarycell': 'lk_pc', 'lymphocyte_cellline': 'lym_cl',
    'lymphocyte_primarycell': 'lym_pc',
    'macrophage_cellline': 'mph_cl', 'macrophage_primarycell': 'mph_pc',
    'melanocyte_primarycell': 'mel_pc', 'mesothelium_all': 'meso_all',
    'myeloidcell_cellline': 'myl_cl',
    'myeloidcell_primarycell': 'myl_pc',
    'oligodendrocyte_primarycell': 'odc_pc', 'platelet_cellline': 'plt_cl',
    'platelet_primarycell': 'plt_pc', 'thymus_primarycell': 'thym_pc',
    'thymus_tissue': 'thym_tis'
}

_SRA_ABBR = _SRA_ADULT.copy()
_SRA_ABBR.update(_SRA_STEMCELLS)
_SRA_ABBR.update(_SRA_DEV)

_FULL_HEATMAP_COLORS = _CANCER_COLORS.copy()
_FULL_HEATMAP_COLORS.update({
    'CSC': 'pale orange', 'ECAD': 'pale orange', 'CASC': 'pale orange',
    'PCHC': 'gold', 'PGG': 'gold',
    'LMS': 'teal', 'UPLS': 'teal', 'DT': 'teal', 'SYNS': 'teal', 'MFS': 'teal',
    'MPNT': 'teal',
    'AC': 'light mauve', 'OAC': 'light mauve', 'ODG': 'light mauve',
    'ESSC': 'nice blue', 'ESAD': 'nice blue'
})
_FULL_HEATMAP_COLORS.update({
    abbr: 'pale grey' for abbr in _SRA_ADULT.values()
})
_FULL_HEATMAP_COLORS.update({
    abbr: 'grey' for abbr in _SRA_DEV.values()
})
_FULL_HEATMAP_COLORS.update({
    abbr: 'light grey' for abbr in _SRA_STEMCELLS.values()
})

#
# def jx_to_motif(jx, reference_fasta, samtools_path):
#     """Given an input junction and reference genome, transcribes RNA sequence.
#
#     Input:
#     jx: a junction in 'chr_;left;right;strand' format, with 0-based fully
#         closed coordinates, (string)
#     reference_fasta: reference genome fasta file, previously sorted and
#         indexed by samtools (path to fasta file, string)
#     samtools_path: to be called from subprocess to collect the sequence (path
#         to samtools executable, string)
#
#     Returns a left-to-right nucleotide sequence on either side of the aberrant
#     junction sufficiently long to generate the desired protein sequence.
#     """
#     chrom, left, right, strand = jx.split(';')
#     if chrom not in _ACCEPTABLE_CHROMS:
#         return ''
#
#     shifted_right = str(int(right) + 2)
#
#     left_start = str(int(left) + 1)
#     left_stop = str(int(left) + 2)
#     left_range = chrom + ':' + left_start + '-' + left_stop
#     left_output = sp.check_output(
#         ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
#          '{}'.format(left_range)]
#     )
#     left_seq = ''.join(left_output.decode("utf-8").splitlines()[1:])
#
#     right_start = str(int(shifted_right) - 2)
#     right_stop = str(int(shifted_right) - 1)
#     right_range = chrom + ':' + right_start + '-' + right_stop
#     right_output = sp.check_output(
#         ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
#          '{}'.format(right_range)]
#     )
#     right_seq = ''.join(right_output.decode("utf-8").splitlines()[1:])
#
#     sequence = left_seq + right_seq
#     if strand == '-':
#         complement = str.maketrans("ATCG", "TAGC")
#         sequence = sequence.translate(complement)[::-1]
#
#     # print(sequence)
#     return sequence
#
#
# def jx_to_seq(chrom, left, right, reference_fasta, samtools_path,
#               protein_length, print_output, dna_len=0, split_seq=False,
#               strand=''):
#     """Given an input junction and reference genome, transcribes RNA sequence.
#
#     Input:
#     chrom: junction chromosome (string)
#     left: junction global left coordinate, 0-based, closed (string)
#     right: junction global right coordinate, 0-based, closed (string)
#     reference_fasta: reference genome fasta file, previously sorted and
#         indexed by samtools (path to fasta file, string)
#     samtools_path: to be called from subprocess to collect the sequence (path
#         to samtools executable, string)
#     protein_length: desired minimum length of returned protein sequence (int)
#     print_output: whether to print DNA sequence to std out (bool)
#     dna_len: to override desired minimum protein length with an exact length
#         for desired DNA sequence returned (int)
#     split_seq: to return two sequences, from the left and right sides of the
#         junction, separately instead of a joined RNA transcript sequence (bool)
#     strand: input required if returning split sequences, in the direction of
#         the transcript instead of global chromosome direction (string)
#
#     Returns a left-to-right nucleotide sequence on either side of the aberrant
#     junction sufficiently long to generate the desired protein sequence.
#     """
#     if dna_len:
#         num_bases = dna_len - 1
#     else:
#         num_bases = int(ceil(protein_length * 1.5) + 2)
#     left_start = str(int(left) - num_bases)
#     left_range = chrom + ':' + left_start + '-' + left
#     left_output = sp.check_output(
#         ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
#          '{}'.format(left_range)])
#
#     left_seq = ''.join(left_output.decode("utf-8").splitlines()[1:])
#
#     right = str(int(right) + 2)
#     right_stop = str(int(right) + num_bases)
#     right_range = chrom + ':' + right + '-' + right_stop
#     right_output = sp.check_output(
#         ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
#          '{}'.format(right_range)]
#     )
#     right_seq = ''.join(right_output.decode("utf-8").splitlines()[1:])
#     if print_output:
#         print('forward sequence is:')
#         sequence = left_seq + '|||' + right_seq
#         print(sequence)
#     if split_seq:
#         if strand == '+':
#             five_pr = left_seq
#             three_pr = right_seq
#         else:
#             complement = str.maketrans("ATCG", "TAGC")
#             five_pr = right_seq.translate(complement)[::-1]
#             three_pr = left_seq.translate(complement)[::-1]
#         return five_pr, three_pr
#     else:
#         sequence = left_seq + right_seq
#         return sequence
#
#
# def seq_to_protein(dna_sequence, strand, length, print_output=False):
#     """Translates an input RNA sequence to (an) amino acid sequence(s).
#
#     Input:
#     dna_sequence: string consisting of A, T, C, and G only (string)
#     strand: '+' or '-' (string)
#     length: target maximum length of protein to return; also determines the
#         minimum acceptable length of protein (int)
#     print_output: whether to print translated sequenc to standard out (boolean)
#
#     Determines the minimum acceptable length (0.75 of target length) of the
#     returned protein.  Translates each codon and appends to the protein in the
#     appropriate reading frame.
#
#     Returns a list of three strings; each is a protein sequence for one of
#     three possible reading frames.
#     """
#     min_pos = ceil(length * 0.75)
#     mid_pos = ceil(length * 1.5)
#     kill_seq = [False, False, False]
#     stop_seq = [False, False, False]
#     amino_acid_seq = ['' for _ in range(3)]
#     orig_seqs = ['' for _ in range(3)]
#     if strand == '-':
#         complement = str.maketrans("ATCG", "TAGC")
#         dna_sequence = dna_sequence.translate(complement)[::-1]
#
#     for k, char in enumerate(dna_sequence[:-2]):
#         window = dna_sequence[k:k+3]
#         frame = k % 3
#         amino_acid = codon_to_aa(window)
#         orig_seqs[frame] += amino_acid
#         if amino_acid == '*':
#             if k <= min_pos:
#                 amino_acid_seq[frame] = ''
#             elif min_pos < k < mid_pos:
#                 kill_seq[frame] = True
#             else:
#                 stop_seq[frame] = True
#         else:
#             if not stop_seq[frame]:
#                 amino_acid_seq[frame] += amino_acid
#
#     for frame, seq in enumerate(amino_acid_seq):
#         if kill_seq[frame] or (len(seq) < ceil(length/2)):
#             amino_acid_seq[frame] = ''
#     if print_output:
#         print('aa sequence is:')
#         print(amino_acid_seq)
#         print('original is:')
#         print(orig_seqs)
#     return amino_acid_seq
#
#
# def codon_to_aa(codon):
#     """Checks a codon for acceptable characters, and returns an amino acid.
#
#     Input:
#         codon: three-nucleotide amino acid sequence (string)
#
#     Returns the appropriate amino acid, stop signal, or not-supported signal.
#     """
#     for char in codon:
#         if char == '-':
#             return '-'
#         if char not in _SUPPORTED_NUCLEOTIDES:
#             return 'X'
#     return _TRANSLATION_DICT[codon]
#
#
# def write_fasta(amino_acid_seq, out_file, length, junction, gene_id,
#                 gene_name, samp_ids, false_positives, true_positives):
#     """Creates defline and writes defline and protein sequence to fasta file.
#
#     Input:
#     amino_acid_sequence: list containing amino acid sequences for three
#         reading frames for the junction (list of strings)
#     out_file: file name for newly written fasta file (string)
#     length: minimum length of protein sequence, for the defline (int)
#     junction: the junction corresponding to the protein sequence (str)
#     gene_id: id gene where these sequences arise, for the defline (string)
#     gene_name: name of gene where these sequences arise, for the defline (str)
#     samp_ids: if the sequences are associated with particular sample id(s),
#         they will be added to the defline (str)
#     false_positives: flag for the defline if the sequence is generated from a
#         known false positive junction (bool)
#     true_positives: flag for the defline if the sequene is generated from a
#         known true positive junction (bool)
#
#     Returns None.
#     """
#     with open(out_file, 'a') as protein_fasta:
#         jx = ';'.join(junction)
#         name = ''
#         id = ''
#         if gene_id:
#             id = 'gene_id:{}'.format(gene_id)
#         if gene_name:
#             name = 'gene_name:{}'.format(gene_name)
#         for k in range(3):
#             frame = k % 3
#             if false_positives:
#                 jx_id = '>FP_jx:{}'.format(jx)
#             elif true_positives:
#                 jx_id = '>TP_jx:{}'.format(jx)
#             else:
#                 jx_id = '>NOVEL_jx:{}'.format(jx)
#             if samp_ids:
#                 jx_id += ':' + samp_ids.strip('"')
#             min_len = 'min_length:{}'.format(length)
#             read_frame = 'frame:{}'.format(frame)
#             defline = ' '.join(filter(None, [jx_id, read_frame, min_len, id,
#                                              name]))
#             if amino_acid_seq[frame] != '':
#                 protein_fasta.write('{}\n'.format(defline))
#                 protein_fasta.write('{}\n'.format(amino_acid_seq[frame]))
#     return
#
#
# def calculate_motif_percents(jx_df, ref_fasta, samtools, samp_length=10000,
#                              out_path=None):
#     """Determines proportion of splice motifs occurring in dataframe jxs.
#
#     Input:
#         jx_df: (pandas DataFrame) dataframe containing junctions (and possibly
#             other information) from database query.
#         ref_fasta: (string) path to reference genome fasta file, previously
#             sorted and indexed by samtools
#         samtools: (string) path to samtools executable, to be called from
#             subprocess to collect the motif sequence
#         samp_length: (int) number of junctions to downsample for motif testing
#         out_path: (string) path to directory to store motif proportions file
#
#     Downsamples dataframe junctions, and collects splice motifs for the
#     sampled junctions. Groups and counts junctions by motifs; prints output to
#     a file.
#
#     Returns None
#     """
#     unique_jxs = pd.DataFrame(jx_df.jx.unique(), columns=['jx'])
#     unique_jxs['jx'] = unique_jxs.jx.astype(str)
#
#     if samp_length < len(jx_df):
#         unique_jxs = unique_jxs.sample(n=samp_length)
#     unique_jxs['motif'] = unique_jxs.jx.apply(
#         lambda x: jx_to_motif(x, ref_fasta, samtools)
#     )
#     num_jxs = len(unique_jxs)
#     motif_groups = unique_jxs.groupby(['motif'])['jx'].top_x()
#     motif_groups = pd.DataFrame(
#         {'motif': motif_groups.index,
#          'jx_count': motif_groups.values}
#     )
#     motif_groups.sort_values(by=['jx_count'], ascending=False, inplace=True)
#     motif_groups['percent'] = motif_groups.jx_count.apply(
#         lambda x: x / num_jxs
#     )
#     print(motif_groups.shape)
#     print(motif_groups.head())
#     print(motif_groups.tail())
#
#     if out_path:
#         out_file = os.path.join(
#             out_path, 'sampled_{}_neojx_motifs.csv'.format(samp_length)
#         )
#         with open(out_file, 'w') as output:
#             motif_groups.to_csv(output, index=False)
#     return
#
#
# def bed_file_to_jxs(jx_bed, regtools=False, shift=True, incl_counts=False):
#     """ Accepts bed file and returns list of junctions
#
#     Input:
#     jx_bed: bed file containing junctions to translate to junction list. (str)
#     regtools: the bed file was/wasn't generated by regtools: should be true if
#         the file was generated by regtools, false if not. (bool)
#     shift: the bed file was/wasn't generated by jx_indexer query: should be
#         False if the file was generated by jx_indexer, True if it is a "normal"
#         bed file that requires a shift to match jx_indexer 0-based closed
#         coordinate system (bool)
#     incl_counts: whether or not to include the number reads supporting each
#         junction; if True, each list item is a tuple of (jx, top_x) instead of
#         a string 'jx' (boolean)
#
#     Returns a list of unique junctions from the bed file.
#     """
#     all_junctions = []
#     with open(jx_bed) as bed:
#         jx_bed = csv.reader(bed, delimiter='\t')
#         for line in jx_bed:
#             chrom, left, right, strand = line[0], line[1], line[2], line[5]
#             if strand == '?':
#                 continue
#             if 'chr' not in chrom:
#                 chrom = 'chr' + chrom
#             if regtools:
#                 block_sizes = line[10].split(',')
#                 left = str(int(left) + int(block_sizes[0]))
#                 right = str(int(right) - int(block_sizes[1]) - 1)
#             elif shift:
#                 right = str(int(right) - 1)
#             jx = ';'.join([chrom, left, right, strand])
#             if incl_counts:
#                 top_x = line[4]
#                 all_junctions.append((jx, int(top_x)))
#             else:
#                 all_junctions.append(jx)
#     return list(set(all_junctions))
#

def create_jx_id_map(db_conn, jx_id_list=()):
    """Queries junction database to map junction recount ids to junctions

    Input:
        db_conn: (sqlite3 database connection) an open connection to the
            database
        jx_id_list: (list) specific junction ids to map to junction
            coordinates; if not provided, all junction recount ids are mapped.

    Queries database for recount ids and junction coordinates.  Loads query
    result into pandas dataframe. Creates dataseries map.

    Returns junction dataseries mapping recount IDs to jx coordinates.
    """
    where = ''
    if jx_id_list:
        jx_id_list = list(map(str, jx_id_list))
        ids = '"' + '", "'.join(jx_id_list) + '"'
        where = ' WHERE jx_id IN ({})'.format(ids)
    jx_map_query = (
        'SELECT jx, jx_id FROM {}'.format(_JX_ANN_TABLE) + where + ';'
    )
    jx_map_response = pd.read_sql_query(jx_map_query, db_conn)
    jx_dict = jx_map_response.set_index('jx_id')['jx']
    return jx_dict


def jx_df_from_file(jx_file, min_per, max_per, chunk_it, glob_form, sample,
                    top_x, pull_list=(), drop_ann=True, cancer=None):
    """Accepts directory containing junction prevalence files and returns df

    Input:
        jx_file: (str) file containing junction prevalences
        min_per: (float) junctions over this minimum percent only will be
            pulled.  (If a junction is pulled for one cancer type, it may have
            a lower prevalence in another cancer type - in which case some
            prevalence values may be below the min_per.)
        max_per: (float) junctions under this maximum percent prevalence only
            will be pulled.
        chunk_it: (bool) if this is called, files will be opened in chunks
            instead of all at once: memory saving.
        glob_form: (str) provides the random sample seed.
        sample: (int) if a value is given, this number of junctions will be
            downsampled using the glob form as the murmurhashed random seed.
        top_x: (int) if a value is given, this number of top-shared junctions
            only will be pulled.  (If a top-x and sample value are both given,
            the top-x will be pulled first and then the sample taken out of
            this set of top junctions.
        pull_list: (list) if this exists, min_per, max_per, top_x, and sample
            will be ignored: only junctions in the pull list will be returned.
        drop_ann: (bool) if selected, annotation flag will not appear
            in output or returned dataframe.

    Returns a pandas dataframe containing junctions & cancer type prevalences.
    """
    if chunk_it:
        max_chunk = 1e6
        jx_df = pd.DataFrame()
        chunks = pd.read_table(jx_file, sep=',', chunksize=max_chunk)
        for chunk in chunks:
            chunk = chunk.fillna(0)
            if pull_list:
                chunk = chunk.loc[chunk['jx'].isin(pull_list)]
            else:
                if cancer:
                    per_col = cancer + _PER
                    chunk = chunk[
                        (chunk[per_col] >= min_per) &
                        (chunk[per_col] <= max_per)
                    ]
                else:
                    if drop_ann and 'annotation' in chunk.columns.values:
                        drop_cols = ['jx', 'annotation']
                    else:
                        drop_cols = ['jx']
                    chunk = chunk[
                        (chunk.drop(drop_cols, axis=1) >= min_per).any(axis=1)
                    ]
                    chunk = chunk[
                        (chunk.drop(drop_cols, axis=1) <= max_per).any(axis=1)
                    ]
            jx_df = pd.concat([jx_df, chunk])
    else:
        jx_df = pd.read_table(jx_file, sep=',').fillna(0)
        if pull_list:
            jx_df = jx_df.loc[jx_df['jx'].isin(pull_list)]
        else:
            if drop_ann and 'annotation' in jx_df.columns.values:
                drop_cols = ['jx', 'annotation']
            else:
                drop_cols = ['jx']
            jx_df = jx_df[
                (jx_df.drop(drop_cols, axis=1) >= min_per).any(axis=1)
            ]
            jx_df = jx_df[
                (jx_df.drop(drop_cols, axis=1) <= max_per).any(axis=1)
            ]
    if not pull_list:
        if top_x and len(jx_df) > top_x:
            cancer = os.path.basename(jx_file).split(glob_form[:7])[0]
            perc = cancer + _PER
            jx_df = jx_df.sort_values(
                by=perc, axis=0, ascending=False
            ).head(top_x)
        if sample and len(jx_df) > sample:
            seed = abs(mmh3.hash(glob_form))
            jx_df = jx_df.sample(n=sample, random_state=seed)

    return jx_df


def jx_dir_to_df(jx_dir, min_per, max_per, chunk_it, print_full='', sample=0,
                 top_x=0, can_order=[], glob_form='_not-in-GTEx*jxs_*',
                 percentile=0, drop_annotation=True):
    """Accepts directory containing junction prevalence files and returns df

    Input:
        jx_dir: (str) directory containing junction prevalence files
        min_per: (float) junctions over this minimum percent only will be
            pulled.  (If a junction is pulled for one cancer type, it may have
            a lower prevalence in another cancer type - in which case some
            prevalence values may be below the min_per.)
        max_per: (float) junctions under this maximum percent prevalence only
            will be pulled.
        chunk_it: (bool) if this is called, files will be opened in chunks
            instead of all at once: memory saving.
        print_full: (str) if no string is given, the full data frame will not
            be written to a file.  If a string is give, the full data frame
            will be written to a file with the string value as its filename.
        sample: (int) if a value is given, this number of junctions will be
            downsampled using the glob form as the murmurhashed random seed.
        top_x: (int) if a value is given, this number of top-shared junctions
            only will be pulled.  (If a top-x and sample value are both given,
            the top-x will be pulled first and then the sample taken out of
            this set of top junctions.
        glob_form: (str) gives the form of prevalence file name to be pulled.
        percentile: (int) if a value is given, all other parameters - min_per,
            max_per, sample, and top_x - will be ignored.  Instead, the top
            percentile of shared junctions will be collected.

    Returns a pandas dataframe containing junctions & cancer type prevalences.
    """
    orig_top_x = top_x
    if percentile:
        min_per = 0.0
        max_per = 1.0
        sample = 0
        top_x = 0
    if can_order:
        jx_files = []
        for abbr in can_order:
            try:
                cancer = _ABBR_TO_CAN[abbr]
            except KeyError:
                cancer = abbr
            file_form = '{}'.format(cancer) + glob_form
            try:
                file = glob.glob(os.path.join(jx_dir, file_form))[0]
                jx_files.append(file)
            except IndexError:
                print('no prevalence file for {}: {}'.format(abbr, cancer))
    else:
        glob_form = '*' + glob_form
        jx_file_path = os.path.join(jx_dir, glob_form)
        jx_files = glob.glob(jx_file_path)

    jx_pull_list = set()
    for file in jx_files:
        temp_df = jx_df_from_file(
            file, min_per, max_per, chunk_it, glob_form, sample=sample,
            top_x=top_x
        )
        if percentile:
            num_jxs = floor(len(temp_df) * percentile / 100)
            cancer = os.path.basename(file).split('_not-in')[0]
            perc = cancer + _PER
            temp_df = temp_df.sort_values(
                by=perc, axis=0, ascending=False
            )
            perc_prev = temp_df.loc[num_jxs][perc]
            if drop_annotation and 'annotation' in temp_df.columns.values:
                drop_cols = ['jx', 'annotation']
            else:
                drop_cols = ['jx']
            temp_df = temp_df[
                (temp_df.drop(drop_cols, axis=1) > perc_prev).any(axis=1)
            ]
            temp_df = temp_df.head(orig_top_x)

        jx_pull_list.update(temp_df.jx.tolist())

    for file in jx_files:
        temp_df = jx_df_from_file(
            file, min_per, max_per, chunk_it, glob_form, sample=sample,
            top_x=top_x, pull_list=list(jx_pull_list)
        )
        try:
            shared_df = pd.merge(
                shared_df, temp_df, on=['jx', 'annotation'], how='outer'
            ).fillna(0)
        except UnboundLocalError:
            shared_df = temp_df

    if print_full:
        with open(os.path.join(jx_dir, 'full_merged_df.csv'), 'w') as output:
            shared_df.to_csv(output, index=False)

    if drop_annotation:
        try:
            shared_df = shared_df.drop(['annotation'], axis=1)
        except UnboundLocalError:
            shared_df = pd.DataFrame()
    return shared_df


def jx_file_to_df(jx_file, min_per, max_per, chunk_it):
    """Loads a pandas df with junction prevalences from a single .csv file

    Input:
    jx_file: .csv file generated by jx_indexer query containing cancer type
        prevalences for a one or more cancer types.
    min_per: filters prevalences to be over this minimum percentage only
    max_per: filters prevalences to be under this max percentage only
    chunk_it: splits the input file into pieces for memory-friendly loading

    Returns a pandas dataframe containing junctions and cancer type prevalences
    """
    if chunk_it:
        max_chunk = 1e6
        shared_df = pd.DataFrame()
        chunks = pd.read_table(jx_file, sep=',', chunksize=max_chunk)
        for chunk in chunks:
            chunk = chunk.fillna(0)
            chunk = chunk[
                (chunk.drop('jx', axis=1) >= min_per).any(axis=1)
            ]
            chunk = chunk[
                (chunk.drop('jx', axis=1) <= max_per).any(axis=1)
            ]
            shared_df = pd.concat([shared_df, chunk])
    else:
        shared_df = pd.read_table(jx_file, sep=',').fillna(0)
    try:
        shared_df['max'] = shared_df.drop(
            ['jx', 'annotation'], axis=1
        ).max(axis=1)
    except KeyError:
        shared_df['max'] = shared_df.drop(['jx'], axis=1).max(axis=1)
    shared_df = shared_df[shared_df['max'] >= min_per]
    shared_df = shared_df[shared_df['max'] <= max_per]
    try:
        shared_df = shared_df.drop(['max', 'annotation'], axis=1)
    except KeyError:
        shared_df = shared_df.drop(['max'], axis=1)

    return shared_df

#
# def load_ont_dict(json_file1, snap_file, jx_indices, shared_jx_set, testing,
#                   abbr, now, out_path):
#     """
#
#     :param json_file1:
#     :param snap_file:
#     :param jx_indices:
#     :param shared_jx_set:
#     :param testing:
#     :param abbr:
#     :param now:
#     :param out_path:
#     :return:
#     """
#     if json_file1:
#         with open(json_file1) as recover:
#             return json.load(recover)
#
#     start = time.time()
#     mutual_jxs = []
#     new_index = []
#     ont_dict = {}
#     ont_counts = {}
#     ont_to_name = {}
#     name_to_ont = {}
#     num_jxs = len(jx_indices)
#     with gzip.open(snap_file) as jx_ont:
#         for j, line in enumerate(jx_ont):
#             line = line.decode('utf-8').strip().split('\t')
#             chrom, left, right, strand = line[0], line[1], line[2], line[3]
#             left = str(int(left) - 1)
#             right = str(int(right) - 1)
#             jx = ';'.join([chrom, left, right, strand])
#             add_jx = False
#
#             if jx not in shared_jx_set:
#                 continue
#
#             ont_terms, counts = line[8], line[9]
#             if len(ont_terms) == 0:
#                 continue
#             ont_terms = ont_terms.split(',')
#             updated_terms = []
#             counts = counts.split(',')
#             for term, top_x in zip(ont_terms, counts):
#                 top_x = int(top_x)
#                 if top_x < 3:
#                     continue
#                 # if top_x < 3 or term in _COMMON_TERMS_SAMPLED:
#                 #     continue
#                 add_jx = True
#
#                 try:
#                     name = ont_to_name[term]
#                 except KeyError:
#                     if term.startswith(_ONT_PREFIXES):
#                         try:
#                             name = get_common_term(term)
#                         except requests.exceptions.ConnectionError:
#                             name = term
#                         ont_to_name[term] = name
#
#                         try:
#                             name_to_ont[name].append(term)
#                         except KeyError:
#                             name_to_ont[name] = [term]
#                     else:
#                         name = term
#
#                 if name not in updated_terms:
#                     try:
#                         ont_counts[name][0] += 1
#                     except KeyError:
#                         ont_counts[name] = [1]
#                     updated_terms.append(name)
#
#                 try:
#                     if top_x > ont_dict[name][jx_indices[jx]]:
#                         ont_dict[name][jx_indices[jx]] = top_x
#                 except KeyError:
#                     ont_dict[name] = np.zeros((num_jxs,), dtype=int)
#                     # ont_dict[term] = np.full((num_jxs,), np.nan)
#                     ont_dict[name][jx_indices[jx]] = top_x
#                 # try:
#                 #     ont_counts[term][0] += 1
#                 # except KeyError:
#                 #     ont_counts[term] = [1]
#
#             if add_jx:
#                 mutual_jxs.append(jx)
#                 new_index.append(jx_indices[jx])
#             if testing:
#                 if j % 200000 == 0 and j > 0:
#                     print('\nstill working:', j)
#                     print('time is:', time.time() - start)
#                     print(len(ont_dict))
#                     break
#             else:
#                 if j % 20000000 == 0:
#                     print('\nstill working:', j)
#                     print('time is:', time.time() - start)
#                     print(len(ont_dict))
#     name = '{}_json_ontdict_{}.txt'.format(abbr, now)
#     out_file = os.path.join(out_path, name)
#     dump_dict = {}
#     for key in ont_dict:
#         dump_dict[key] = ont_dict[key].tolist()
#     with open(out_file, 'w') as output:
#         dump_list = [
#             ont_counts, dump_dict, mutual_jxs, new_index, ont_to_name,
#             name_to_ont
#         ]
#         json.dump(dump_list, output)
#
#     return dump_list

#
# def get_common_term(ontology_term):
#     """Translates ontology term into common-language term
#
#     Input:
#         ontology_term: (string) DOID, OBERON, or other ontology term in the
#             form 'abbreviation:integer'
#
#     Queries the metasra site and parses the returned JSON file for the common
#     term. If the query is unsuccessful, substitutes an empty string.
#
#     Returns the metasra-translated common term for the input ontology term if
#     the query was successful, or an empty string if not.
#     """
#     query_site = (
#         'http://metasra.biostat.wisc.edu/api/v01/terms?id={}'
#         ''.format(ontology_term)
#     )
#     headers = {'Content-Type': 'application/json'}
#     response = requests.get(query_site, headers=headers)
#     if response.status_code == 200:
#         term_results = json.loads(response.content.decode('utf-8'))
#         try:
#             common_term = term_results['terms'][0]['name']
#         except IndexError:
#             common_term = ''
#     else:
#         print('status code is wrong!')
#         common_term = ''
#
#     return common_term


def collect_metasra_count(name_tag, list_dir):
    """Finds file containing list of experiments and counts the expts.

    Input:
        name_tag: (string) the sample type flag for the desired expt file.
        list_dir: (string) the path to the directory containing all experiment
            list files.

    Finds the correct experiment list file using glob. Counts the number of
    experiments (lines) in the file.

    Returns an integer: the number of lines/expts. in the file.
    """
    file_path = os.path.join(list_dir, '*{}*exptlist.txt'.format(name_tag))
    file = glob.glob(file_path)
    count = 0
    with open(file[0]) as expt_list:
        for line in expt_list:
            if line != '\n':
                count += 1
    return count


def snaptron_results_to_jxs(snaptron_result_lines, total_count=0,
                            min_sample_count=1, min_read_count=1):
    """

    :param snaptron_result_lines:
    :return:
    """
    expt_jxs = {}
    for line in snaptron_result_lines:
        items = line.split('\t')
        if not items[0].endswith('I'):
            continue

        samp_count = int(items[13])
        if samp_count < min_sample_count:
            continue

        read_count = int(items[14])
        if read_count < min_read_count:
            continue

        chrom, left, right, strand = items[2], items[3], items[4], items[6]
        left = str(int(left) - 1)
        right = str(int(right) - 1)
        jx = ';'.join([chrom, left, right, strand])
        if total_count:
            expt_jxs[jx] = samp_count / total_count
        else:
            expt_jxs[jx] = samp_count

    return expt_jxs


def collect_mutual_jxs(db_jx_df, expt_jxs, key, ret_type='percentage',
                       range=''):
    """

    :param db_jx_df:
    :param expt_jxs:
    :param key:
    :param ret_type:
    :param range:
    :return:
    """
    overlap = 0
    db_jxs = db_jx_df['jx'].tolist()
    mutual_jxs = set(db_jxs).intersection(set(expt_jxs))
    logging.info('\n{} mutual junctions'.format(len(mutual_jxs)))
    logging.info(
        'out of {} total shared database junctions'.format(len(db_jx_df))
    )
    if range:
        logging.info('with prevalences within range {}'.format(range))
    logging.info('and {} total {} junctions'.format(len(expt_jxs), key))
    if ret_type == 'percentage':
        overlap = (len(mutual_jxs) / len(db_jxs))
    elif ret_type == 'top_x':
        overlap = len(mutual_jxs)
    elif ret_type == 'ratio':
        overlap = '{}/{}'.format(len(mutual_jxs), len(db_jxs))
    elif ret_type == 'jxs':
        overlap = mutual_jxs
    return overlap


def logscale_heatmap(ax, df_to_plot, base_cmap=cm.Blues, min_val=0, max_val=0,
                     set_bad_to_white=False):
    """

    :param ax:
    :param df_to_plot:
    :param base_cmap:
    :param min_val:
    :param max_val:
    :param set_bad_to_white:
    :return:
    """
    df_to_plot.replace(to_replace=0, value=np.nan, inplace=True)
    if not min_val:
        min_val = df_to_plot.min().min()

    if not max_val:
        max_val = df_to_plot.max().max()

    plot_norms = matplotlib.colors.LogNorm(vmin=min_val, vmax=max_val)
    min_cmap_pt = 0.5
    max_cmap_pt = 1.0
    listcommand = 'trunc({name}, {a:.2f}, {b:.2f})'.format(
        name=base_cmap, a=min_cmap_pt, b=max_cmap_pt
    )
    hm_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        listcommand,
        base_cmap(np.linspace(min_cmap_pt, max_cmap_pt, base_cmap.N))
    )
    if set_bad_to_white:
        hm_cmap.set_bad(color='white')

    hm_mesh = ax.pcolormesh(df_to_plot, cmap=hm_cmap, norm=plot_norms)
    return hm_mesh


def heatmap_colorbar(plot_mesh, bar_label, log_plot, cbar_fontsize):
    """

    :param plot_mesh:
    :param bar_label:
    :param log_plot:
    :return:
    """
    cbar = plt.colorbar(
        plot_mesh, fraction=0.046, pad=0.04,
        format=ticker.FuncFormatter(
            lambda y, _: '{}%'.format('{:g}'.format(100*y))
        )
    )
    cbar.ax.text(
        0.55, 0.1, bar_label, fontsize=cbar_fontsize, rotation=90, ha='center',
        va='center', transform=cbar.ax.transAxes, color='black'
    )
    cbar_ticklabels = cbar.ax.get_yticklabels()
    num_ticks = len(cbar_ticklabels)
    new_labels = []
    for i, label in enumerate(cbar_ticklabels):
        text = label.get_text()
        if i == num_ticks - 1 or text.endswith('1'):
            new_labels.append(label)
        elif text.endswith('1%') or text.startswith('1'):
            new_labels.append(label)
        else:
            new_labels.append('')

    cbar.ax.set_yticklabels(new_labels)
    if log_plot:
        cbar.ax.tick_params(labelsize=5, width=1, length=4)
    return


def add_bottom_colorbar(ax, colorbar, plot_df):
    """

    :param ax:
    :param colorbar:
    :param plot_df:
    :return:
    """
    label_ranges = colorbar['ranges']
    colors = colorbar['colors']
    texts = colorbar['labels']
    try:
        bar_height = colorbar['height_percent']
    except KeyError:
        bar_height = '3.75%'
    # create axes next to plot
    divider = make_axes_locatable(ax)
    axb = divider.append_axes("bottom", bar_height, pad=0.0, sharex=ax)
    barkw = dict(
        color=colors, clip_on=False, align='edge', ec="k"
    )
    axb.bar(
        label_ranges[:, 0], 1,
        width=np.diff(label_ranges, axis=1).flatten(), **barkw
    )
    # set margins to zero again
    ax.margins(0)
    ax.tick_params(
        axis="both", bottom=0, left=0, labelbottom=0, labelleft=0
    )
    if texts:
        # Label the boxes
        textkw = dict(
            ha="center", va="center", fontsize=5.75, weight='heavy'
        )
        for (k, l), label in zip(label_ranges, texts):
            fontcolor = colorbar['fontcolors'][label]
            axb.text((k + l) / 2., 0.5, label, color=fontcolor, **textkw)

    axb.set(
        xticks=np.arange(0.5, plot_df.shape[1], 1),
        xticklabels=plot_df.columns.values.tolist(),
        yticklabels=[]
    )
    plt.setp(
        axb.xaxis.get_majorticklabels(), rotation=90, fontsize=5
    )
    return


def masked_double_heatmap(df_to_plot, cols_to_mask, fig_file, colorbar={},
                          masked_cmap=cm.Blues, other_cmap=cm.Greens,
                          title='', vline_pos=[], log_plot=True,
                          other_cbar_label='', size=None):
    """

    :param df_to_plot:
    :param cols_to_mask:
    :param fig_file:
    :param colorbar:
    :param masked_cmap:
    :param other_cmap:
    :param title:
    :param vline_pos:
    :param log_plot:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    if size:
        plt.rcParams['figure.figsize'] = size
        cbar_fontsize = min(size)
    else:
        cbar_fontsize = 8

    fig, ax = plt.subplots()

    white_df = df_to_plot.mask(df_to_plot > 0)
    masked = df_to_plot.mask(
        ((df_to_plot == df_to_plot) | df_to_plot.isnull())
        & (~df_to_plot.columns.isin(cols_to_mask))
    )
    other = df_to_plot.mask(
        ((df_to_plot == df_to_plot) | df_to_plot.isnull())
        & (df_to_plot.columns.isin(cols_to_mask))
    )
    other_mesh = ax.pcolormesh(
        white_df, cmap='binary'
    )

    if log_plot:
        masked.replace(to_replace=0, value=np.nan, inplace=True)
        other.replace(to_replace=0, value=np.nan, inplace=True)
        vmin = min(masked.min().min(), other.min().min())
        vmax = max(masked.max().max(), other.max().max())
        other_mesh = logscale_heatmap(
            ax, other, base_cmap=other_cmap, min_val=vmin, max_val=vmax
        )
        masked_mesh = logscale_heatmap(
            ax, masked, base_cmap=masked_cmap, min_val=vmin, max_val=vmax
        )

    else:
        other_mesh = ax.pcolormesh(other, cmap=other_cmap)
        masked_mesh = ax.pcolormesh(masked, cmap=masked_cmap)

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, fontsize=5)
    if title:
        ax.set_title(title, y=1.07)

    if vline_pos:
        ax.vlines(vline_pos, *ax.get_ylim())

    bar_label = '                           TCGA cancer type prevalence'
    heatmap_colorbar(masked_mesh, bar_label, log_plot, cbar_fontsize)

    if cols_to_mask != df_to_plot.columns.values.tolist():
        bar_label = '                      ' + other_cbar_label
        heatmap_colorbar(other_mesh, bar_label, log_plot, cbar_fontsize)

    ax.set_xticks(np.arange(0, df_to_plot.shape[1], 1), minor=True)
    ax.grid(axis='x', ls='-', color='white', which='minor')

    if colorbar:
        add_bottom_colorbar(ax, colorbar, df_to_plot)
    else:
        ax.set(
            xticks=np.arange(0.5, df_to_plot.shape[1], 1),
            xticklabels=df_to_plot.columns.values.tolist(),
            yticklabels=[]
        )
        ax.tick_params(direction='out')

    fig = plt.gcf()
    fig.savefig(fig_file)
    fig.clf()
    plt.clf()
    return


def get_jx_prev_filename(jx_dir, cancer):
    """Tests jx filename options; returns correct one based on input directory

    Input:
        jx_dir: (string) path to directory containing junction prevalence files
            of interest
        cancer: (string) TCGA primary cancer type or cancer subtype; this is
            the cancer type of interest for which to collect the junction file.

    Junction file may contain non-GTEx jxs, developmental jxs, or GTEx- and
    SRA-unexplained jxs; each of these has a different filename format. Tests
    each filename format in the given directory.

    Returns the filename, the appropriate flag, and the glob form of the file.
    """
    glob_options = {
        'all': '{}_not-in-GTEx*jxs_*'.format(cancer),
        'unexpl': '{}_piechart_annotation_unexplained*.csv'.format(cancer),
        'dev': '{}_piechart_annotation_developmental*.csv'.format(cancer),
        'sets': '{}_piechart_annotation*.csv'.format(cancer),
        'annotated': '{}*gene_antisense_ann*.csv'.format(cancer)
    }

    for flag, glob_form in glob_options.items():
        try:
            prev_glob = os.path.join(jx_dir, glob_form)
            jx_file = glob.glob(prev_glob)[0]
            break
        except IndexError:
            jx_file = ''
            flag = ''
            prev_glob = ''

    return  jx_file, flag, prev_glob


def make_id_name_dict(gtf_file):
    """Creates dictionary mapping gene IDs to gene names

    Input: gtf_file containing gene names and gene IDs

    Parses each line of the gtf file, and adds gene name-ID pairs to the dict.

    Returns the created dictionary.
    """
    id_name_dict = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()
            try:
                values = line.split('\t')[-1]
            except ValueError:
                continue
            values_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')
            if 'gene_id' not in values_dict or 'gene_name' not in values_dict:
                continue
            if values_dict['gene_id'] not in id_name_dict:
                id_name_dict[values_dict['gene_id']] = values_dict['gene_name']
    return id_name_dict


def gtf_to_cds(gtf_file):
    """ References cds_dict to get coding region (gene) bounds
    Keys in the dictionary are transcript IDs, while entries are lists of
        relevant CDS/stop codon data
        Data: [chromosome, left, right, +/- strand]
    Writes cds_dict as a pickled dictionary
    gtf_file: input gtf file to process
    dictdir: path to directory to store pickled dicts

    NOTE: from neoepiscope transcript.py

    Return value: dictionary
    """
    cds_dict = defaultdict(list)
    # Parse GTF to obtain CDS/stop codon info
    with open(gtf_file) as f:
        for line in f:
            if line[0] != '#':
                tokens = line.strip().split('\t')
                if tokens[2] == 'gene' and 'protein_coding' in line:
                    gene_id = re.sub(r'.*gene_id \"([A-Z0-9._]+)\"[;].*',
                                     r'\1', tokens[8])
                    gene_type = re.sub(r'.*gene_type \"([a-z_]+)\"[;].*',
                                       r'\1', tokens[8])
                    if gene_type == 'protein_coding':
                        # Create new dictionary entry for new each gene
                        cds_dict[gene_id].append(
                            [tokens[0], int(tokens[3]), int(tokens[4]),
                             tokens[6]]
                        )
    return cds_dict


def cds_to_tree(cds_dict):
    """ Creates searchable tree of chromosome intervals from CDS dictionary

    Each chromosome is stored in the dictionary as an interval tree object
        Intervals are added for each CDS, with the associated transcript ID
        Assumes transcript is all on one chromosome - does not work for
            gene fusions
    Writes the searchable tree as a pickled dictionary
    cds_dict: CDS dictionary produced by gtf_to_cds()

    NOTE: from neoepiscope transcript.py

    Return value: searchable tree
    """
    searchable_tree = {}
    # Add genomic intervals to the tree for each transcript
    for gene_id in cds_dict:
        gene = cds_dict[gene_id]
        chrom = gene[0][0]
        # Add new entry for chromosome if not already encountered
        if chrom not in searchable_tree:
            searchable_tree[chrom] = {}
        # Add CDS interval to tree with transcript ID
        for cds in gene:
            left = cds[1]
            right = cds[2]
            strand = cds[3]
            # Interval coordinates are inclusive of left, exclusive of right
            if strand not in searchable_tree[chrom]:
                searchable_tree[chrom][strand] = IntervalTree()
            if right > left:
                searchable_tree[chrom][strand][left:right+1] = gene_id
    return searchable_tree


def cds_to_antisense(cds_dict):
    """ Creates searchable tree of chromosome intervals from CDS dictionary

    Each chromosome is stored in the dictionary as an interval tree object
        Intervals are added for each CDS, with the associated transcript ID
        Assumes transcript is all on one chromosome - does not work for
            gene fusions
    Writes the searchable tree as a pickled dictionary
    cds_dict: CDS dictionary produced by gtf_to_cds()

    NOTE: from neoepiscope transcript.py

    Return value: searchable tree
    """
    searchable_tree = {}
    strand_flip = {'+': '-', '-': '+'}
    # Add genomic intervals to the tree for each transcript
    for gene_id in cds_dict:
        gene = cds_dict[gene_id]
        chrom = gene[0][0]
        # Add new entry for chromosome if not already encountered
        if chrom not in searchable_tree:
            searchable_tree[chrom] = {}
        # Add CDS interval to tree with transcript ID
        for cds in gene:
            left = cds[1]
            right = cds[2]
            orig_strand = cds[3]
            strand = strand_flip[orig_strand]
            # Interval coordinates are inclusive of left, exclusive of right
            if strand not in searchable_tree[chrom]:
                searchable_tree[chrom][strand] = IntervalTree()
            if right > left:
                searchable_tree[chrom][strand][left:right+1] = gene_id
    return searchable_tree


def jx_ends_in_cds(junction, cds_tree, id_name_dict):
    """Check found junctions for coding region overlap

    Input:
        junction information: chromosome, left and right boundaries of the
            junction, and strand
        CDS tree containing coding regions, created by gtf_to_cds
        dictionary mapping gene ids from .gtf file to gene names

    Checks to see whether either junction end overlaps coding regions.  If
    either or both do, collects gene ids and names for the overlapping genes.
    If both sides overlap, checks to see if any of the genes are the same.

    Returns eight entries comprising new column information for the junction
    database; columns contain the following information:
    1) binary: whether both ends of the junction overlap no known gene regions
    2) binary: whether one end of the junction only overlaps gene regions
    3) binary: whether two ends of the junction overlap different genes
    4) binary: whether two ends of the junction overlap the same gene
    5) comma-separated list of gene ids overlapped by the 5' junction end
    6) comma-separated list of gene ids overlapped by the 3' junction end
    7) comma-separated list of gene names overlapped by the 5' junction end
    8) comma-separated list of gene names overlapped by the 3' junction end
    """
    no_overlap = 0
    one_overlap = 0
    both_same = 0
    both_diff = 0
    left_genes = []
    right_genes = []
    left_names = []
    right_names = []
    chrom, left, right, strand = junction.split(';')
    left = int(left)
    right = int(right)
    try:
        jx_start = cds_tree[chrom][strand].overlaps(left)
        jx_stop = cds_tree[chrom][strand].overlaps(right)
    except KeyError:
        return (no_overlap, one_overlap, both_diff, both_same, left_genes,
                right_genes, left_names, right_names)

    if jx_start or jx_stop:
        for start_set in list(cds_tree[chrom][strand][left]):
            if start_set[2] not in left_genes:
                left_genes.append(start_set[2])
            name = id_name_dict[start_set[2]]
            if name not in left_names:
                left_names.append(name)
        for stop_set in list(cds_tree[chrom][strand][right]):
            if stop_set[2] not in right_genes:
                right_genes.append(stop_set[2])
            name = id_name_dict[stop_set[2]]
            if name not in right_names:
                right_names.append(name)
        if jx_start and jx_stop:
            num_same_genes = len(set(left_genes) & set(right_genes))
            if num_same_genes > 0:
                both_same = 1
            if ((len(right_genes) - num_same_genes > 0) or
                    (len(left_genes) - num_same_genes > 0)):
                both_diff = 1
        else:
            one_overlap = 1

    if strand == '+':
        fivepr_genes = ','.join(left_genes)
        threepr_genes = ','.join(right_genes)
        fivepr_names = ','.join(left_names)
        threepr_names = ','.join(right_names)
    else:
        fivepr_genes = ','.join(right_genes)
        threepr_genes = ','.join(left_genes)
        fivepr_names = ','.join(right_names)
        threepr_names = ','.join(left_names)

    return (no_overlap, one_overlap, both_diff, both_same, fivepr_genes,
            threepr_genes, fivepr_names, threepr_names)


def jx_gene_overlap(junction, cds_tree, id_name_dict):
    """Check found junctions for coding region overlap

    Input:
        junction information: chromosome, left and right boundaries of the
            junction, and strand
        CDS tree containing coding regions, created by gtf_to_cds
        dictionary mapping gene ids from .gtf file to gene names

    Checks to see whether either junction end overlaps coding regions.  If
    either or both do, collects gene ids and names for the overlapping genes.
    If both sides overlap, checks to see if any of the genes are the same.

    Returns a semicolon-joined string of gene names overlapping the 5' and 3'
    junction ends
    """
    left_genes = []
    right_genes = []
    left_names = []
    right_names = []
    chrom, left, right, strand = junction.split(';')
    left = int(left)
    right = int(right)
    try:
        jx_start = cds_tree[chrom][strand].overlaps(left)
        jx_stop = cds_tree[chrom][strand].overlaps(right)
    except KeyError:
        return ';'

    if jx_start or jx_stop:
        for start_set in list(cds_tree[chrom][strand][left]):
            if start_set[2] not in left_genes:
                left_genes.append(start_set[2])
            name = id_name_dict[start_set[2]]
            if name not in left_names:
                left_names.append(name)
        for stop_set in list(cds_tree[chrom][strand][right]):
            if stop_set[2] not in right_genes:
                right_genes.append(stop_set[2])
            name = id_name_dict[stop_set[2]]
            if name not in right_names:
                right_names.append(name)

    if strand == '+':
        fivepr_names = ','.join(left_names)
        threepr_names = ','.join(right_names)
    else:
        fivepr_names = ','.join(right_names)
        threepr_names = ','.join(left_names)

    return ';'.join([fivepr_names, threepr_names])


def check_annotations(junction, junction_dict):
    """Adds annotated splice junctions from .gtf file to the junction list.

    Input:
        junction: (string) 0-based closed coordinates; from database query.
        junction_dict: (dict) created by extract_splice_sites function.

    Junction column key:
        0 = neither junction side is annotated
        1 = one junction side is annotated
        2 = both junction sides are annotated, but not together
        3 = junction is fully annotated

    NOTE: currently, db junctions use 0-based closed coordinates and GENCODE
          1-based open: so the right end of a junction is 2 smaller than what
          is listed in GENCODE annotation. To accommodate this, all junction
          right ends have 2 added before being checked against the junction
          annotation dictionary.

    Returns an int corresponding to junction column key above.
    """
    tokens = junction.split(';')
    chrom, left, right, strand = tokens[0], tokens[1], tokens[2], tokens[3]
    right = str(int(right) + 2)
    junction = chrom + ';' + left + ';' + right + ';' + strand
    try:
        if junction in junction_dict[chrom][strand]['full']:
            annotated_col = 3
        else:
            annotated_col = 0
            if strand == '+':
                five_site = left
                three_site = right
            else:
                five_site = right
                three_site = left
            if five_site in junction_dict[chrom][strand]['fivepr']:
                annotated_col += 1
            if three_site in junction_dict[chrom][strand]['threepr']:
                annotated_col += 1
    except KeyError:
        annotated_col = 0
    return annotated_col


def extract_splice_sites(gtf_file):
    """Extracts splice site anns_same_same from .gtf file

    This function is a modified version of one that is part of HISAT.
    Copyright 2014, Daehwan Kim <infphilo@gmail.com>

    HISAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HISAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
    """
    genes = defaultdict(list)
    trans = {}

    annotations = {}

    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                item = line.split('\t')
                chrom, left, right, strand = item[0], item[3], item[4], item[6]
                feature, values = item[2], item[8]
            except ValueError:
                continue

            left, right = int(left), int(right)

            if feature != 'exon' or left >= right:
                continue

            values_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

            if 'gene_id' not in values_dict or \
                    'transcript_id' not in values_dict:
                continue

            transcript_id = values_dict['transcript_id']
            if transcript_id not in trans:
                trans[transcript_id] = [chrom, strand, [[left, right]]]
                genes[values_dict['gene_id']].append(transcript_id)
            else:
                trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, exons[i - 1][1], exons[i][0], strand))
    junctions = sorted(junctions)

    for chrom, left, right, strand in junctions:
        if chrom not in annotations:
            annotations[chrom] = {}
            annotations[chrom]['+'] = {'full': [], 'fivepr': [], 'threepr': []}
            annotations[chrom]['-'] = {'full': [], 'fivepr': [], 'threepr': []}
        annotations[chrom][strand]['full'].append(
            ';'.join([chrom, str(left), str(right), strand]))
        if strand == '+':
            five_site = str(left)
            three_site = str(right)
        else:
            five_site = str(right)
            three_site = str(left)
        if five_site not in annotations[chrom][strand]['fivepr']:
            annotations[chrom][strand]['fivepr'].append(five_site)
        if three_site not in annotations[chrom][strand]['threepr']:
            annotations[chrom][strand]['threepr'].append(three_site)
    return annotations


def grouped_boxplots_with_table(data_dict, plot_dict, fig_file,
                                fig_size=(3.0, 5.0), logscale=True,
                                y_label='cohort prevalence', percent=True,
                                right_lim_shift=2):
    """

    :param data_dict:
    :param out_path:
    :param now:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = fig_size
    sns.set_context("paper")
    sns.set_style("whitegrid")
    sns.set(style="ticks")

    f, (ax, ax2) = plt.subplots(
        nrows=2, ncols=1, gridspec_kw={'height_ratios': [3, 1]}
    )
    plt.sca(ax)
    light_cols = plot_dict['light colors']
    dark_cols = plot_dict['dark colors']
    num_boxes = len(light_cols)
    adjust_val = num_boxes + 1
    num_table_rows = 0
    for values in data_dict.values():
        num_table_rows = max(num_table_rows, len(values['table_data']))

    bp_pos = list(range(1, num_boxes + 1))

    columns = []
    col_cols = []
    label_locs = []
    labels = []
    table_vals = [[] for _ in range(num_table_rows)]
    curr_label_loc = sum(bp_pos) / len(bp_pos)

    for abbr, values in data_dict.items():
        columns.append(abbr)
        try:
            col_cols.append('xkcd:{}'.format(_CANCER_COLORS[abbr]))
        except KeyError:
            for can_abbr in _CANCER_COLORS.keys():
                if can_abbr in abbr:
                    col_cols.append('xkcd:{}'.format(_CANCER_COLORS[can_abbr]))
                    break

        table_data = values['table_data']
        for full_table_set, curr_table_val in zip(table_vals, table_data):
            try:
                full_table_set.append('{:,}'.format(curr_table_val))
            except ValueError:
                full_table_set.append(curr_table_val)

        data = values['data']

        boxes = plt.boxplot(
            data, positions=bp_pos, widths=0.6, patch_artist=True
        )
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=4.0, linestyle='none', linewidth=0.15,
                markeredgecolor=dark_cols[i]
            )
            plt.setp(
                med, color=dark_cols[i], linewidth=1.5
            )
            plt.setp(box, facecolor='white')

        label_locs.append(curr_label_loc)

        bp_pos = [x + adjust_val for x in bp_pos]
        curr_label_loc += adjust_val

    ax.set_xticklabels(labels)
    ax.set_xticks(label_locs)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.set_xlim(left=0, right=curr_label_loc-right_lim_shift)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('gray')
    ax.spines['right'].set_color('gray')
    ax.spines['left'].set_color('gray')

    if logscale:
        plt.yscale('log')

    plt.ylabel(y_label, fontsize=7)
    plt.setp(
        ax.xaxis.get_majorticklabels(), rotation=90, fontsize=5, color='black'
    )
    plt.setp(
        ax.yaxis.get_majorticklabels(), fontsize=7, color='black'
    )
    if percent:
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(
                lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
            )
        )
    else:
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
        )
    plt.xticks([])

    rows = plot_dict['row labels']
    whitefont_cols = []
    for i, abbr in enumerate(columns):
        try:
            font_color = _FONT_COLORS[abbr]
        except KeyError:
            for can_abbr in _FONT_COLORS.keys():
                if can_abbr in abbr:
                    font_color = _FONT_COLORS[can_abbr]
                    break

        if font_color == 'white':
            whitefont_cols.append(i)

    row_label_cols = plot_dict['row colors']
    the_table = ax.table(
        cellText=table_vals, rowLabels=rows, colLabels=columns, loc='bottom',
        cellLoc='center', colColours=col_cols, rowColours=row_label_cols
    )
    the_table.set_fontsize(5)

    for (row, col), cell in the_table.get_celld().items():
        if (row == 0):
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=5.5)
            )
            if len(columns[0]) > 5:
                cell.set_height(cell.get_height() * 1.5)
            if col in whitefont_cols:
                cell._text.set_color('white')
        if col == -1:
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=5.5)
            )
            cell._text.set_color(plot_dict['row font color'][row - 1])

    ax2.yaxis.grid(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.patch.set_alpha(0)

    plt.sca(ax2)
    plt.plot([])
    plt.xticks([], [])
    plt.yticks([], [])

    fig = plt.gcf()
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return
#
#
# def scatter_with_regression_and_size(scatter_df, cancer, out_path, now,
#                                      y_axis_label='median antisense coverage',
#                                      flag=''):
#     """
#
#     :param kdf_init:
#     :param dbdf_init:
#     :param cancer:
#     :param can_count:
#     :param out_path:
#     :return:
#     """
#     with open(os.path.join(out_path, 'initscatter.csv'), 'w') as output:
#         scatter_df.to_csv(output, index=False)
#
#     x = scatter_df['x']
#     y = scatter_df['y']
#     slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
#     line = slope * x + intercept
#     print('linear regression stats:')
#     print(slope, intercept, r_value, p_value, std_err)
#     line_label = (
#         'y = {}x + {}\nR = {}'
#         ''.format(round(slope, 3), round(intercept, 3), round(r_value, 2))
#     )
#
#     groups = scatter_df.groupby(['x', 'y', 'cancer'])
#     scatter_dict = {'x': [], 'y': [], 'top_x': [], 'color': []}
#     for group_index, group in groups:
#         scatter_dict['x'].append(group_index[0])
#         scatter_dict['y'].append(group_index[1])
#         # scatter_dict['top_x'].append(6 * (group.item.top_x() ** (1.0 / 1.5)))
#         scatter_dict['top_x'].append(group.item.top_x() * 6)
#         # scatter_dict['top_x'].append((math.sqrt(group.item.top_x()))**1.3)
#         # scatter_dict['top_x'].append(10 * math.sqrt(group.item.top_x()))
#         can_color = _CANCER_COLORS[_TCGA_ABBR[group_index[2]]]
#         scatter_dict['color'].append('xkcd:{}'.format(can_color))
#     scatter_df = pd.DataFrame(scatter_dict)
#
#     with open(os.path.join(out_path, 'init2scatter.csv'), 'w') as output:
#         scatter_df.to_csv(output, index=False)
#
#     # print(scatter_df[scatter_df['k_prev'] > 0.5].head(10))
#     num_jxs = len(scatter_df)
#     # ax_max = max(scatter_df.x.max(), scatter_df.y.max()) * 1.1
#     ax_max = max(scatter_df.x.max(), scatter_df.y.max()) * 1.1
#     # ax_min = 0.1 * ax_max
#
#     point_dict = {'fake_x': [], 'fake_y': [], 'color': []}
#     for color in scatter_dict['color']:
#         point_dict['fake_x'].append(-ax_max * 10)
#         point_dict['fake_y'].append(-ax_max * 10)
#         point_dict['color'].append(color)
#     fake_df = pd.DataFrame(point_dict)
#
#     plt.rcParams.update({'figure.autolayout': True})
#     plt.rcParams['figure.figsize'] = 5.0, 5.0
#     matplotlib.style.use('seaborn-whitegrid')
#     plt.scatter(
#         x='x', y='y', c='color', data=scatter_df, s='top_x',
#         edgecolors='xkcd:light grey', linewidths=0.05, label=None
#     )
#     # plt.scatter(
#     #     x='x', y='y', c='color', data=scatter_df, s='top_x', linewidths=0.05,
#     #     label=None
#     # )
#     # plt.scatter(
#     #     x='k_prev', y='rail_prev', c='color', data=scatter_df, s='top_x',
#     #     edgecolors='xkcd:light grey', linewidths=0.05, label='shared jx'
#     # )
#     legend_label = '{} sample'.format(_TCGA_ABBR[cancer])
#     plt.scatter(
#         x='fake_x', y='fake_y', c='color', data=fake_df, s=40,
#         edgecolors='xkcd:light grey', linewidth=0.01, label=legend_label
#     )
#     plt.plot(
#         x, line, c='xkcd:dark grey', linewidth=0.8, label=line_label
#     )
#     # plt.plot(
#     #     x, line, c='xkcd:grey', linewidth=1, label=line_label
#     # )
#     # plt.legend(fontsize='small', markerscale=.4, frameon=True);
#     plt.legend(fontsize='small', frameon=True)
#     # lgnd = plt.legend(fontsize='small', scatterpoints=1, frameon=True)
#     # lgnd.legendHandles[0]._legmarker.set_markersize(6)
#
#     # lgnd = plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
#     # lgnd.legendHandles[0]._sizes = [30]
#     # # lgnd.legendHandles[1]._sizes = [30]
#
#     fig_name = (
#         '{}_antisense_coverage_withregression_{}_{}.pdf'
#         ''.format(cancer, now, flag)
#     )
#     # try:
#     #     title = (
#     #         '{}: {} TCGA samples, {} mutual jxs'
#     #         ''.format(_TCGA_ABBR[cancer], can_count, num_jxs)
#     #     )
#     # except KeyError:
#     #     title = (
#     #         'All TCGA (except LAML): {} samples, {} mutual jxs'
#     #         ''.format(can_count, num_jxs)
#     #     )
#     # plt.title(title)
#     plt.xlabel('median sense coverage')
#     # plt.xticks(x_locs, x_labels)
#     # plt.setp(ax.get_xticklabels(), rotation=90, fontsize=6)
#     # plt.yscale('log')
#     # plt.yscale('symlog', linthresy=2500)
#     ax = plt.gca()
#     # ax.set_ylim(ymin=-ax_min, ymax=ax_max)
#     # ax.set_xlim(xmin=-ax_min, xmax=ax_max)
#     # ax.set_ylim(ymin=-ax_min)
#     # ax.set_xlim(xmin=-ax_min)
#     y_max = scatter_df.y.max() * 1.1
#     x_max = scatter_df.x.max() * 1.1
#     ax.set_ylim(ymin=-(0.05 * y_max), ymax=y_max)
#     ax.set_xlim(xmin=-(0.05 * x_max), xmax=x_max)
#
#     plt.ylabel(y_axis_label)
#     fig = plt.gcf()
#     fig_file = os.path.join(out_path, fig_name)
#     print(fig_file)
#     fig.savefig(fig_file)
#     plt.clf()
#     return
