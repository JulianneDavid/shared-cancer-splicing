from collections import defaultdict
import glob
from intervaltree import IntervalTree
import logging
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
import seaborn as sns; sns.set(color_codes=True)

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
    'neonate_cellline': 'neonate cellline',
    'neonate_primarycell': 'neonate primary cell',
    'neonate_tissue': 'neonate tissue',
}

_SRA_ZYG = {
    'zygote_primarycell': 'zygote primary cell',
}

_SRA_OO = {
    'oocyte_primarycell': 'oocyte primary cell'
}

_SRA_PLC = {
    'placenta_cellline': 'placenta cell line',
    'placenta_primarycell': 'placenta primary cell',
    'placenta_tissue': 'placenta tissue'
}

_SRA_EMB_ECT = {
    'ectoderm_cellline': 'ectoderm cell line'
}

_SRA_EMB_EMB = {
    'embryo_cellline': 'embryo cell line',
    'embryo_primarycell': 'embryo primary cell',
    'embryo_stemcells': 'embryo stem cells',
    'embryo_tissue': 'embryo tissue'
}

_SRA_EMB_LATE = {
    'lateembryo_cellline': 'late embryo cell line',
    'lateembryo_primarycell': 'late embryo primary cell',
    'lateembryo_tissue': 'late embryo tissue'
}

_SRA_EMB_MES = {
    'mesenchyme_primarycells': 'mesenchyme primary cells'
}

_SRA_EMB_MYO = {
    'myoblast_cellline': 'myoblast cell line',
    'myoblast_primarycell': 'myoblast primary cell'
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
    'inducedpluripotentstemcell_cellline': 'iPS cell line',
    'pluripotentstemcell_cellline': 'pluripotent stem cell line',
    'pluripotentstemcell_stemcells': 'pluripotent stem cells',
    'somaticstemcell_stemcells': 'somatic stem cells',
    'mesenchymalstemcell_stemcells': 'mesenchymal stem cells',
    'hematopoietic_stemcells': 'hematopoietic stem cell stem cells',
    'glialcell_stemcells': 'glial cell stem cells',
    'astrocyte_stemcells': 'astrocyte stem cells',
    'epithelialcell_stemcells': 'epithelial cell stem cells',
    'fibroblast_stemcells': 'fibroblast stem cells',
    'mesenchyme_stemcells': 'mesenchyme stem cells'
}

_SRA_ADULT = {
    'melanocyte_cellline': 'melanocyte cell line',
    'aorta_cellline': 'aorta cell line',
    'aorta_tissue': 'aorta tissue',
    'astrocyte_cellline': 'astrocyte cell line',
    'astrocyte_primarycell': 'astrocyte primary cell',
    'biliarytree_all': 'biliary tree',
    'bone_cellline': 'bone cell line',
    'bone_primarycell': 'bone primary cell',
    'bone_tissue': 'bone tissue',
    'epithelialcell_cellline': 'epithelial cell line',
    'epithelialcell_primarycell': 'epithelial primary cell',
    'eye_cellline': 'eye cell line',
    'eye_primarycell': 'eye primary cell',
    'eye_tissue': 'eye tissue',
    'fallopiantube_cellline': 'fallopian tube cell line',
    'fallopiantube_tissue': 'fallopian tube tissue',
    'fibroblast_cellline': 'fibroblast cell line',
    'fibroblast_primarycell': 'fibroblast primary cell',
    'gallbladder_primarycell': 'gall bladder primary cell',
    'glialcell_cellline': 'glial cell line',
    'glialcell_primarycell': 'glial primary cell',
    'hematopoietic_cellline': 'hematopoietic cell line',
    'hematopoietic_primarycell': 'hematopoietic primary cell',
    'hepatocyte_cellline': 'hepatocyte cell line',
    'hepatocyte_primarycell': 'hepatocyte primary cell',
    'isletoflangerhans_cellline': 'pancreatic islet cell line',
    'isletoflangerhans_primarycell': 'pancreatic islet primary cell',
    'leukocyte_cellline': 'leukocyte cell line',
    'leukocyte_primarycell': 'leukocyte primary cell',
    'lymphocyte_cellline': 'lymphocyte cell line',
    'lymphocyte_primarycell': 'lymphocyte primary cell',
    'macrophage_cellline': 'macrophage cell line',
    'macrophage_primarycell': 'macrophage primary cell',
    'melanocyte_primarycell': 'melanocyte primary cell',
    'mesothelium_all': 'mesothelium',
    'myeloidcell_cellline': 'myeloid cell line',
    'myeloidcell_primarycell': 'myeloid primary cell',
    'oligodendrocyte_primarycell': 'oligodendrocyte primary cell',
    'platelet_cellline': 'platelet cell line',
    'platelet_primarycell': 'platelet primary cell',
    'thymus_primarycell': 'thymus primary cell',
    'thymus_tissue': 'thymus tissue'
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


def heatmap_colorbar(plot_mesh, bar_label, log_plot, cbar_fontsize, fraction,
                     pad, cbar_labels):
    """

    :param plot_mesh:
    :param bar_label:
    :param log_plot:
    :return:
    """
    cbar = plt.colorbar(
        plot_mesh, fraction=fraction, pad=pad,
        format=ticker.FuncFormatter(
            lambda y, _: '{}%'.format('{:g}'.format(100*y))
        )
    )
    cbar.ax.text(
        0.55, 0.1, bar_label, fontsize=cbar_fontsize, rotation=90, ha='center',
        va='center', transform=cbar.ax.transAxes, color='black'
    )
    if not cbar_labels:
        cbar_ticklabels = cbar.ax.get_yticklabels()
        num_ticks = len(cbar_ticklabels)
        cbar_labels = []
        for i, label in enumerate(cbar_ticklabels):
            text = label.get_text()
            if i == num_ticks - 1 or text.endswith('1'):
                cbar_labels.append(label)
            elif text.endswith('1%') or text.startswith('1'):
                cbar_labels.append(label)
            else:
                cbar_labels.append('')

    cbar.ax.set_yticklabels(cbar_labels)
    if log_plot:
        cbar.ax.tick_params(labelsize=4, width=1, length=3, pad=0.75)

    return


def add_bottom_colorbar(ax, colorbar, plot_df, pad_dist, xlabel_fontsize):
    """

    :param ax:
    :param colorbar:
    :param plot_df:
    :param pad_dist:
    :param xlabel_fontsize:
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
    if pad_dist:
        axb.tick_params(axis='x', pad=pad_dist)
    plt.setp(
        axb.xaxis.get_majorticklabels(), rotation=90, fontsize=xlabel_fontsize
    )
    return


def masked_double_heatmap(df_to_plot, cols_to_mask, fig_file, colorbar={},
                          masked_cmap=cm.Blues, other_cmap=cm.Greens,
                          title='', vline_pos=[], log_plot=True,
                          other_cbar_label='', size=None, bottom_pad=0,
                          xlabel_fontsize=5, cbar_fraction=0.046,
                          cbar_pad=0.032, cbar_labels=[], cbar_font_adder=1.5):
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
    :param other_cbar_label:
    :param size:
    :param bottom_pad:
    :param xlabel_fontsize:
    :param cbar_fraction:
    :param cbar_pad:
    :param cbar_labels:
    :param cbar_font_adder:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    if size:
        plt.rcParams['figure.figsize'] = size
        cbar_fontsize = min(size) + cbar_font_adder
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
    ax.pcolormesh(white_df, cmap='binary')

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
    heatmap_colorbar(
        masked_mesh, bar_label, log_plot, cbar_fontsize, cbar_fraction,
        cbar_pad, cbar_labels
    )

    if cols_to_mask != df_to_plot.columns.values.tolist():
        bar_label = '                      ' + other_cbar_label
        heatmap_colorbar(
            other_mesh, bar_label, log_plot, cbar_fontsize, cbar_fraction,
            cbar_pad, cbar_labels
        )

    ax.set_xticks(np.arange(0, df_to_plot.shape[1], 1), minor=True)
    ax.grid(axis='x', ls='-', color='white', which='minor')

    if colorbar:
        add_bottom_colorbar(
            ax, colorbar, df_to_plot, bottom_pad, xlabel_fontsize
        )
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
                                right_lim_shift=2, back_cols={}, font_cols={},
                                tab_fontsize=5.5, intab_fontsize=0):
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
        if back_cols:
            col_cols.append(back_cols[abbr])
        else:
            try:
                col_cols.append('xkcd:{}'.format(_CANCER_COLORS[abbr]))
            except KeyError:
                for can_abbr in _CANCER_COLORS.keys():
                    if can_abbr in abbr:
                        col_cols.append(
                            'xkcd:{}'.format(_CANCER_COLORS[can_abbr])
                        )
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
        if font_cols:
            font_color = font_cols[abbr]
        else:
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
    if intab_fontsize:
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(intab_fontsize)

    for (row, col), cell in the_table.get_celld().items():
        if (row == 0):
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=tab_fontsize)
            )
            if len(columns[0]) > 5:
                cell.set_height(cell.get_height() * 1.5)
            if col in whitefont_cols:
                cell._text.set_color('white')
        if col == -1:
            cell.set_text_props(
                fontproperties=FontProperties(weight='bold', size=tab_fontsize)
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
