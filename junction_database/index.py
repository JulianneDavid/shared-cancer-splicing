import csv
import gc
import os
import pandas as pd
import sqlite3 as sql
import sys
import time
try:
    from utilities.utilities import _PHEN_TABLE, _JX_ANN_TABLE, _JX_SAMP_TABLE
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _PHEN_TABLE, _JX_ANN_TABLE, _JX_SAMP_TABLE
from utilities.utilities import _CONSTANT_CANCER_TYPES
from utilities.utilities import gtf_to_cds, cds_to_tree, extract_splice_sites
from utilities.utilities import check_annotations, jx_ends_in_cds
from utilities.utilities import make_id_name_dict

class JXIndexError(Exception):
    pass

def primary_cancer(main_tcga_type, brain_cervical_histology,
                   esca_sarc_pcpg_meso_uter_histology):
    subtype = ''
    cesc = 'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma'
    esca = 'Esophageal_Carcinoma'
    lgg = 'Brain_Lower_Grade_Glioma'
    sarc = 'Sarcoma'
    pcpg = 'Pheochromocytoma_and_Paraganglioma'

    if main_tcga_type == lgg:
        types = ['Oligoastrocytoma', 'Astrocytoma', 'Oligodendroglioma']
        if brain_cervical_histology in types:
            subtype = brain_cervical_histology
    elif main_tcga_type == cesc:
        adenosquamous = 'Adenosquamous'
        squamous = 'Cervical Squamous Cell Carcinoma'
        adenocarc = [
            'Endometrioid Adenocarcinoma of Endocervix',
            'Endocervical Adenocarcinoma of the Usual Type',
            'Mucinous Adenocarcinoma of Endocervical Type',
            'Endocervical Type of Adenocarcinoma'
        ]
        if brain_cervical_histology == adenosquamous:
            subtype = 'Cervical_Adenosquamous'
        elif brain_cervical_histology in adenocarc:
            subtype = 'Endocervical_Adenocarcinoma'
        elif brain_cervical_histology == squamous:
            subtype = brain_cervical_histology.replace(' ', '_')
        else:
            print(brain_cervical_histology)
    elif main_tcga_type == sarc:
        lms = 'Leiomyosarcoma (LMS)'
        myxo = 'Myxofibrosarcoma'
        mpnst = 'Malignant Peripheral Nerve Sheath Tumors (MPNST)'
        des = 'Desmoid Tumor'
        dediff = 'Dedifferentiated liposarcoma'
        syn = [
            'Sarcoma; synovial; poorly differentiated',
            'Synovial Sarcoma - Biphasic', 'Synovial Sarcoma - Monophasic'
        ]
        ups = [
            'Undifferentiated Pleomorphic Sarcoma (UPS)',
            "Pleomorphic 'MFH' / Undifferentiated pleomorphic sarcoma",
            "Giant cell 'MFH' / Undifferentiated pleomorphic sarcoma with "
            "giant cells"
        ]
        if esca_sarc_pcpg_meso_uter_histology == lms:
            subtype = 'Leiomyosarcoma'
        elif esca_sarc_pcpg_meso_uter_histology == myxo:
            subtype = myxo
        elif esca_sarc_pcpg_meso_uter_histology == des:
            subtype = des.replace(' ', '_')
        elif esca_sarc_pcpg_meso_uter_histology == mpnst:
            subtype = 'Malignant_Peripheral_Nerve_Sheath_Tumors'
        elif esca_sarc_pcpg_meso_uter_histology == dediff:
            subtype = dediff.replace(' ', '_')
        elif esca_sarc_pcpg_meso_uter_histology in syn:
            subtype = 'Synovial_Sarcoma'
        elif esca_sarc_pcpg_meso_uter_histology in ups:
            subtype = 'Undifferentiated_Pleomorphic_Sarcoma'
    elif main_tcga_type == esca:
        squame = 'Esophagus Squamous Cell Carcinoma'
        adeno = 'Esophagus Adenocarcinoma, NOS'
        if esca_sarc_pcpg_meso_uter_histology == adeno:
            subtype = adeno.split(', NOS')[0].replace(' ', '_')
        elif esca_sarc_pcpg_meso_uter_histology == squame:
            subtype = squame.replace(' ', '_')
    elif main_tcga_type == pcpg:
        pheo = 'Pheochromocytoma'
        para = [
            'Paraganglioma; Extra-adrenal Pheochromocytoma', 'Paraganglioma'
        ]
        if esca_sarc_pcpg_meso_uter_histology in para:
            subtype = 'Paraganglioma'
        elif esca_sarc_pcpg_meso_uter_histology == pheo:
            subtype = 'Pheochromocytoma'
    elif main_tcga_type in _CONSTANT_CANCER_TYPES:
        subtype = main_tcga_type
    return subtype


def collect_stage(raw_stage_data):
    stage = ''
    try:
        stage_check = raw_stage_data.split('stage ')[-1]
        if stage_check[-1] in ['i', 'a', 'b', 'c', '0', 'v']:
            stage = stage_check
        return stage
    except AttributeError:
        return stage


def seminoma_info(cancer_type, seminoma_label):
    label = ''
    if cancer_type == 'Testicular_Germ_Cell_Tumors':
        if type(seminoma_label) != str:
            return label
        if seminoma_label.startswith('Seminoma'):
            label = 0
        elif seminoma_label != 'NA':
            label = 1
    return label


def gleason_score(cancer_type, prostate_info):
    gleason = ''
    if cancer_type == 'Prostate_Adenocarcinoma':
        if type(prostate_info) == str:
            gleason = prostate_info[0]
    return gleason


def mesothelioma_subtype(cancer_type, mesothelioma_info):
    meso_subtype = ''
    if cancer_type == 'Mesothelioma' and type(mesothelioma_info) == str:
        meso_subtype = mesothelioma_info.split(' - NOS')[0]
    return meso_subtype


def collect_hpv_status(cancer_type, hpv_info_1, hpv_info_2):
    hpv_status = ''
    if cancer_type == 'Head_and_Neck_Squamous_Cell_Carcinoma':
        if hpv_info_1 == 'Positive' or hpv_info_2 == 'Positive':
            hpv_status = 0
        elif hpv_info_1 == 'Negative' or hpv_info_2 == 'Negative':
            hpv_status = 1
    return hpv_status


def collect_hnsc_loc(cancer_type, location):
    hnsc_loc = ''
    locations = {
        'Alveolar Ridge': 'Oral_Cavity', 'Oral Tongue': 'Oral_Cavity',
        'Hypopharynx': 'Hypopharynx', 'Tonsil': 'Oropharynx',
        'Larynx': 'Larynx', 'Base of tongue': 'Oropharynx',
        'Oral Cavity': 'Oral_Cavity', 'Floor of mouth': 'Oral_Cavity',
        'Buccal Mucosa': 'Oral_Cavity', 'Lip': 'Oral_Cavity',
        'Oropharynx': 'Oropharynx', 'Hard Palate': 'Oral_Cavity'
    }
    if cancer_type == 'Head_and_Neck_Squamous_Cell_Carcinoma':
        if type(location) == str:
            hnsc_loc = locations[location]
    return hnsc_loc


def hepatitis_status(cancer_type, hepatitis_label):
    hep_status = ''
    if cancer_type == 'Liver_Hepatocellular_Carcinoma':
        if type(hepatitis_label) != str:
            return hep_status
        if 'Hepatitis' in hepatitis_label:
            hep_status = 0
        elif hepatitis_label != 'NA':
            hep_status = 1
    return hep_status


def collect_ucs_subtype(cancer_type, ucs_info):
    ucs_subtype = ''
    if cancer_type == 'Uterine_Carcinosarcoma':
        if ucs_info == 'Uterine Carcinosarcoma/ MMMT: Heterologous Type':
            ucs_subtype = 0
        elif ucs_info == 'Uterine Carcinosarcoma/MMMT: Homologous Type':
            ucs_subtype = 1
    return ucs_subtype


def tumor_normal(tn_flag):
    tumor_stat = ''
    if type(tn_flag) == str:
        if tn_flag == 'Solid Tissue Normal':
            tumor_stat = 1
        else:
            tumor_stat = 0
    return tumor_stat


def vit_stat_map(vital_status):
    status = ''
    if type(vital_status) == str:
        if vital_status == 'alive':
            status = 1
        elif vital_status == 'dead':
            status = 0
    return status


def gender_map(gender):
    label = ''
    if gender == 'male':
        label = 0
    elif gender == 'female':
        label = 1
    return label


def create_phenotype_table(cancer_phens, tissue_phens, id_file, index_db,
                           conn):
    """Parses TCGA phenotype file and stores info in new database table

    Input:
    cancer_phenotypes: .tsv file of TCGA phenotype information from recount2.
    id_file: file mapping recount2 ids to TCGA project ids.
    _PHEN_IND: the name to assign to the phenotype table within the database.
    index_db: the SQL cursor for executing sqlite3 commands.

    Loads desired columns from the phenotype .tsv into a pandas dataframe. For
    each cancer sample, collects and cleans column contents. Creates a new row
    entry SQL command and executes it. Commits the new table in SQL.

    Keys:
    gender: 0 = male, 1 = female
    hep_status: 0 = positive, 1 = negative
    hpv_status: 0 = positive, 1 = negative
    project: 0 = TCGA, 1 = GTEx
    seminoma_status: 0 = seminoma, 1 = non seminoma
    tumor_normal: 0 = tumor, 1 = normal
    ucs_subtype: 0 = heterologous, 1 = homologous
    vital_status: 0 = dead, 1 = alive

    Returns None
    """
    index_db.execute('DROP TABLE IF EXISTS {};'.format(_PHEN_TABLE))
    phen_command = (
        'CREATE TABLE {}(tcga_id VARCHAR, case_id VARCHAR, '
        'gender UNSIGNED TINYINT, gleason UNSIGNED TINYINT, '
        'hep_status UNSIGNED TINYINT, hnsc_loc VARCHAR, '
        'hpv_status UNSIGNED TINYINT, meso_subtype VARCHAR, '
        'primary_type VARCHAR, project UNSIGNED TINYINT, '
        'project_type_label VARCHAR, recount_id UNSIGNED INT, '
        'sample_origin VARCHAR, seminoma_status UNSIGNED TINYINT, '
        'stage VARCHAR, tumor_normal UNSIGNED TINYINT, '
        'ucs_subtype UNSIGNED TINYINT, universal_id VARCHAR, '
        'vital_status UNSIGNED TINYINT);'.format(_PHEN_TABLE)
    )
    index_db.execute(phen_command)
    conn.commit()

    sample_ids = pd.read_table(
        id_file, header=None, usecols=[0, 2],
        names=['recount_id', 'universal_id']
    )
    sample_ids = sample_ids.set_index('universal_id')['recount_id']

    uppercase_and_spaces = {
        'gdc_file_id':
            lambda ident: ident.upper(),
        'gdc_cases.project.name':
            lambda name: name.replace('-', '_').replace(' ', '_')
    }
    tcga_phen = pd.read_table(
        cancer_phens,
        usecols=['gdc_file_id', 'gdc_cases.case_id', 'gdc_cases.submitter_id',
                 'project', 'gdc_cases.project.name',
                 'gdc_cases.diagnoses.tumor_stage',
                 'gdc_cases.diagnoses.vital_status',
                 'gdc_cases.demographic.gender',
                 'gdc_cases.samples.sample_type',
                 'cgc_case_histological_diagnosis',
                 'xml_primary_pathology_histological_type',
                 'xml_primary_pathology_histology_list',
                 'xml_stage_event_gleason_grading',
                 'xml_hpv_status_by_p16_testing',
                 'xml_hpv_status_by_ish_testing',
                 'xml_anatomic_neoplasm_subdivision',
                 'xml_history_hepato_carcinoma_risk_factors'
                 ],
        converters=uppercase_and_spaces, dtype=str
    )
    new_names = {
        'gdc_file_id': 'universal_id', 'gdc_cases.submitter_id': 'TCGA_id',
        'gdc_cases.case_id': 'case_id', 'project': 'project',
        'gdc_cases.project.name': 'project_type_label',
        'gdc_cases.diagnoses.tumor_stage': 'raw_stage',
        'gdc_cases.diagnoses.vital_status': 'vital_status',
        'gdc_cases.demographic.gender': 'gender',
        'cgc_case_histological_diagnosis': 'brain_cervical_info',
        'xml_primary_pathology_histological_type':
            'esca_sarc_pcpg_meso_uter_info',
        'xml_primary_pathology_histology_list': 'testicular_info',
        'xml_stage_event_gleason_grading': 'prostate_info',
        'xml_hpv_status_by_p16_testing': 'hnsc_hpv1',
        'xml_hpv_status_by_ish_testing': 'hnsc_hpv2',
        'xml_anatomic_neoplasm_subdivision': 'hnsc_loc',
        'xml_history_hepato_carcinoma_risk_factors': 'hepatitis_info',
        'gdc_cases.samples.sample_type': 'sample_origin'
    }
    tcga_phen.rename(new_names, axis='columns', inplace=True)
    tcga_phen = tcga_phen[tcga_phen.project == 'TCGA']

    tcga_phen['primary_type'] = tcga_phen.apply(
        lambda x: primary_cancer(
            x['project_type_label'], x['brain_cervical_info'],
            x['esca_sarc_pcpg_meso_uter_info']
        ), axis=1
    )
    tcga_phen['seminoma'] = tcga_phen.apply(
        lambda x: seminoma_info(
            x['project_type_label'], x['testicular_info']
        ), axis=1
    )
    tcga_phen['gleason'] = tcga_phen.apply(
        lambda x: gleason_score(
            x['project_type_label'], x['prostate_info']
        ), axis=1
    )
    tcga_phen['meso_subtype'] = tcga_phen.apply(
        lambda x: mesothelioma_subtype(
            x['project_type_label'], x['esca_sarc_pcpg_meso_uter_info']
        ), axis=1
    )
    tcga_phen['hpv_status'] = tcga_phen.apply(
        lambda x: collect_hpv_status(
            x['project_type_label'], x['hnsc_hpv1'], x['hnsc_hpv2']
        ), axis=1
    )
    tcga_phen['hnsc_location'] = tcga_phen.apply(
        lambda x: collect_hnsc_loc(
            x['project_type_label'], x['hnsc_loc']
        ), axis=1
    )
    tcga_phen['hepatitis_status'] = tcga_phen.apply(
        lambda x: hepatitis_status(
            x['project_type_label'], x['hepatitis_info']
        ), axis=1
    )
    tcga_phen['ucs_subtype'] = tcga_phen.apply(
        lambda x: collect_ucs_subtype(
            x['project_type_label'], x['esca_sarc_pcpg_meso_uter_info']
        ), axis=1
    )
    tcga_phen['recount_id'] = tcga_phen.universal_id.apply(
        lambda x: sample_ids[x]
    )
    tcga_phen['stage'] = tcga_phen.raw_stage.apply(
        lambda x: collect_stage(x)
    )
    tcga_phen['project'] = 0
    tcga_phen['tumor_normal'] = tcga_phen.sample_origin.apply(
        lambda x: tumor_normal(x)
    )
    tcga_phen.drop(
        ['raw_stage', 'brain_cervical_info', 'esca_sarc_pcpg_meso_uter_info',
         'testicular_info', 'prostate_info', 'hnsc_hpv1', 'hnsc_hpv2',
         'hnsc_loc', 'hepatitis_info'], inplace=True, axis=1
    )

    gtex_phen = pd.read_table(
        tissue_phens, usecols=['sampid', 'smts', 'run', 'smtsd']
    )
    new_names = {
        'sampid': 'gtex_id', 'smts': 'project_type_label',
        'run': 'universal_id', 'smtsd': 'sample_origin'
    }
    gtex_phen.rename(new_names, axis='columns', inplace=True)
    gtex_phen = gtex_phen[gtex_phen.gtex_id.str.startswith('GTEX')]
    gtex_phen.dropna(subset=['project_type_label'], inplace=True)
    gtex_phen['stage'] = ''
    gtex_phen['primary_type'] = gtex_phen['project_type_label']
    gtex_phen['tumor_normal'] = 1
    gtex_phen['project'] = 1
    gtex_phen['recount_id'] = gtex_phen.universal_id.apply(
        lambda x: sample_ids[x]
    )

    gtex_phen.drop(['gtex_id'], inplace=True, axis=1)

    full_df = pd.concat([tcga_phen, gtex_phen], ignore_index=True).fillna('')

    for ind, row in full_df.iterrows():
        # jx_info = list(map(str, jx_info + cds_cols))

        items = "'" + "', '".join(list(map(str, row.tolist()))) + "'"
        row_command = (
            'INSERT INTO {} VALUES({});'.format(_PHEN_TABLE, items)
        )
        try:
            index_db.execute(row_command)
        except sql.IntegrityError:
            continue

    phen_index1 = (
        'CREATE INDEX primary_type_index ON {}(primary_type);'
        ''.format(_PHEN_TABLE)
    )
    index_db.execute(phen_index1)

    phen_index2 = (
        'CREATE INDEX project_index ON {}(project);'.format(_PHEN_TABLE)
    )
    index_db.execute(phen_index2)

    phen_index3 = (
        'CREATE INDEX recount_id_index ON {}(recount_id);'.format(_PHEN_TABLE)
    )
    index_db.execute(phen_index3)

    phen_index4 = (
        'CREATE INDEX project_label_index ON {}(project_type_label);'
        ''.format(_PHEN_TABLE)
    )
    index_db.execute(phen_index4)

    conn.commit()
    print('phenotype table creating complete, moving to indexing')
    return


def parse_phenotypes(tissue_phen, cancer_phen, id_file):
    """Parses phenotype and sample id files to generate a sample id translator.

    Input:
        tissue_phen: Recount-generated file containing sample SRA numbers
            and tissue types for normal tissue samples.
        cancer_phen: Recount-generated file containing experiment id numbers
            and cancer types for cancer samples.
        id_file: Recount-generated file for translating Recount sample IDs to
            phenotype experiment ids (SRA or TCGA).

    Returns sorted lists of the unique tissue and cancer types present in the
    phenotype files, and a id_decoder data series for translating
    Recount IDs into the appropriate tissue and cancer types.

    """
    gtex_phen = pd.read_table(
        tissue_phen, usecols=['sampid', 'smts', 'run', 'auc']
    )
    new_names = {
        'sampid': 'gtex_id', 'smts': 'sample_type', 'run': 'universal_id'
    }
    gtex_phen.rename(new_names, axis='columns', inplace=True)
    gtex_phen = gtex_phen[gtex_phen.gtex_id.str.startswith('GTEX')]
    gtex_phen.dropna(subset=['sample_type'], inplace=True)
    gtex_phen['tumor_stage'] = 'NA'
    tissue_id_index = gtex_phen.set_index('universal_id')['sample_type']
    tissue_auc_index = gtex_phen.set_index('universal_id')['auc']
    tissue_stage_index = gtex_phen.set_index('universal_id')['tumor_stage']
    tissue_counts = tissue_id_index.value_counts().to_dict()
    uppercase_fixer = {'gdc_file_id': lambda ident: ident.upper()}
    tcga_phen = pd.read_table(
        cancer_phen,
        usecols=['gdc_file_id', 'project', 'gdc_cases.project.name', 'auc',
                 'gdc_cases.diagnoses.tumor_stage'],
        converters=uppercase_fixer,
        dtype={
            'project': str, 'auc': str, 'gdc_cases.project.name': str,
            'gdc_cases.diagnoses.tumor_stage': str
        }
    )
    new_names = {
        'gdc_file_id': 'universal_id', 'project': 'tcga',
        'gdc_cases.project.name': 'sample_type',
        'gdc_cases.diagnoses.tumor_stage': 'tumor_stage'
    }
    tcga_phen.rename(new_names, axis='columns', inplace=True)
    tcga_phen = tcga_phen[tcga_phen.tcga == 'TCGA']
    tcga_phen['auc'] = tcga_phen['auc'].astype(int)
    cancer_id_index = tcga_phen.set_index('universal_id')['sample_type']
    cancer_auc_index = tcga_phen.set_index('universal_id')['auc']
    cancer_stage_index = tcga_phen.set_index('universal_id')['tumor_stage']
    cancer_counts = cancer_id_index.value_counts().to_dict()

    id_index = pd.concat([tissue_id_index, cancer_id_index])
    auc_index = pd.concat([tissue_auc_index, cancer_auc_index])
    stage_index = pd.concat([tissue_stage_index, cancer_stage_index])
    sample_ids = pd.read_table(
        id_file, header=None, usecols=[0, 2],
        names=['recount_id', 'universal_id']
    )
    sample_ids = sample_ids.set_index('recount_id')['universal_id']
    id_decoder = sample_ids.map(id_index).dropna().to_dict()
    auc_decoder = sample_ids.map(auc_index).dropna().to_dict()
    stage_decoder = sample_ids.map(stage_index).dropna().to_dict()

    del sample_ids
    del tissue_id_index
    del tissue_auc_index
    del cancer_id_index
    del cancer_auc_index
    del id_index
    del gtex_phen
    del tcga_phen
    gc.collect()
    return tissue_counts, cancer_counts, id_decoder, auc_decoder, stage_decoder


def add_jxs_to_db(cov_file, bed_file, auc_decoder, index_db, cds_interval_tree,
                  id_name_dict, jx_annotations, conn, start_time, jxs_added):
    """

    :param cov_file:
    :param bed_file:
    :param auc_decoder:
    :param index_db:
    :param cds_interval_tree:
    :param id_name_dict:
    :param jx_annotations:
    :param conn:
    :param start_time:
    :param jxs_added:
    :return:
    """
    fill_int = time.time()
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

            jx = ';'.join([chrom, left, right, strand])
            ids = ids.split(',')
            covs = covs.split(',')
            for id, cov in zip(ids, covs):
                try:
                    auc = auc_decoder[int(id)]
                except KeyError:
                    continue
                cov = str(1e10 * float(cov) / auc)
                items = "'" + "', '".join([jx_id, id, cov]) + "'"
                row_command = (
                    'INSERT INTO {} VALUES ({});'.format(_JX_SAMP_TABLE, items)
                )
                try:
                    index_db.execute(row_command)
                except sql.IntegrityError:
                    continue

            if jx not in jxs_added:
                cds_cols = jx_ends_in_cds(jx, cds_interval_tree, id_name_dict)
                annotation_col = check_annotations(jx, jx_annotations)
                cds_cols = list(cds_cols) + [annotation_col]
                jx_info = [jx_id, jx, chrom, left, right, strand]
                jx_info = list(map(str, jx_info + cds_cols))
                jxs_added.add(jx)

                values = "'" + "', '".join(jx_info) + "'"
                row_command = (
                    'INSERT INTO {} VALUES ({});'.format(_JX_ANN_TABLE, values)
                )
                try:
                    index_db.execute(row_command)
                except sql.IntegrityError:
                    continue

            if (i % 1000000) == 0:
                fill_2000 = time.time()
                print('\n{}th entry, writing'.format(i))
                print('intermediate fill time is', fill_2000 - fill_int)
                print('total fill time is', fill_2000 - start_time)
                conn.commit()
                fill_int = time.time()
    return jxs_added


def fill_new_table(tcga_cov, tcga_bed, gtex_cov, gtex_bed, index_db,
                   cds_interval_tree, jx_annotations, id_name_dict,
                   auc_decoder, conn):
    """Creates a new SQL database to store ID and coverage data for junctions.

    Input:
        excluded: the sorted list of unique excluded generated by
            parse_phenotype
        cancers: the sorted list of unique cancers generated by parse_phenotype
        tissue_ids_table: the name to assign to the tissue junction index table
            within the database.
        cancer_ids_table: the name to assign to the cancer juntion index table
            within the database.
        index_db: the SQL cursor for executing sqlite3 commands.

    Creates column names for each tissue and cancer type, for both sample IDs
    that contain a given junction and their associated coverages.  Then creates
    and executes a SQL command for generating tissue and cancer junction tables
    within the database.

    Returns None
    """
    fill_start = time.time()

    # index_db.execute('DROP TABLE IF EXISTS {};'.format(_JX_SAMP_TABLE))
    jx_samp_command = (
        'CREATE TABLE {}(jx_id UNSIGNED INT, recount_id UNSIGNED INT, '
        'coverage DOUBLE);'.format(_JX_SAMP_TABLE)
    )
    index_db.execute(jx_samp_command)

    # index_db.execute('DROP TABLE IF EXISTS {};'.format(_JX_ANN_TABLE))
    jx_ann_command = (
        'CREATE TABLE {}(jx_id UNSIGNED INT, jx VARCHAR, chrom VARCHAR, '
        'left UNSIGNED INT, right UNSIGNED INT, strand VARCHAR, '
        'no_cds UNSIGNED TINYINT, one_cds UNSIGNED TINYINT, '
        'two_cds_diff UNSIGNED TINYINT, two_cds_same UNSIGNED TINYINT, '
        'fiveprime_genes VARCHAR, threeprime_genes VARCHAR, '
        'fiveprime_aliases VARCHAR, threeprime_aliases VARCHAR, '
        'annotation UNSIGNED TINYINT);'.format(_JX_ANN_TABLE)
    )
    index_db.execute(jx_ann_command)

    previously_added_jxs = set()
    previously_added_jxs = add_jxs_to_db(
        tcga_cov, tcga_bed, auc_decoder, index_db, cds_interval_tree,
        id_name_dict, jx_annotations, conn, fill_start, previously_added_jxs
    )
    add_jxs_to_db(
        gtex_cov, gtex_bed, auc_decoder, index_db, cds_interval_tree,
        id_name_dict, jx_annotations, conn, fill_start, previously_added_jxs
    )
    conn.commit()
    print('all junctions added to db!  adding db indexes.')
    fill_int = time.time()
    jx_ann_index = (
        'CREATE INDEX jx_ann_id_index ON {}(jx_id);'.format(_JX_ANN_TABLE)
    )
    try:
        index_db.execute(jx_ann_index)
    except sql.OperationalError:
        print('jx_ann_index operational error')
        pass
    conn.commit()
    fill_2000 = time.time()
    print('\nfirst index done')
    print('intermediate time is', fill_2000 - fill_int)
    print('total time is', fill_2000 - fill_start)
    fill_int = time.time()

    jx_ann_index2 = (
        'CREATE INDEX jx_index ON {}(jx);'
        ''.format(_JX_ANN_TABLE)
    )
    try:
        index_db.execute(jx_ann_index2)
    except sql.OperationalError:
        print('jx ann_index2 operational error')
        pass
    conn.commit()
    fill_2000 = time.time()
    print('\nsecond index done')
    print('intermediate time is', fill_2000 - fill_int)
    print('total time is', fill_2000 - fill_start)
    fill_int = time.time()

    jx_samp_index = (
        'CREATE INDEX samp_id_index ON {}(recount_id);'
        ''.format(_JX_SAMP_TABLE)
    )
    try:
        index_db.execute(jx_samp_index)
    except sql.OperationalError:
        print('jxsampleindex operational error')
        pass
    conn.commit()
    fill_2000 = time.time()
    print('\nthird index done')
    print('intermediate time is', fill_2000 - fill_int)
    print('total time is', fill_2000 - fill_start)
    fill_int = time.time()

    jx_samp_index2 = (
        'CREATE INDEX jx_samp_id_index ON {}(jx_id); '.format(_JX_SAMP_TABLE)
    )
    try:
        index_db.execute(jx_samp_index2)
    except sql.OperationalError:
        print('jxsamp2 operational error')
        pass
    conn.commit()
    fill_2000 = time.time()
    print('\nfourth index done')
    print('intermediate time is', fill_2000 - fill_int)
    print('FINAL total time is', fill_2000 - fill_start)
    return


def main(args, conn, index_db):
    tcga_c = args.tcga_jx_cov
    tcga_b = args.tcga_jx_bed
    tcga_p = args.tcga_phenotype
    gtex_c = args.gtex_jx_cov
    gtex_b = args.gtex_jx_bed
    gtex_p = args.gtex_phenotype
    id_file = args.sample_id_file
    gtf_path = args.gtf_file

    results = parse_phenotypes(gtex_p, tcga_p, id_file)
    tissues, cancers, id_decoder, auc_decoder, stage_decoder = results
    create_phenotype_table(tcga_p, gtex_p, id_file, index_db, conn)
    coding_regions = gtf_to_cds(gtf_path)
    print('coding regions discovered')
    CDS_interval_tree = cds_to_tree(coding_regions)
    print('CDS tree created')
    jx_annotations = extract_splice_sites(gtf_path)
    print('splice sites extracted')
    id_name_dict = make_id_name_dict(gtf_path)

    fill_new_table(
        tcga_c, tcga_b, gtex_c, gtex_b, index_db, CDS_interval_tree,
        jx_annotations, id_name_dict, auc_decoder, conn
    )
