import junction_database.index as index
import numpy as np
import unittest

class TestJxtapose(unittest.TestCase):
    def setUp(self):
        pass

    def test_is_not_tumor(self):
        normal_status = {
            'Solid Tissue Normal': True, 'Metastatic': False,
            'Additional - New Primary': False, 'Additional Metastatic': False,
            'Primary Tumor': False, 'Recurrent Tumor': False,
            'Primary Blood Derived Cancer - Peripheral Blood': False,
            np.nan: ''
        }
        for tissue, status in normal_status.items():
            self.assertEqual(index.tumor_normal(tissue), status)

    def test_collect_stage(self):
        raw_stage_vals = {
            'i/ii nos': '', 'is': '', 'NA': '', 'not reported': '',
            'stage 0': '0', 'stage i': 'i', 'stage ia': 'ia', 'stage ib': 'ib',
            'stage ii': 'ii', 'stage iia': 'iia', 'stage iib': 'iib',
            'stage iic': 'iic', 'stage iii': 'iii', 'stage iiia': 'iiia',
            'stage iiib': 'iiib', 'stage iiic': 'iiic', 'stage iv': 'iv',
            'stage iva': 'iva', 'stage ivb': 'ivb', 'stage ivc': 'ivc',
            'stage x': '', 'other': '', np.nan: ''
        }
        for raw_stage, final_stage in raw_stage_vals.items():
            self.assertEqual(index.collect_stage(raw_stage), final_stage)

    def test_seminoma_info(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        seminoma_options = {
            'Non-Seminoma; Embryonal Carcinoma100': 1, 'Seminoma; NOS100': 0,
            'Seminoma; NOS80Non-Seminoma; Choriocarcinoma20': 0,
            'Non-Seminoma; Yolk Sac Tumor100': 1, 'NA': '', 0.0: '', np.nan: ''
        }
        for label, answer in seminoma_options.items():
            self.assertEqual(
                index.seminoma_info('Testicular_Germ_Cell_Tumors', label),
                answer
            )
        for cancer in cancers:
            for label in seminoma_options:
                self.assertEqual(index.seminoma_info(cancer, label), '')

    def test_gleason_score(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        prostate_info_options = {
            '633': '6', '734': '7', '853': '8', '7435': '7', '9543': '9',
            np.nan: ''
        }
        for label, answer in prostate_info_options.items():
            self.assertEqual(
                index.gleason_score(
                    'Prostate_Adenocarcinoma', label
                ), answer
            )
        for cancer in cancers:
            for label in prostate_info_options:
                self.assertEqual(
                    index.seminoma_info(cancer, label), ''
                )

    def test_mesothelioma_subtype(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        subtype_options = {
            'Diffuse malignant mesothelioma - NOS':
                'Diffuse malignant mesothelioma',
            'Epithelioid mesothelioma': 'Epithelioid mesothelioma',
            'Biphasic mesothelioma': 'Biphasic mesothelioma',
            'Sarcomatoid mesothelioma': 'Sarcomatoid mesothelioma',
            np.nan: ''
        }
        for label, answer in subtype_options.items():
            self.assertEqual(
                index.mesothelioma_subtype('Mesothelioma', label), answer
            )
        for cancer in cancers:
            for label in subtype_options:
                self.assertEqual(
                    index.mesothelioma_subtype(cancer, label), ''
                )

    def test_hpv_status(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        hpv_sets = {
            ('NA', 'NA'): '', ('Positive', 'Positive'): 0,
            ('NA', 'Positive'): 0, ('Positive', 'NA'): 0,
            ('Negative', 'Positive'): 0, ('Positive', 'Negative'): 0,
            ('NA', 'Negative'): 1, ('Negative', 'NA'): 1,
            ('Negative', 'Negative'): 1,
            ('NA', np.nan): '', (np.nan, 'NA'): '', (np.nan, np.nan): ''
        }
        for (l1, l2), answer in hpv_sets.items():
            self.assertEqual(
                index.collect_hpv_status(
                    'Head_and_Neck_Squamous_Cell_Carcinoma', l1, l2
                ), answer
            )
        for cancer in cancers:
            for (l1, l2) in hpv_sets:
                    self.assertEqual(
                        index.collect_hpv_status(cancer, l1, l2), ''
                    )

    def test_hnsc_locs(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        locations = {
            'Alveolar Ridge': 'Oral_Cavity', 'Oral Tongue': 'Oral_Cavity',
            'Hypopharynx': 'Hypopharynx', 'Tonsil': 'Oropharynx',
            'Larynx': 'Larynx', 'Base of tongue': 'Oropharynx',
            'Oral Cavity': 'Oral_Cavity', 'Floor of mouth': 'Oral_Cavity',
            'Buccal Mucosa': 'Oral_Cavity', 'Lip': 'Oral_Cavity',
            'Oropharynx': 'Oropharynx', 'Hard Palate': 'Oral_Cavity',
            np.nan: ''
        }
        for label, answer in locations.items():
            self.assertEqual(
                index.collect_hnsc_loc(
                    'Head_and_Neck_Squamous_Cell_Carcinoma', label
                ), answer
            )
        for cancer in cancers:
            for loc in locations:
                self.assertEqual(
                    index.collect_hnsc_loc(cancer, loc), ''
                )

    def test_hepatitis_status(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        hep_options = {
            'Hepatitis B': 0, 'Alcohol consumptionHepatitis C': 0,
            'Hepatitis C': 0, 'Alcohol consumptionHepatitis B': 0,
            'Hepatitis CCirrhosis, Alcohol abuse': 0,
            'Hepatitis BHepatitis C': 0, 'Hepatitis Cno': 0,
            'Alcohol consumptionHepatitis BHepatitis C': 0,
            'Alcohol consumptionHepatitis BHepatitis CHemochromatosis': 0,
            'No History of Primary Risk Factors': 1, 'Hemochromatosis': 1,
            'OtherHBcAB total (+), HBsAg (-), HBcIgm Ab (-)': 1,
            'OtherWorker at steel factory': 1,
            'NA': '', np.nan: '', 0.0: ''
        }
        for label, answer in hep_options.items():
            self.assertEqual(
                index.hepatitis_status(
                    'Liver_Hepatocellular_Carcinoma', label
                ), answer
            )
        for cancer in cancers:
            for label in hep_options.items():
                self.assertEqual(
                    index.hepatitis_status(cancer, label), ''
                )

    def test_ucs_subtype(self):
        cancers = ['Skin_Cutaneous_Melanoma', 'NA', np.nan, 0.0]
        subtypes = {
            'Uterine Carcinosarcoma/ '
            'Malignant Mixed Mullerian Tumor (MMMT): NOS':
                '',
            'Uterine Carcinosarcoma/ MMMT: Heterologous Type': 0,
            'Uterine Carcinosarcoma/MMMT: Homologous Type': 1,
            'other': '', np.nan: ''
        }
        for label, answer in subtypes.items():
            self.assertEqual(
                index.collect_ucs_subtype(
                    'Uterine_Carcinosarcoma', label
                ), answer
            )
        for cancer in cancers:
            for subt in subtypes:
                self.assertEqual(
                    index.collect_ucs_subtype(cancer, subt), ''
                )

    def test_vital_status(self):
        options = {
            'alive': 1, 'dead': 0, 'not reported': '', 'NA': '', 0: '',
            np.nan: ''
        }
        for option, answer in options.items():
            self.assertEqual(index.vit_stat_map(option), answer)

    def test_gender_map(self):
        options = {
            'female': 1, 'male': 0, 'NA': '', 'NX': '', np.nan: ''
        }
        for option, answer in options.items():
            self.assertEqual(index.gender_map(option), answer)

    def test_primary_type(self):
        cesc = (
            'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma'
        )
        esca = 'Esophageal_Carcinoma'
        lgg = 'Brain_Lower_Grade_Glioma'
        sarc = 'Sarcoma'
        pcpg = 'Pheochromocytoma_and_Paraganglioma'

        tcga_options = {
            'Skin_Cutaneous_Melanoma': 'Skin_Cutaneous_Melanoma',
            'NA': '', np.nan: '', 0.0: '', 'Skin Cutaneous Melanoma': '',
            'Brain_Lower_Grade_Glioma': ''
        }
        for cancer, result in tcga_options.items():
            self.assertEqual(index.primary_cancer(cancer, '', ''), result)

        cesc_options = {
            'Adenosquamous': 'Cervical_Adenosquamous',
            'Cervical Squamous Cell Carcinoma':
                'Cervical_Squamous_Cell_Carcinoma',
            'Endocervical Adenocarcinoma of the Usual Type':
                'Endocervical_Adenocarcinoma',
            'Endocervical Type of Adenocarcinoma':
                'Endocervical_Adenocarcinoma',
            'Endometrioid Adenocarcinoma of Endocervix':
                'Endocervical_Adenocarcinoma',
            'Mucinous Adenocarcinoma of Endocervical Type':
                'Endocervical_Adenocarcinoma'
        }
        for option, answer in cesc_options.items():
            self.assertEqual(
                index.primary_cancer(cesc, option, ''), answer
            )

        lgg_options = {
            'Astrocytoma': 'Astrocytoma', 'NA': '',
            'Oligoastrocytoma': 'Oligoastrocytoma',
            'Oligodendroglioma': 'Oligodendroglioma'
        }
        for option, answer in lgg_options.items():
            self.assertEqual(
                index.primary_cancer(lgg, option, ''), answer
            )
        esca_options = {
            'Esophagus Adenocarcinoma, NOS': 'Esophagus_Adenocarcinoma',
            'Esophagus Squamous Cell Carcinoma':
                'Esophagus_Squamous_Cell_Carcinoma'
        }
        for option, answer in esca_options.items():
            self.assertEqual(
                index.primary_cancer(esca, '', option), answer
            )
        pcpg_options = {
            'Paraganglioma': 'Paraganglioma',
            'Paraganglioma; Extra-adrenal Pheochromocytoma': 'Paraganglioma',
            'Pheochromocytoma': 'Pheochromocytoma'
        }
        for option, answer in pcpg_options.items():
            self.assertEqual(
                index.primary_cancer(pcpg, '', option), answer
            )
        sarc_options = {
            'Dedifferentiated liposarcoma': 'Dedifferentiated_liposarcoma',
            'Desmoid Tumor': 'Desmoid_Tumor',
            "Giant cell 'MFH' / Undifferentiated pleomorphic sarcoma "
            "with giant cells":
                'Undifferentiated_Pleomorphic_Sarcoma',
            'Leiomyosarcoma (LMS)': 'Leiomyosarcoma',
            'Malignant Peripheral Nerve Sheath Tumors (MPNST)':
                'Malignant_Peripheral_Nerve_Sheath_Tumors',
            'Myxofibrosarcoma': 'Myxofibrosarcoma',
            "Pleomorphic 'MFH' / Undifferentiated pleomorphic sarcoma":
                'Undifferentiated_Pleomorphic_Sarcoma',
            'Sarcoma; synovial; poorly differentiated': 'Synovial_Sarcoma',
            'Synovial Sarcoma - Biphasic': 'Synovial_Sarcoma',
            'Synovial Sarcoma - Monophasic': 'Synovial_Sarcoma',
            'Undifferentiated Pleomorphic Sarcoma (UPS)':
                'Undifferentiated_Pleomorphic_Sarcoma'
        }
        for option, answer, in sarc_options.items():
            self.assertEqual(
                index.primary_cancer(sarc, '', option), answer
            )

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
