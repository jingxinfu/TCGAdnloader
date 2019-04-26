
import os

BASE_DIR = os.path.dirname(os.path.realpath(__file__))+'/'
# The map info was directly downloaded from
# https://gdc.xenahubs.net/download/probeMaps/gencode.v22.annotation.gene.probeMap.gz
ANNO_PATH = '/'.join([BASE_DIR, 'data/gencode.v22.annotation.gene.probeMap'])

# TCGA BRCA Article (Cancer Genome Atlas Network. Nature, 2012)
## Supplementary files: http://www.nature.com/nature/journal/v490/n7418/extref/nature11412-s2.zip
### SuppTable1, contains TCGA patients with PAM50
PAM50_PATH = '/'.join([BASE_DIR, 'data/PAM50.txt'])

# TCGA BRCA Article (Ciriello G, Gatza ML, Beck AH, et al. Cell. 2015)
## Supplementary files: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4603750/bin/NIHMS724218-supplement-2.xlsx
### Suppl. Table 1: Triple negative: All of ['ER IHC','PR IHC','HER2 IHC'] are negative


Biospecimen_INFO = {
    "histology": {
        "slide":[
                'bcr_sample_barcode', 'percent_lymphocyte_infiltration',
                'percent_monocyte_infiltration', 'percent_necrosis',
                'percent_neutrophil_infiltration',
                'percent_normal_cells', 'percent_stromal_cells',
                'percent_tumor_cells', 'percent_tumor_nuclei'
              ],
        
    },
    "phenotype":{
        "sample": ['bcr_sample_barcode', 'sample_type'],
        "auxilary": [
                'bcr_patient_barcode', 'hpv_status',
                'mononucleotide_and_dinucleotide_marker_panel_analysis_status'
                ],
        
        "patient": ["bcr_patient_barcode",'race']
    }
    
}

Biospecimen_MAP ={
    'bcr_sample_barcode': 'sample',
    'bcr_patient_barcode': 'patient',
    'mononucleotide_and_dinucleotide_marker_panel_analysis_status':'MSI_Status'
}

CLIN_INFO = {
    "submitter_id",
    # "diagnoses.primary_diagnosis",
    "diagnoses.age_at_diagnosis",
    "diagnoses.vital_status",
    "diagnoses.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.tumor_stage",
    "demographic.gender",
    # "demographic.race",
    # "demographic.ethnicity"
}

CLIN_MAP = {
    'submitter_id': 'patient',
    'diagnoses.0.tumor_stage': 'stage',
    'diagnoses.0.vital_status': 'OS_Event',
    'diagnoses.0.age_at_diagnosis': 'age',
    'diagnoses.0.days_to_death':'OS_D',
    'diagnoses.0.days_to_last_follow_up':'OS_F',
    'demographic.gender':'gender'
}


DRUG_MAP = {
    "bcr_patient_barcode": "patient",
    "bcr_drug_barcode":"drug_barcode",
    "pharmaceutical_therapy_drug_name":"drug_name",
    "pharmaceutical_therapy_type":"drug_type",
    "pharmaceutical_tx_started_days_to":"started_days",
    "pharmaceutical_tx_ongoing_indicator": "ongoing_indicator",
    "pharmaceutical_tx_ended_days_to":"ended_days",
    "pharmaceutical_tx_dose_units":"dose_units",
    "pharmaceutical_tx_total_dose_units": "total_dose_units",
    "prescribed_dose": "prescribed_dose",
    "regimen_number": "regimen_number",
    "total_dose":"total_dose"
}

CANCER_LIST = [
    'ACC','BRCA','CHOL','ESCA','KICH','KIRP','LGG','LUAD','MESO','PAAD','PRAD',
    'SARC','STAD','TGCT','THYM','UCS','BLCA','CESC','COAD','DLBC','GBM','HNSC',
    'KIRC','LAML','LIHC','LUSC','OV','PCPG','READ','SKCM','THCA','UCEC','UVM',
]

# CLIN_VERSION = {
#     'ACC': '4.0', 'BRCA': '4.0', 'CHOL': '4.0', 'ESCA': '4.0', 'KICH': '4.4',
#     'KIRP': '1.0', 'LGG': '1.0', 'LUAD': '1.0', 'MESO': '4.0', 'PAAD': '4.4',
#     'PRAD': '1.0', 'SARC': '4.0', 'STAD': '1.0', 'TGCT': '4.0', 'THYM': '4.0',
#     'UCS': '4.0', 'BLCA': '4.0', 'CESC': '4.0', 'COAD': '1.0', 'DLBC': '4.4',
#     'GBM': '1.0', 'HNSC': '4.8', 'KIRC': '1.0', 'LAML': None, 'LIHC': '4.0',
#     'LUSC': '1.0', 'OV': '1.0', 'PCPG': '4.0', 'READ': '1.0', 'SKCM': '2.0',
#     'THCA': '4.0', 'UCEC': '4.0', 'UVM': '4.0',
# }


