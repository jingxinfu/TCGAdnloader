
import os

BASE_DIR = os.path.dirname(os.path.realpath(__file__))+'/'
# The map info was directly downloaded from
# https://gdc.xenahubs.net/download/probeMaps/gencode.v22.annotation.gene.probeMap.gz
ANNO_PATH = '/'.join([BASE_DIR, 'data/gencode.v22.annotation.gene.probeMap'])

Biospecimen_INFO = {
    "slide": [
        'bcr_sample_barcode', 'percent_lymphocyte_infiltration',
        'percent_monocyte_infiltration', 'percent_necrosis',
        'percent_neutrophil_infiltration',
        'percent_normal_cells', 'percent_stromal_cells',
        'percent_tumor_cells', 'percent_tumor_nuclei'
              ],
    "auxilary": [
        'bcr_patient_barcode', 'hpv_status',
        'mononucleotide_and_dinucleotide_marker_panel_analysis_status'
                 ]
}
Biospecimen_MAP ={
    "slide": {'bcr_sample_barcode':'patient'},
    "auxilary": {
        'bcr_patient_barcode': 'patient',
        'mononucleotide_and_dinucleotide_marker_panel_analysis_status':'MSI_Status'
    }
}

CLIN_INFO = {
    "patient": [
        'bcr_patient_barcode', 'gender','race',
        'ajcc_pathologic_tumor_stage',
        'histologic_diagnosis',
    ],
    "survival": [
        'bcr_patient_barcode', 'vital_status',
        'days_to_death','days_to_last_followup',
        'age_at_initial_pathologic_diagnosis'
        
    ],

}


# CANCER_LIST = [
#     'ACC','BRCA','CHOL','ESCA','KICH','KIRP','LGG','LUAD','MESO','PAAD','PRAD',
#     'SARC','STAD','TGCT','THYM','UCS','BLCA','CESC','COAD','DLBC','GBM','HNSC',
#     'KIRC','LAML','LIHC','LUSC','OV','PCPG','READ','SKCM','THCA','UCEC','UVM',
# ]
CANCER_LIST = ['ACC']

# RNASEQ_DIR = '../RNASeq'
# CNV_DIR = '../CNV'
# RPPA_DIR = '../RPPA'
