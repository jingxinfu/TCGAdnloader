
import os

BASE_DIR = os.path.dirname(os.path.realpath(__file__))+'/'
# The map info was directly downloaded from
# https://gdc.xenahubs.net/download/probeMaps/gencode.v22.annotation.gene.probeMap.gz
ANNO_PATH = '/'.join([BASE_DIR, 'data/gencode.v22.annotation.gene.probeMap'])


# CANCER_LIST = [
#     'ACC','BRCA','CHOL','ESCA','KICH','KIRP','LGG','LUAD','MESO','PAAD','PRAD',
#     'SARC','STAD','TGCT','THYM','UCS','BLCA','CESC','COAD','DLBC','GBM','HNSC',
#     'KIRC','LAML','LIHC','LUSC','OV','PCPG','READ','SKCM','THCA','UCEC','UVM',
# ]
CANCER_LIST = ['ACC']

# RNASEQ_DIR = '../RNASeq'
# CNV_DIR = '../CNV'
# RPPA_DIR = '../RPPA'
