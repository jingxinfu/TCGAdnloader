## Introduction

This package is for downloading TCGA project data on a managable data structure (see following output structure). Both of hg19 and hg38 version are provided. Please feel free to open up issues if you have any questions.

***

> *Noticce: The TPM/FPKM range between hg19 and hg38 has hug difference.* 

***



## Install:

```
git clone git@github.com:jingxinfu/TCGAdnloader.git
cd dataCollection
pip install .
```

## Usage:
```
usage: tcgaDnloader [-h] -o OUTPUT [-r REF [REF ...]] [-c CANCER [CANCER ...]]
                    [--no-meta] [--no-drug] [-t DATATYPE [DATATYPE ...]]

Tools to download public genomic data

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (default: None)
  -r REF [REF ...], --ref REF [REF ...]
                        Aligned Reference (default: ['hg19', 'hg38'])
  -c CANCER [CANCER ...], --cancer CANCER [CANCER ...]
                        Cancer type included (default: ['ACC', 'BRCA', 'CHOL',
                        'ESCA', 'KICH', 'KIRP', 'LGG', 'LUAD', 'MESO', 'PAAD',
                        'PRAD', 'SARC', 'STAD', 'TGCT', 'THYM', 'UCS', 'BLCA',
                        'CESC', 'COAD', 'DLBC', 'GBM', 'HNSC', 'KIRC', 'LAML',
                        'LIHC', 'LUSC', 'OV', 'PCPG', 'READ', 'SKCM', 'THCA',
                        'UCEC', 'UVM'])
  --no-meta             Do not install meta information (default: True)
  --no-drug             Do not install drug information (default: True)
  -t DATATYPE [DATATYPE ...], --datatype DATATYPE [DATATYPE ...]
                        Data type included (default: ['rnaseq', 'cnv', 'rppa',
                        'snv'])
```
##  Output Data Structure:
### Overview
```
(output directory)
├── hg19
│   ├── CNV
│   │   ├── all
│   │   │   └── segment
│   │   └── somatic
│   │       ├── gene
│   │       │   ├── broad_focal
│   │       │   ├── focal
│   │       │   └── threds
│   │       └── segment
│   ├── RNASeq
│   │   ├── count
│   │   │   └── origin
│   │   ├── fpkm
│   │   │   ├── origin
│   │   │   └── zscore_tumor
│   │   ├── norm_count
│   │   │   └── origin
│   │   └── tpm
│   │       ├── origin
│   │       └── zscore_tumor
│   └── RPPA
└── hg38
    ├── CNV
    │   ├── all
    │   │   └── segment
    │   └── somatic
    │       ├── gene
    │       │   └── focal
    │       └── segment
    ├── RNASeq
    │   ├── count
    │   │   └── origin
    │   ├── fpkm
    │   │   ├── origin
    │   │   └── zscore_tumor
    │   ├── fpkm_uq
    │   │   └── origin
    │   └── tpm
    │       ├── origin
    │       └── zscore_tumor
    ├── sample_pheno
    │   ├── origin
    |   ├── normal
    │   └── tumor
    ├── SNV
    │   ├── muse
    │   ├── mutect2
    │   ├── SomaticSnipe
    │   └── VarScan2
    └── Surv
```
### Note:
- CNV

  - broad_focal : Region reaches significance due to either broad or focal event.
- RNA-Seq
  - origin: Files directly downloaded by API

  - zscore_paired : 
    1. log2(x+1)

    2. minus the average across normal samples

  - zscore_tumor : 

    1. pick up tumor samples

    2. log2(x+1)

    3. minus the average across samples

  - fpkm_uq :

    - A modified version of the FPKM formula in which the 75th percentile read count is used as the denominator in place of the total number of protein-coding reads

- RPPA

  - Only available for hg19

- SNV:
  - Mutation files for hg19 can be downloaded separatly by this [url](https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc). This is a MC3 results.