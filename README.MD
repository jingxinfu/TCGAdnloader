## Introduction

This package is for downloading TCGA project data on a managable data structure. 
### Notice:
The scale of TPM/FPKM range between hg19 and hg38 has hug difference.
## Install:

```
git clone git@github.com:jingxinfu/TCGAdnloader.git
cd dataCollection
pip install .
```

## Usage:
```
usage: tcgaDnloader [-h] -o OUTPUT [-r REF [REF ...]] [-c CANCER [CANCER ...]]
               [-t DATATYPE [DATATYPE ...]]

Tools to download public genomic data

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output directory (default: None)
  -r REF [REF ...], --ref REF [REF ...]
                        Aligned Reference (default: ['hg19', 'hg38'])
  -c CANCER [CANCER ...], --cancer CANCER [CANCER ...]
                        Cancer type included (default: ['BRCA'])
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
    ├── SNV
    │   ├── muse
    │   ├── mutect2
    │   ├── SomaticSnipe
    │   └── VarScan2
    └── Surv
```
### Detial
- CNV
  - broad_focal : Region reaches significance due to either broad or focal event
- RNA-Seq
  - origin:
  - zscore_paired : 
    - log2(x+1)
    - minus the average across normal samples
  - zscore_tumor : 
    - pick up tumor samples
    - log2(x+1)
    - minus the average across samples
  - fpkm_uq :
    - A modified version of the FPKM formula in which the 75th percentile read count is used as the denominator in place of the total number of protein-coding reads
  
- RPPA
  - Only available for hg19