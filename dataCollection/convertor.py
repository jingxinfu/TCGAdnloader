# !/usr/bin/env python3

import pandas as pd

def mergeSampleToPatient(df):
    '''
    A inplace function.
    Changes samples level profile into patient level, but keep tumor and normal information.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame index by gene/entrez_id, column by TCGA barcode

    Inplace
    -------
    pandas.DataFrame
        Data frame index by short TCGA barcode.
        Eg.
        TCGA-OR-A5J1-01 tumor
        TCGA-OR-A5J1-11 normal
    '''
    df.columns = df.columns.map(lambda x: '-'.join(x.split('-')[:4])[:-1])
    df = df.T.groupby(level=0).mean().T

def tpmToFpkm(df,reverse=False):
    ''' Conversion between TPM and FPKM/RPKM

    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with variable(e.g. gene) on rows and observarion on columns.
    reverse : bool, optional
        Convert direction.
        True: TPM --> FPKM/RPKM
        False: FPKM/RPKM --> TPM
        (the default is False, which means convert TPM to FPKM/RPKM)


    Returns
    -------
    pandas.DataFrame
        A data frame with variable(e.g. gene) on rows and observarion on columns.
    '''

    row_sum = df.sum(axis=0)

    if reverse:
        tpm = (df/row_sum) * 10e6
        return tpm
    else:
        fpkm = (df * row_sum) / 10e3
        return fpkm


def segmentToGene(sgment_path, ref_version='hg19'):
    # TODO segment to gene level, figure out the parameter meaning.
    pass
