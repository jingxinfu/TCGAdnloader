# !/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from .setting import ANNO_PATH

def rmEntrez(df):
    ''' Format gene expression profile with gene symbol index
    Specific function for expression data dowloaded from firebrowse.
    
    Parameters
    ----------
    df : pandas.DataFrame
        df with genes in rows and samples in columns.
        The row index is like: AATK|9625.

    Returns
    -------
    pandas.DataFrame
        Same format as input, but modified index. 
        Turn AATK|9625 --> AATK
    '''

    df.index = df.index.map(lambda x: x.split('|')[0])
    df = df.iloc[~ df.index.isin(['?','']),:]
    df = df.groupby(level=0).mean()
    return df

def mapEm2Gene(df, anno_path=ANNO_PATH):
    name_map = pd.read_table(anno_path, index_col=0)['gene'].to_dict()
    df.rename(index=name_map, inplace=True)
    df = df.groupby(level=0).mean()
    return df
    
def pick(df, source='tumor'):
    source_map = {'tumor': '0[0-9]$',
            'normal': '1[0-9]$'}

    if not source in source_map.keys():
        raise KeyError("""
        {0} is not a valid type of source, only accept following input: {1}
        """.format(source, ','.join(source_map.keys())))
    
    return df.loc[:, df.columns.str.contains(source_map[source])]

def calTNzcore(df,pair_TN = True):
    ''' Ways to normalize tumor expression data on TCGA
    
    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with gene indentifier on rows and samples on columns.
        Value on the data frame ARE log scaled!
    pair_TN : bool, optional
        Tell whether normalize tumor expression profile by paired tumor samples (the default is True)
    
    Raises
    ------
    ValueError
        cannot find enough paired normal samples
    
    Returns
    -------
    pandas.DataFrame
        Normlized tumor gene expression profile
    '''

    tumor = pick(df, source='tumor')

    if pair_TN:
        normal = pick(df, source='normal')
        if normal.shape[1] > 30:
            norm_factor = normal.mean(axis=1)
        else:
            raise ValueError('Cannot find enough paired normal samples (<10)')
    else:
         norm_factor = tumor.mean(axis=1)

    result = tumor.subtract(norm_factor, axis=0)
    result.columns = result.columns.map(lambda x: '-'.join(x.split('-')[:3]))
    return result


def mergeSampleToPatient(df):
    '''
    A inplace function.
    Changes samples level profile into patient level, but keep tumor and normal information.
    Data from the same sample but from different vials/portions/analytes/aliquotes is averaged.

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


