# !/usr/bin/env python3

import pandas as pd
import numpy as np
import os,re
from .setting import ANNO_PATH, CLIN_MAP


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
    
def pick(df, source='tumor',transpose=False):
    ''' Picking samples by where it dissect from. 
    
    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with gene indentifier on rows and samples on columns.
    source : str, optional
        The region to pick (the default is 'tumor')
    transpose : bool, optional
        Whether the samples id is on the index (the default is False)
    
    Raises
    ------
    KeyError
        if input a invalid data type
    
    Returns
    -------
    pandas.DataFrame
        A data frame with gene indentifier on rows and picked samples on columns.
    '''

    source_map = {'tumor': '0[0-9]$',
            'normal': '1[0-9]$'}

    if not source in source_map.keys():
        raise KeyError("""
        {0} is not a valid type of source, only accept following input: {1}
        """.format(source, ','.join(source_map.keys())))

    if transpose:
        return df.loc[df.index.str.contains(source_map[source]),:]
    else:
        return df.loc[:, df.columns.str.contains(source_map[source])]

def calTNzcore(df,pair_TN = True):
    ''' Ways to normalize tumor expression data on TCGA
    
    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with gene indentifier on rows and samples on columns.
        NOTICE: Value on the data frame ARE already been log scaled!

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

    result = tumor.sub(norm_factor, axis='index')
    result = result.div(result.std(axis=1), axis='index')
    
    return result


def mergeToSample(df,transpose=False):
    '''
    Changes aliquot level profile into sample level, but keep tumor and normal information.
    Data from the same sample but from different vials/portions/analytes/aliquotes is averaged.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame index by gene/entrez_id, column by TCGA barcode
    transpose : bool (option)
        Whether TCGA barcode is on index.(default is False)

    Inplace
    -------
    pandas.DataFrame
        Data frame index by short TCGA barcode.
        Eg.
        TCGA-OR-A5J1-01 tumor
        TCGA-OR-A5J1-11 normal
    '''
    if transpose is True:
        df.index = df.index.map(lambda x: '-'.join(x.split('-')[:4])[:-1])
        df = df.groupby(level=0).mean()
    else:
        df.columns = df.columns.map(lambda x: '-'.join(x.split('-')[:4])[:-1])
        df = df.T.groupby(level=0).mean().T
    return df

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


def formatClin(df):
        
    df['OS_Event'] = df['OS_Event'].map({'dead': 1, 'alive': 0})
    df['OS'] = df[['OS_D', 'OS_F', ]].max(axis=1)
    df = df.drop(['OS_D', 'OS_F'],axis=1)

    df['gender'] = pd.to_numeric(df['gender'].map({'female': 1, "male": 2}))
    df['stage'] = df['stage'].map(lambda x: re.sub(
        '[a-f]$', '', x) if isinstance(x, str) else np.nan)
    df['stage'] =df['stage'].map(
                            {
                                'stage i': 1,
                                'stage ii': 2,
                                'stage iii': 3,
                                'stage iv':4,
                                'stage v': 5,
                                'stage vi': 6,
                                'stage vii': 7,
                                'stage viii': 8,
                                'stage ix': 9,
                                'stage x': 10,
                            }
                            )

    return df.groupby('patient').max()
