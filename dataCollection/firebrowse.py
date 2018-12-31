#!/usr/bin/env python3

import subprocess,os
import pandas as pd
from .convertor import mergeSampleToPatient, segmentToGene
from .setting import CANCER_LIST

# TODO compared gistic result with xenas



def fget(cancer, data_type, store_dir,
        base_url="http://gdac.broadinstitute.org/runs",
        release_time="2016_01_28"):
    ''' Download level 3 data from FireBrowse

    Parameters
    ----------
    cancer : str
        Cancer type included in TCGA project
    data_type : str
        Level 3 data type provided by FireBrowse
    store_dir : str
        Output directory
    base_url : str, optional
        URL prefix (the default is "http://gdac.broadinstitute.org/runs", which is the prefix provided by FireBrowse)
    release_time : str, optional
        Release version and this release recored by date. (the default is "2016_01_28", which is the latest available release for now.)

    Raises
    ------
    KeyError
        if the input parameter is out of provided list.

    Returns
    -------
    str
        Run messages. Return 'Success' if no error occurs.
    '''

    data_type_dict = {
        "rna_raw" : "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3",
        "rna_norm": "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3",
        "rppa": "RPPA_AnnotateWithGene.Level_3",
        "cnv_all": "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3",
        "cnv_minus_germline": "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
    }
    keep_suffix_dict = {
        "rna_raw": "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
        "rppa" : "rppa.txt",
        "rna_norm": "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
        "cnv_all": "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt",
        "cnv_minus_germline": "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
    }

    if not data_type in data_type_dict.keys():
        raise KeyError("""
        {0} is not a valid data type, only accept following input: {1}
        """.format(data_type,','.join(data_type_dict.keys())))

    short_release_time = "".join(release_time.split('_'))

    release = "stddata__{release_time}"
    sub_folder = "data/{cancer}/{short_release_time}"
    file_name = "gdac.broadinstitute.org_{cancer}.{data_type}.{short_release_time}00.0.0.tar.gz"

    url = "/".join([base_url,release,sub_folder,file_name])
    url = url.format(**dict(
        cancer = cancer,
        data_type = data_type_dict[data_type],
        release_time = release_time,
        short_release_time = short_release_time,
        )
    )

    cmd ="""
    set -x
    [ -d {store_dir}/{cancer}_tmp ] || mkdir -p {store_dir}/{cancer}_tmp
    wget -v -O {store_dir}/{cancer}.gz {url}
    tar -xvvf {store_dir}/{cancer}.gz -C {store_dir}/{cancer}_tmp --strip-components=1
    rm {store_dir}/{cancer}.gz
    mv {store_dir}/{cancer}_tmp/*{keep_suffix} {store_dir}/{cancer}
    rm -rf {store_dir}/{cancer}_tmp
    """.format(**dict(
        store_dir=store_dir,
        cancer=cancer,
        url=url,
        keep_suffix=keep_suffix_dict[data_type]
        )
    )

    try:
        subprocess.call(cmd,shell=True)
        log = 'Success'
    except subprocess.CalledProcessError as e:
        log = e

    return log


def splitCountTPM(raw_rnaseq_path):
    ''' Split one data frame with both count and scaled_estiamte into two data frames and
    merge the sample level data frame into pateint level data frame, but keep separating tumor and normal samples.
    Then, based on the scaled_estimate column, calculate TPM and RPKM information.

    Parameters
    ----------
    raw_rnaseq_path : str
        Path to raw rnaseq data download from FireBrowse

    Returns
    -------
    Dict
        A dict that contains three pandas.DataFrame, which are raw count, TPM and RPKM.
        All of those data frame are index by both Entrez ID and gene symbol and colum named by four digits TCGA barcode.
    '''


    df = pd.read_table(raw_rnaseq_path, index_col=0,skiprows=[1])
    col_selector = pd.read_table(raw_rnaseq_path, index_col=0, nrows=2)

    raw_count = df.loc[:, col_selector.iloc[0, :] =='raw_count']
    raw_count = round(raw_count)
    mergeSampleToPatient(raw_count)

    transcipt_fraction = df.loc[:,col_selector.iloc[0, :] == 'scaled_estimate']
    mergeSampleToPatient(transcipt_fraction)

    tpm = transcipt_fraction * 10e6
    normalize_factor = transcipt_fraction.sum(axis=0)
    fpkm = transcipt_fraction * normalize_factor * 10e9

    return dict(count=raw_count,tpm=tpm,fpkm=fpkm)


def storeData(df, parental_dir, sub_folder, cancer):
    store_folder = '/'.join([parental_dir, sub_folder])
    if not os.path.exists(store_folder):
            os.makedirs(store_folder)

    df.to_csv('/'.join([store_folder, cancer]), sep='\t')


def rnaseqWorkflow(parental_dir, cancer):
    '''
    Workflow for downloading RNAseq data from FireBrowse and preprocessing data format.

    Parameters
    ----------
    parental_dir : str
        Path to parental folder that you want to store the whole RNAseq data
    cancer : str
        Cancer name you want to download from FireBrowse, it must be a cancer type included in TCGA project.

    '''


    ########################## Raw count and Scale Estimate ##########################
    # 1. Fetch raw count and RSEM information from FireBrowse
    # 2. Split fetched data frame into raw count and RSEM separatly.
    # 3. Merge sample level data into pateint level data, but still separate tumor and normal sample.
    # 4. Calculate TPM and RPKM based on RSEM results.
    ##################################################################################

    log = fget(
        cancer=cancer, data_type='rna_raw',
        store_dir='/'.join([parental_dir, 'raw'])
        )

    if log != 'Success':
        return 'rna_raw:    '+log

    raw_rnaseq = splitCountTPM(
        raw_rnaseq_path='/'.join([parental_dir, 'raw', cancer])
        )
    for name, df in raw_rnaseq.items():
        storeData(df=df, parental_dir=parental_dir,
                  sub_folder=name, cancer=cancer)

    subprocess.call('rm -rf {}'.format('/'.join([parental_dir, 'raw'])), shell=True)
    ########################## Raw count and Scale Estimate ##########################
    # 1. Fetch normalized count from FireBrowse
    # 2. remove the second row, which only  indicate the normalized count
    # 3. Merge sample level data into pateint level data, but still separate tumor and normal sample.
    ##################################################################################

    log = fget(
        cancer=cancer,data_type='rna_norm',
        store_dir='/'.join([parental_dir,'norm'])
        )

    if log != 'Success':
        return 'rna_norm:    '+log

    rnaseq_norm = pd.read_table('/'.join([parental_dir, 'norm', cancer]), index_col=0, skiprows=[1])
    mergeSampleToPatient(rnaseq_norm)
    storeData(df=rnaseq_norm, parental_dir=parental_dir,
              sub_folder='norm_count', cancer=cancer)

    subprocess.call('rm -rf {}'.format('/'.join([parental_dir, 'norm'])), shell=True)
    return 'Success'

def cnvWorkflow(parental_dir,cancer):
    '''
    Workflow for downloading copy number variation data from FireBrowse and preprocessing data format.

    Parameters
    ----------
    parental_dir : str
        Path to parental folder that you want to store the whole copy number variation data
    cancer : str
        Cancer name you want to download from FireBrowse, it must be a cancer type included in TCGA project.
    '''

    log = fget(
        cancer=cancer, data_type='cnv_minus_germline',
        store_dir='/'.join([parental_dir, 'cnv_minus_germline','segment'])
        )

    if log != 'Success':
        return 'cnv_minus_germline:    '+log

    cnv_gene = segmentToGene(
        sgment_path='/'.join([parental_dir, 'cnv_minus_germline','segment',cancer]),
        ref_version='hg19',
    )
    parental_dir = '/'.join([parental_dir,'gene'])
    for name, df in cnv_gene.items():
        storeData(
            df=df, parental_dir=parental_dir,
            sub_folder=name, cancer=cancer
            )

    return 'Success'

def rppaWorkflow(parental_dir,cancer):
    '''
    Workflow for downloading RPPA data from FireBrowse and preprocessing data format.

    Parameters
    ----------
    parental_dir : str
        Path to parental folder that you want to store the whole RPPA data
    cancer : str
        Cancer name you want to download from FireBrowse, it must be a cancer type included in TCGA project.
    '''

    log = fget(
        cancer=cancer,data_type='rppa',
        store_dir=parental_dir
    )

    if log != 'Success':
        return 'RPPA:    '+log

    rppa = pd.read_table('/'.join(parental_dir, cancer), index_col=0)
    mergeSampleToPatient(rppa)
    rppa.to_csv('/'.join([parental_dir, cancer]), sep='\t')

    return 'Success'


def main():
    RNASEQ_DIR = '../RNASeq'
    CNV_DIR = ''
    RPPA_DIR = ''
    log_out = pd.Series(
        [rnaseqWorkflow(parental_dir=RNASEQ_DIR, cancer=cancer)
         for cancer in CANCER_LIST],
        index=CANCER_LIST)
    log_out.to_frame().to_csv('/'.join(RNASEQ_DIR, 'log.txt'))


if __name__ == "__main__":
    main()
