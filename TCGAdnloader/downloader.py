#!/usr/bin/env python3

import subprocess, os,time
import pandas as pd
import numpy as np
from .convertor import mergeSampleToPatient, calTNzcore, rmEntrez, tpmToFpkm, mapEm2Gene, formatClin, pick
from .outformat import storeData
import requests,json,re,io
from .setting import CLIN_INFO, Biospecimen_INFO, Biospecimen_MAP, PAM50_PATH, CLIN_VERSION

class GdcApi(object):
    ''' 
    API for download files from GDC
    '''

    __slot__ = ["files_endpt", "data_endpt", "cancer", "parental_dir"]

    def __init__(self, cancer, parental_dir,data_endpt="https://api.gdc.cancer.gov/data", files_endpt="https://api.gdc.cancer.gov/files", **kwargs):
        ''' Intialize instance parameters
        
        Parameters
        ----------
        cancer : str
            Cancer type
        parental_dir : str
            Path to store datas
        data_endpt : str, optional
            [Endpoint for files id searching] (the default is "https://api.gdc.cancer.gov/data")
        files_endpt : str, optional
            [Endpoint for files downloading] (the default is "https://api.gdc.cancer.gov/files")
        
        '''

        self.files_endpt = files_endpt
        self.data_endpt = data_endpt
        self.cancer = cancer
        self.parental_dir = parental_dir

    def _projFilter(self, data_type):
        dtype_dict = {
            "cnv_segment_somatic": "Masked Copy Number Segment",
            "cnv_segment_all": "Copy Number Segment",
        }
        filters = {
            "op": "and",
            "content":[
                {
                    "op": "in",
                    "content": {
                        "field": "files.data_type",
                        "value": [
                            dtype_dict[data_type]
                        ]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": [
                            "TCGA-"+self.cancer.upper()
                        ]
                    }
                },
            ]
        }

        params = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "size": "3000"
        }

        return params

    def _nameFilter(self, data_type):
        dtype_dict = {
            'gistic': '{}.focal_score_by_genes.txt'.format(self.cancer.upper()),
            'survival': "nationwidechildrens.org_clinical_follow_up_v{0}_{1}.txt".format(CLIN_VERSION[self.cancer], self.cancer.lower()),
            'patient': "nationwidechildrens.org_clinical_patient_{}.txt".format(self.cancer.lower()),
            'aliquot': "nationwidechildrens.org_biospecimen_aliquot_{}.txt".format(self.cancer.lower()),
            'slide': "nationwidechildrens.org_biospecimen_slide_{}.txt".format(self.cancer.lower()),
            'sample': "nationwidechildrens.org_biospecimen_sample_{}.txt".format(self.cancer.lower()),
            'auxilary': "nationwidechildrens.org_auxiliary_{}.txt".format(self.cancer.lower()),
        }
        filters = {
            "op": "in",
            "content": {
                "field": "files.file_name",
                "value": [
                    dtype_dict[data_type]
                ]
            }
        }

        params = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "size": "1"
        }

        return params

    def _fetchFileID(self, data_type, by_name=True):
        ''' Get files id by upstream filter parameters
        
        Parameters
        ----------
        data_type : str
            Data type to be download. eg. gistic
        by_name : bool, optional
            Whether getting files id by matching file names (the default is True).
            If not, we will use project filtering options to get file id list.
        
        Returns
        -------
        list
            A list contains file ids.
        '''

        if by_name is True:
            
            file_uuid_list = []
            params = self._nameFilter(data_type)
            response = requests.get(self.files_endpt, params=params)

            for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
                file_uuid_list.append(file_entry["file_id"])

        else:
            file_uuid_list = []
            params = self._projFilter(data_type)
            response = requests.get(self.files_endpt, params=params)
            for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
                file_uuid_list.append(file_entry["file_id"])

        if len(file_uuid_list) == 0:
            return None,'Not found'
        else:
            return file_uuid_list,None
       
    def getTable(self, data_type, by_name=True, **kwargs):
       
        ''' 
        Merging tables downloaded by a list of file ids

        '''
        try:
            file_uuid_list, error = self._fetchFileID(
                data_type=data_type, by_name=by_name)
        except requests.exceptions.SSLError:
            time.sleep(10)
            file_uuid_list, error = self._fetchFileID(
                data_type=data_type, by_name=by_name)

        if error != None:
            return None, error
        ready_to_merge = []

        for ids in file_uuid_list:
            params = {"ids": [ids]}
            try:
                response = requests.post(self.data_endpt, data=json.dumps(
                    params), headers={"Content-Type": "application/json"})
            except requests.exceptions.SSLError:
                time.sleep(10)
                response = requests.post(self.data_endpt, data=json.dumps(
                    params), headers={"Content-Type": "application/json"})

            df = pd.read_table(io.StringIO(
                response.content.decode("utf-8")), **kwargs)
            ready_to_merge.append(df)
        
        return pd.concat(ready_to_merge,axis=0),None

    def clin(self):
        '''
        Downloading clinical information
        '''

        read_to_merge=[]
        for k,v in CLIN_INFO.items():
            meta,_ = self.getTable(data_type=k)
            # The reason why use the second raw as row name is that
            # the name of follow_up info columns is consistent with different version
            meta.columns = meta.iloc[0,:]
            meta = meta.iloc[2:,]
            meta = meta[meta.columns.intersection(v)]
            non_info = pd.Index(v).difference(meta.columns)
            
            for c in non_info:
                meta[c] = '[Not Applicable]'
            read_to_merge.append(meta)

        basic_clin = pd.concat(read_to_merge,join='outer',axis=0,sort=True)
        basic_clin = formatClin(basic_clin)
        basic_clin.index.name = 'patient'

        storeData(basic_clin,parental_dir=self.parental_dir,
                  sub_folder='Surv',cancer=self.cancer)

    def biospecimen(self):
        '''
        Downloading  biopecimen information
        '''
        for sub_folder,files in Biospecimen_INFO.items():
            read_to_merge = []
            for k, v in files.items():
                meta, errors = self.getTable(data_type=k)
                if errors == None:
                    meta = meta[meta.columns.intersection(v)]
                    non_info = pd.Index(v).difference(meta.columns)
                    for c in non_info:
                        meta[c] = np.nan

                    meta.replace('[Not Available]', np.nan, inplace=True)
                    meta.replace('[Not Applicable]', np.nan, inplace=True)
                    meta.rename(columns=Biospecimen_MAP,inplace=True)

                    ## header process
                    if 'bcr_sample_barcode' in v:
                        meta.rename(
                            columns={'bcr_sample_barcode': 'patient'}, inplace=True)
                        meta = meta.drop(0, axis=0).set_index('patient')
                       
                    elif 'hpv_status' in v:
                        meta = meta.drop(0,axis=0).set_index('patient')
                    else:
                        meta = meta.drop([0,1],axis=0).set_index('patient')
                 
                    ## additional info
                    if k == 'slide':
                        meta = meta.apply(pd.to_numeric)
                        meta = mergeSampleToPatient(meta,transpose=True)

                    if k == "patient" and self.cancer == 'BRCA':
                        pam50 = pd.read_table(PAM50_PATH, index_col=0).rename(columns={
                            "PAM50 mRNA":'Subtype'})['Subtype'].to_frame()
                        meta = meta.merge(pam50, left_index=True,right_index=True,how='left')

                    read_to_merge.append(meta)


            result = pd.concat(read_to_merge, axis=1, join='outer', sort=True)

            ## Store tumor and normal info separatelly
            if sub_folder == "sample_pheno":
                for s in ['tumor','normal']:
                    sub_result = pick(result, source=s, transpose=True)
                    sub_result.index = sub_result.index.map(
                        lambda x: '-'.join(x.split('-')[:3]))
                        
                    storeData(sub_result,
                              parental_dir=self.parental_dir,
                              sub_folder='/'.join([sub_folder,s]), cancer=self.cancer)
                   
                sub_folder += '/origin'

            result.index.name = 'patient'

            storeData(pd.concat(read_to_merge, axis=1,join='outer',sort=True),
                     parental_dir=self.parental_dir,
                     sub_folder=sub_folder,cancer=self.cancer)



class Workflow(object):
    __slot__ = ['cancer', 'parental_dir', 'workflow']
    def __init__(self,cancer,parental_dir,workflow):
        self.cancer = cancer
        self.parental_dir = parental_dir
        self.workflow = workflow


    def run(self):
        for n in self.workflow:
            self.__getattribute__(n)()
            

class FireBrowseDnloader(Workflow):
    __slot__ = ['release_time']
    def __init__(self, release_time="2016_01_28", base_url="http://gdac.broadinstitute.org/runs",**kwargs):
        super(FireBrowseDnloader, self).__init__(**kwargs)
        self.release_time = release_time
        self.base_url = base_url
        
    def _fget(self,data_type, store_dir):
        
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
        # modifition to adapt CNV data on the function
        if data_type == 'cnv_gene_somatic':
            release_prefix = 'analyses'
            cancer_suffix = '-TP'
            if self.cancer == 'SKCM':
                cancer_suffix = '-TM'
        else:
            cancer_suffix = ''
            release_prefix = 'stddata'

        data_type_dict = {
            "rna_raw" : "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3",
            "rna_norm": "Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3",
            "rppa": "RPPA_AnnotateWithGene.Level_3",
            "cnv_gene_somatic": "CopyNumber_Gistic2.Level_4",
            "cnv_segment_somatic": "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
            "cnv_segment_all": "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3",
        }
        keep_suffix_dict = {
            "rna_raw": "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
            "rppa" : "rppa.txt",
            "rna_norm": "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
            "cnv_gene_somatic": "by_genes.txt",
            "cnv_segment_somatic": "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
            "cnv_segment_all": "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt",
        }


        if not data_type in data_type_dict.keys():
            raise KeyError("""
            {0} is not a valid data type, only accept following input: {1}
            """.format(data_type,','.join(data_type_dict.keys())))

        short_release_time = "".join(self.release_time.split('_'))
    
    
        release = release_prefix+"__{release_time}" 

        sub_folder = "data/{cancer}/{short_release_time}"
        file_name = "gdac.broadinstitute.org_{cancer}.{data_type}.{short_release_time}00.0.0.tar.gz"

        url = "/".join([self.base_url, release, sub_folder, file_name])
        
        url = url.format(**dict(
            cancer=self.cancer+cancer_suffix,
            data_type=data_type_dict[data_type],
            release_time=self.release_time,
            short_release_time=short_release_time,
            )
        )
        cmd ="""
        set -x
        [[ -d {store_dir}_{cancer}_{data_type}_tmp ]] || mkdir -p {store_dir}_{cancer}_{data_type}_tmp
        wget -q -O {store_dir}_{cancer}_{data_type}.gz {url}
        tar -xvvf {store_dir}_{cancer}_{data_type}.gz -C {store_dir}_{cancer}_{data_type}_tmp --strip-components=1
        rm {store_dir}_{cancer}_{data_type}.gz
        if [ $(ls {store_dir}_{cancer}_{data_type}_tmp/*{keep_suffix}| wc -l) -gt 1 ];then
            [[ -d {store_dir}_{cancer} ]] || mkdir {store_dir}_{cancer}
        fi
        mv {store_dir}_{cancer}_{data_type}_tmp/*{keep_suffix} {store_dir}_{cancer}
        rm -rf {store_dir}_{cancer}_{data_type}_tmp
        """.format(**dict(
            store_dir=store_dir,
            cancer=self.cancer,
            url=url,
            keep_suffix=keep_suffix_dict[data_type],
            data_type=data_type,
            )
        )

        try:
            subprocess.call(cmd,shell=True)
            log = 'Success'
        except subprocess.CalledProcessError as e:
            log = e

        return log

    def _splitCountTPM(self, raw_rnaseq_path):

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
        raw_count = mergeSampleToPatient(raw_count)
        raw_count = round(raw_count)

        ## Get fpkm and tpm information from transcript fractions
        transcipt_fraction = df.loc[:,col_selector.iloc[0, :] == 'scaled_estimate']
        tpm = transcipt_fraction * 10e6
        normalize_factor = transcipt_fraction.sum(axis=0)
        fpkm = transcipt_fraction * normalize_factor * 10e9

        tpm = mergeSampleToPatient(tpm)
        fpkm = mergeSampleToPatient(fpkm)

        return dict(count=raw_count,tpm=tpm,fpkm=fpkm)

    
    def _formatGistic(self, gistic_path):
        ''' Formating GISTIC results and sepratate files into segment and gene level
        
        Parameters
        ----------
        gistic_path : str
            Path to the folder of gistic output
        
        Returns
        -------
        dict
            Dictionary with files output name as key and pandas.DataFrame as value 
        '''

        f_dict = {
            "broad_focal": '{}/all_data_by_genes.txt',
            "focal": '{}/focal_data_by_genes.txt',
            "threds": '{}/all_thresholded.by_genes.txt'
        }
        result = {}
        for k, v in f_dict.items():
            if os.path.isfile(v.format(gistic_path)):
                result[k] = pd.read_table(v.format(gistic_path),index_col=0).drop(['Locus ID', 'Cytoband'],axis=1)
        
        return result

    def rnaseq(self):
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
        store_dir = '/'.join([self.parental_dir, 'RNASeq'])
        store_dir_raw = '_'.join([store_dir, 'raw'])
        store_dir_norm = '_'.join([store_dir, 'norm'])

        log = self._fget(data_type='rna_raw',store_dir=store_dir_raw)

        if log != 'Success':
            return 'rna_raw:    '+log

        raw_rnaseq = self._splitCountTPM(
            raw_rnaseq_path='_'.join([store_dir_raw, self.cancer])
            )
        for name, df in raw_rnaseq.items():
            df = rmEntrez(df)
            if name in ['fpkm','tpm']:
                log_df = np.log2( 1+ df )
                tumor_zscore = calTNzcore(log_df, pair_TN=False)
                storeData(df=tumor_zscore, parental_dir=store_dir,
                          sub_folder=name+'/zscore_tumor/', cancer=self.cancer)
                try:
                    paired_zscore = calTNzcore(log_df, pair_TN=True)
                    storeData(df=paired_zscore, parental_dir=store_dir,
                            sub_folder=name+'/zscore_paired/', cancer=self.cancer)
                except ValueError:
                    pass

            name += '/origin'
            storeData(df = df, parental_dir = store_dir,
                    sub_folder=name, cancer=self.cancer)

        subprocess.call(
            'rm -rf {}'.format('_'.join([store_dir_raw, self.cancer])), shell=True)

        ########################## Raw count and Scale Estimate ##########################
        # 1. Fetch normalized count from FireBrowse
        # 2. remove the second row, which only  indicate the normalized count
        # 3. Merge sample level data into pateint level data, but still separate tumor and normal sample.
        ##################################################################################

        log = self._fget(data_type='rna_norm',store_dir=store_dir_norm)

        if log != 'Success':
            return 'rna_norm:    '+log

        rnaseq_norm = pd.read_table(
            '_'.join([store_dir_norm, self.cancer]), index_col=0, skiprows=[1])
        rnaseq_norm = mergeSampleToPatient(rnaseq_norm)
        rnaseq_norm = rmEntrez(rnaseq_norm)
        storeData(df=rnaseq_norm, parental_dir=store_dir,
                  sub_folder='norm_count/origin', cancer=self.cancer)

        subprocess.call(
            'rm -rf {}'.format('_'.join([store_dir_norm, self.cancer])), shell=True)
        return 'Success'

    def cnv(self):
        '''
        Workflow for downloading copy number variation data from FireBrowse and preprocessing data format.

        Parameters
        ----------
        parental_dir : str
            Path to parental folder that you want to store the whole copy number variation data
        cancer : str
            Cancer name you want to download from FireBrowse, it must be a cancer type included in TCGA project.
        '''
        ## Gene 
        store_dir = '/'.join([self.parental_dir, 'CNV/somatic', 'gene'])
        log = self._fget( data_type='cnv_gene_somatic',store_dir=store_dir)

        if log != 'Success':
            return 'cnv:    '+log

        cnv_gene = self._formatGistic(
            gistic_path='_'.join([store_dir, self.cancer]))
        for name, df in cnv_gene.items():
            df = mergeSampleToPatient(df)
            storeData(df=df, parental_dir=store_dir,
                      sub_folder=name, cancer=self.cancer)
        subprocess.call(
            'rm -rf {}'.format('_'.join([store_dir, self.cancer])), shell=True)

        
        ## Segment
        for lv in ['somatic','all']:
            store_dir = '/'.join([self.parental_dir, 'CNV/'+lv, 'segment'])
            log = self._fget(data_type='cnv_segment_'+lv, store_dir=store_dir)

            if log != 'Success':
                return 'cnv:    '+log

            if not os.path.exists(store_dir):
                    os.makedirs(store_dir)
            subprocess.call(
                'mv {0} {1}'.format('_'.join([store_dir, self.cancer]),
                                    '/'.join([store_dir, self.cancer])
                                    ),
                    shell=True)
        
        return 'Success'

    def rppa(self):
        '''
        Workflow for downloading RPPA data from FireBrowse and preprocessing data format.

        Parameters
        ----------
        parental_dir : str
            Path to parental folder that you want to store the whole RPPA data
        cancer : str
            Cancer name you want to download from FireBrowse, it must be a cancer type included in TCGA project.
        '''
        store_dir = '/'.join([self.parental_dir, 'RPPA'])
       
        log=self._fget(data_type='rppa',store_dir=store_dir)

        if log != 'Success':
            return 'RPPA:    '+log

        rppa = pd.read_table(
            '_'.join([store_dir,self.cancer]), index_col=0)

        rppa = rmEntrez(rppa)
        rppa = mergeSampleToPatient(rppa)

        storeData(df=rppa, parental_dir=store_dir,
                  sub_folder='', cancer=self.cancer)

        subprocess.call(
            'rm -rf {}'.format('_'.join([store_dir, self.cancer])), shell=True)

        return 'Success'

    def snv(self):
        print('''
            Please use MC3 downloader to fetch the SNV result for all cancer in TCGA, 
            which is more robust.
        ''')


class GdcDnloader(GdcApi, Workflow):
    __slot__ = ['type_available', 'base_url']
    def __init__(self, base_url="https://gdc.xenahubs.net/download/",**kwargs):
        Workflow.__init__(self,**kwargs)
        GdcApi.__init__(self, cancer=self.cancer,parental_dir=self.parental_dir)
        # super(GdcDnloader, self).__init__(data_endpt="https://api.gdc.cancer.gov/data",files_endpt="https://api.gdc.cancer.gov/files",**kwargs)
        # data-release-80
        self.base_url = base_url
        self.type_available = {
                    'RNASeq': ['fpkm','count','fpkm_uq'],
                    'SNV': ['muse',"mutect2","VarScan2","SomaticSnipe"],
                    'cnv': ['somatic','all']
                }

    def _fget(self, data_type, store_dir):
        '''Download level 3 data from Xenas
        
        Parameters
        ----------
        data_type : str
            Data type to be downloaded
        store_dir : str
            Path to store the data
        
        Raises
        ------
        KeyError
            If cannot fetching the files
        
        Returns
        -------
        str
            Tell if the downloading is successful or not
        '''

        data_type_dict = {
                'fpkm': "htseq_fpkm",
                'count':"htseq_counts",
                'fpkm_uq': "htseq_fpkm-uq",
                'muse': "muse_snv",
                "mutect2": "mutect2_snv",
                "VarScan2": "varscan2_snv",
                "SomaticSnipe":"somaticsniper_snv",
                }

        if not data_type in data_type_dict.keys():
            raise KeyError("""
            {0} is not a valid data type, only accept following input: {1}
            """.format(data_type, ','.join(data_type_dict.keys())))
        # https: // gdc.xenahubs.net/download/TCGA-CHOL/Xena_Matrices/TCGA-CHOL.htseq_fpkm.tsv.gz
        subpath = 'TCGA-{cancer}/Xena_Matrices/TCGA-{cancer}.{data_type}.tsv.gz'
        url = "/".join([self.base_url, subpath])

        url = url.format(**dict(
                cancer=self.cancer,
                data_type=data_type_dict[data_type]
                )
            )
        cmd = """
        set -x
        [[ -d {store_dir} ]] || mkdir -p {store_dir}
        wget -q -O {store_dir}/{cancer}.gz {url}
        gunzip {store_dir}/{cancer}.gz
        """.format(**dict(
                store_dir=store_dir,
                cancer=self.cancer,
                url=url,
            )
        )
        try:
            subprocess.call(cmd, shell=True)
            log = 'Success'
        except subprocess.CalledProcessError as e:
                log = e

        return log

    def rnaseq(self):
        store_parental = '/'.join([self.parental_dir, 'RNASeq'])
        for name in self.type_available['RNASeq']:
            store_dir = '/'.join([store_parental, name])
            log = self._fget(data_type=name, store_dir=store_dir)

            if log != 'Success':
                return name+':    '+log

            df = pd.read_table('/'.join([store_dir,self.cancer]),index_col=0)
            df = np.exp2(df) - 1  # since all matrix download from xenas have been log transformed
            df = mergeSampleToPatient(df)
            df = mapEm2Gene(df)
            
            if name == 'fpkm':
                tpm = tpmToFpkm(df, reverse=True)
                for raw_name,raw_df in {'tpm':tpm,'fpkm':df}.items():
                    log_df = np.log2(1 + raw_df)
                    tumor_zscore = calTNzcore(log_df, pair_TN=False)
                    storeData(df=tumor_zscore, parental_dir=store_parental,
                            sub_folder=raw_name+'/zscore_tumor/', cancer=self.cancer)
                    try:
                        paired_zscore = calTNzcore(log_df, pair_TN=True)
                        storeData(df=paired_zscore, parental_dir=store_parental,
                                sub_folder=raw_name+'/zscore_paired/', cancer=self.cancer)
                    except ValueError:
                        pass

                    storeData(df=df, parental_dir=store_parental,
                              sub_folder=raw_name+'/origin', cancer=self.cancer)

            else:
                if name == 'count':
                    df = df.round(0)
                storeData(df=df, parental_dir=store_parental,
                          sub_folder=name+'/origin', cancer=self.cancer)

            subprocess.call(
                'rm -rf {}'.format('/'.join([store_dir, self.cancer])), shell=True)
                
    def snv(self):
        store_parental = '/'.join([self.parental_dir, 'SNV'])
        for name in self.type_available['SNV']:
            store_dir = '/'.join([store_parental, name])

            log = self._fget(data_type=name, store_dir=store_dir)

            if log != 'Success':
                return name+':    '+log


    def cnv(self):
        store_parental = '/'.join([self.parental_dir, 'CNV'])
        
        # meta data
        ## map uuid to barcode
        meta, errors = self.getTable(data_type='aliquot')
        if errors != None:
            return errors

        meta = meta.dropna(
            axis=0).set_index('bcr_aliquot_uuid')
        meta.index = meta.index.map(lambda x: x.lower())
        meta = meta['bcr_sample_barcode'].to_dict()

        # focal data
        df,errors = self.getTable(data_type='gistic')
        if errors == None:
            df = df.set_index('Gene Symbol').drop(['Gene ID', 'Cytoband'],axis=1)
            df.columns = df.columns.map(meta)
            df = mergeSampleToPatient(df)
            df = mapEm2Gene(df)
            storeData(df=df, parental_dir=store_parental,
                    sub_folder='somatic/gene/focal', cancer=self.cancer)

        # Segment data
        ## somatic 
        df, errors = self.getTable(data_type='cnv_segment_somatic', by_name=False)
        if errors == None:
            df['GDC_Aliquot'] = df['GDC_Aliquot'].map(meta)
            storeData(df=df, parental_dir=store_parental,
                    sub_folder='somatic/segment', cancer=self.cancer,index=False)

        # all 
        df, errors = self.getTable(data_type='cnv_segment_all', by_name=False)
        if errors == None:
            df['GDC_Aliquot'] = df['GDC_Aliquot'].map(meta)
            storeData(df=df, parental_dir=store_parental,
                    sub_folder='all/segment', cancer=self.cancer, index=False)

        return 'Success'
       

    
    def rppa(self):
        print('RPPA data for hg38 is not available.')

# class mc3Dnloader(object):

#     def __init__(self, store_dir, ref='hg19'):
#         self.ref = ref
#         self.store_dir = store_dir

#     def download(self):
#         cmd = ''' 
#         set - x
#         [-d {store_dir}/] | | mkdir - p {store_dir}/
#         wget -v -O {store_dir}/mc3.gz {url} 
#         tar -xvvf {store_dir}/mc3.gz -C {store_dir}/ --strip-components=1
#         rm {store_dir}/mc3.gz
#         '''

#         if self.ref == 'hg19':
#             url = 'https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc'
#         else:
#             url = ''

#         cmd = cmd.format(**dict(
#             store_dir=self.store_dir,
#             url=url,
#         )
#         )
#         subprocess.call()
        
