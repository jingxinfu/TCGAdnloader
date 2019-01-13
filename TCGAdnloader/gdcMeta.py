#!/usr/bin/env python
import requests
import json
import re,os
import pandas as pd 
import io

class GdcApi(object):
    __slot__ = ["files_endpt", "data_endpt", "cancer"]

    def __init__(self, cancer,data_endpt="https://api.gdc.cancer.gov/data",files_endpt="https://api.gdc.cancer.gov/files", **kwargs):
        self.files_endpt = files_endpt
        self.data_endpt = data_endpt
        self.cancer = cancer

    def _contructFilter(self,data_type):
        dtype_dict = {
            'gistic': '{}.focal_score_by_genes.txt'.format(self.cancer.upper()),
            'aliquot': "nationwidechildrens.org_biospecimen_aliquot_{}.txt".format(self.cancer.lower()),
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
        self._params = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "size": "1"
        }
    
    def _fetchFileID(self):
        self.file_uuid_list = []
        response = requests.get(self.files_endpt, params=self._params)
        for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
            self.file_uuid_list.append(file_entry["file_id"])
        
    def get(self, data_type):
        self._contructFilter(data_type)
        self._fetchFileID()
        params = {"ids": self.file_uuid_list}
        response = requests.post(self.data_endpt, data=json.dumps(
            params), headers={"Content-Type": "application/json"})
        df = pd.read_table(io.StringIO(response.content.decode("utf-8"))).dropna(axis=0)
        print(df.head())
        

metas = GdcApi(cancer='THCA')
metas.get('aliquot')

# This step populates the download list with the file_ids from the previous query


