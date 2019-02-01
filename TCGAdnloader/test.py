import requests
import json
import pandas as pd
import io
cases_endpt = 'https://api.gdc.cancer.gov/cases'

# The 'fields' parameter is passed as a comma-separated string of single names.
fields = [
    "submitter_id",
    # "diagnoses.age_at_diagnosis",
    # "diagnoses.vital_status",
    # "diagnoses.days_to_death",
    # "diagnoses.days_to_last_follow_up",
    # "diagnoses.tumor_stage",
    "demographic.gender",
   
]
filters = {

        "op": "in",
        "content": {
            "field": "cases.project.project_id",
            "value": [
                "TCGA-SKCM"
            ]
        }
    }
fields = ','.join(fields)

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "TSV",
    "size": "3000"
}

response = requests.get(cases_endpt, params=params)
result = pd.read_table(io.StringIO(response.content.decode("utf-8")))
print(result.columns)
print(result.head())
print(result.shape)
