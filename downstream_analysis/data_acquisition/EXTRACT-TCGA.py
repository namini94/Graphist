import requests
import json
import pandas as pd
import numpy as np


def extract_tcga_data(project="TCGA-BRCA", data_category="Transcriptome Profiling", data_type="Gene Expression Quantification"):
    files_endpt = "https://api.gdc.cancer.gov/files"

    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [project]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [data_category]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [data_type]
                }
            }
        ]
    }

    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.submitter_id,data_category,data_type",
        "format": "JSON",
        "size": "100"
    }

    try:
        response = requests.get(files_endpt, params=params)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        
        data = response.json()
        print("API Response:")
        print(json.dumps(data, indent=2))  # Print the full response for debugging
        
        if 'data' in data and 'hits' in data['data']:
            hits = data['data']['hits']
            if hits:
                return pd.DataFrame(hits)
            else:
                print("No data found matching the criteria.")
                return None
        else:
            print("Unexpected response structure.")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error making request: {e}")
        return None
    
def download_file_content(file_id):
    data_endpt = f"https://api.gdc.cancer.gov/data/{file_id}"
    response = requests.get(data_endpt)

metadata = extract_tcga_data()
print(metadata)

if metadata is not None and not metadata.empty:
    # Download and process the first file as an example
    first_file_id = metadata.iloc[0]['file_id']
    gene_expression_data = download_file_content(first_file_id)

    if gene_expression_data is not None:
        gene_expression_data.to_csv('/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/outs/tcga_gene_expression.csv')
        