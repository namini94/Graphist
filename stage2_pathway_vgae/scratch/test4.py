import pandas as pd 
exp_paths = '/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/processed/Normalized_5k_BRCA_PACSI.csv'
X = pd.read_csv(exp_paths,header=0,index_col=0)
print(X)