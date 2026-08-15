import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import stlearn as st
from PIL import Image
import numpy as np
import seaborn as sns

morph_file = '/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/Morph-ResNet50.csv'

#data = st.ReadOldST(count_matrix_file=counts_file,
#                    spatial_file=pos_file,
#                    image_file=image_file)

#print(data)


counts = pd.read_csv(morph_file, index_col=0)
#counts = counts.T
print(counts)

adata = sc.AnnData(counts)

annot_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/Nami-Annotations.csv'

anno_df = pd.read_csv(annot_file, index_col=0)


anno_df['V1'] = anno_df['V1'].replace({'fat': 'Fat'})
anno_df['V1'] = anno_df['V1'].replace({'necrosis': 'Necrosis'})
anno_df['V1'] = anno_df['V1'].replace({'tumor': 'Invasive Carcinoma'})
anno_df['V1'] = anno_df['V1'].replace({'immune': 'Immune Cell'})
anno_df['V1'] = anno_df['V1'].replace({'fibrous': 'Fibrous Tissue'})
anno_df['V1'] = anno_df['V1'].replace({'Not Annotated': 'na'})


#adata.obs = adata.obs.join(anno_df)
#adata.obs['annotations'] = adata.obs.V1

# 1. Reset the index of the counts DataFrame
counts_reset = counts.reset_index(drop=True)

# 2. Create a new index from anno_df
new_index = anno_df.index

# 3. Assign this new index to counts_reset
counts_reset.index = new_index

# 4. Create a new AnnData object with the correctly indexed counts
adata_new = sc.AnnData(counts_reset)

# 5. Now the obs DataFrame will have the correct index
print(adata_new.obs.head())

# 6. Join with the annotations
adata_new.obs = adata_new.obs.join(anno_df)

# 7. Assign annotations
adata_new.obs['annotations'] = adata_new.obs['V1'].astype('category')

# 8. Verify the results
print(adata_new.obs['annotations'].value_counts())
print(adata_new.obs['annotations'].isna().sum())

# If you want to update your original adata object:
adata = adata_new


print(adata)
print(adata.obs['annotations'])
#sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=100)
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'annotations')



