import scanpy as sc
import pandas as pd 
from pathlib import Path

data_root = Path('/Users/naminiyakan/Documents/BulkToST/Dataset/Maynard/151673/outs')

adata = sc.read_visium(data_root)
adata.var_names_make_unique()



adata.layers['count'] = adata.X.toarray()
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer='count', n_top_genes=5000)
adata = adata[:, adata.var['highly_variable'] == True]
sc.pp.scale(adata)

# sklearn PCA is used because PCA in scanpy is not stable. 
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable. 
adata_X = PCA(n_components=200, random_state=42).fit_transform(adata.X)
adata.obsm['X_pca'] = adata_X

#print(adata)
#print(adata.obs)
#print(adata.obsm['spatial'])
print(adata.X)

GT = pd.read_csv(
    '/Users/naminiyakan/Documents/BulkToST/Dataset/Maynard/151673/processed/GT-labels.csv',
    header=0,
    index_col=0
)
#adata.obs = adata.obs.join(GT)
adata.obs['Ground-truth'] = pd.Categorical(GT.iloc[:,0])
print(adata.obs['Ground-truth'])


#idx_g1 = dosage.index.values[adata.iloc[:,27] == g]

# image_key = None (Remove Background Image)
sc.pl.spatial(adata, color='Ground-truth',frameon=False, alpha_img=0.9)

#Remove Immune Cells
#cell_types_of_int = ["Hepatocytes - central", "Hepatocytes - portal", "Cholangiocytes", "Stellate Cells", "Portal Fibroblasts", "Endothelial Cells"]
#adata = adata[adata.obs['celltype'].isin(cell_types_of_int)]



#cell = 'Hepatocytes - portal'
#portal_adata = adata[(adata.obs["celltype"] == cell)]
#train_adata, test_adata = prepare_cont_data(adata, "celltype", "dose", "Dose", cell, 0, normalized=True)
#print(portal_adata)
#print(portal_adata.obs)


count = adata.to_df()
print(count)
pd.DataFrame(count).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/Maynard/151673/processed/Normalized_5k_151673.csv",index=True)
#pd.DataFrame(portal_adata.obs).to_csv("/Users/naminiyakan/Documents/VEGA_code/TCDD/metadata_portal.csv",index=True)
pd.DataFrame(adata.var['highly_variable']).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/Maynard/151673/processed/5k_HVG.csv",index=True)
