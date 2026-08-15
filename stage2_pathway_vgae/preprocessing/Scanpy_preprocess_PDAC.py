import scanpy as sc
import pandas as pd 
from pathlib import Path


counts_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/PDAC-A-ST1-filtered.txt'
meta_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/STAGE-PDAC-A-ST1-meta.txt'
pos_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/Pos-A-ST1.csv'
image_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/PDAC-A-ST1-HE.jpg'




counts = pd.read_csv(counts_file, sep='\t', index_col=0)
counts = counts.T
meta_df = pd.read_csv(meta_file, sep='\t', index_col=0)
print(counts.shape, meta_df.shape)

# Coordinates (coor_x, coor_y) and label (human_anno_region)
adata = sc.AnnData(counts)
coor_df = meta_df.loc[adata.obs_names, ["coor_x", "coor_y"]]
adata.obsm["coord"] = coor_df.to_numpy()
spatial = meta_df.loc[adata.obs_names, ["coor_x", "coor_y"]]
adata.obsm["spatial"] = spatial.to_numpy()
adata.obs[meta_df.columns] = meta_df.loc[adata.obs_names, meta_df.columns]


#adata.layers['count'] = adata.X.toarray()
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e6)
#sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
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
pd.DataFrame(count).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Processed-A/Normalized_5k_PDAC.csv",index=True)
#pd.DataFrame(portal_adata.obs).to_csv("/Users/naminiyakan/Documents/VEGA_code/TCDD/metadata_portal.csv",index=True)
pd.DataFrame(adata.var['highly_variable']).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Processed-A/5k_HVG.csv",index=True)
