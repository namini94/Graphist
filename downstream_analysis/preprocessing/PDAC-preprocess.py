import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import stlearn as st
from PIL import Image
import numpy as np
import seaborn as sns

counts_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/PDAC-A-ST1-filtered.txt'
meta_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/STAGE-PDAC-A-ST1-meta.txt'
pos_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/Pos-A-ST1.csv'
image_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/PDAC-A-ST1-HE.jpg'

#data = st.ReadOldST(count_matrix_file=counts_file,
#                    spatial_file=pos_file,
#                    image_file=image_file)

#print(data)


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

sc.set_figure_params(dpi=80, figsize=(4, 4))
sc.pl.embedding(adata, basis="coord", color="human_anno_region", title='Pancreatic Cancer', s=120, show=True, frameon=False)

sc.pp.filter_genes(adata, min_cells=10)
print('After flitering: ', adata.shape)

# Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


#sc.tl.pca(adata)
#sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
# Compute UMAP
#sc.tl.umap(adata)
#sc.pl.umap(adata,s=20, color="human_anno_region")

selected = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-PDAC/DDX60L_selected_al5e-1.csv',
    index_col=0, header=0
)


#selected = pd.read_csv('/Users/naminiyakan/Documents/Multi-Agent-Design/BioMaster/outputs/PAAD_output/data.csv',
#    index_col=0, header=0
#)


selected['x'] = selected['x'].replace({0: 'Background'})
selected['x'] = selected['x'].replace({1: 'Graphist (+)'})
selected['x'] = selected['x'].replace({2: 'Graphist (-)'})
adata.obs = adata.obs.join(selected)
adata.obs['selection'] = adata.obs.x


n_cancer = 0
n_stroma = 0 
n_pancreatic = 0
n_duct = 0
n_tot = 0
for i in range(adata.n_obs):
    if((adata.obs['selection'].iloc[i] == 'Graphist (-)')):
        n_tot = n_tot + 1 
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['human_anno_region'].iloc[i]== 'Cancer region')):
        n_cancer = n_cancer + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['human_anno_region'].iloc[i]== 'Stroma')):
        n_stroma = n_stroma + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['human_anno_region'].iloc[i]== 'Pancreatic tissue')):
        n_pancreatic = n_pancreatic + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['human_anno_region'].iloc[i]== 'Duct epithelium')):
        n_duct = n_duct + 1
    
    
print('n_tot',n_tot)
print('n_cancer:', n_cancer) 
print('n_stroma:', n_stroma)  
print('n_duct:', n_duct)  
print('n_pancreatic', n_pancreatic)  
    




sc.pl.embedding(adata, basis="coord", color="selection", title='Phenotype: DDX60L Knockdown', s=120, show=True, frameon=False,
                palette=["gray", "blue", "red"])

# Select marker genes
#show_gene = ["S100A6","KRT19","LAMC2","REG3A","CTRB2","MUC5B"]
#show_gene=["LAMC2","KRT17","MROH6"]
#show_gene = ["KRT19","S100A6","GAPDH","KRT17","PIGR","C3","MUC5B"]
show_gene = ["S100A6","TMSB4X","KRT19"]

sc.pl.embedding(adata, basis="coord", color=show_gene, s=120, show=True, frameon=False)
#category_order = ['Graphist (+)','Graphist (-)', 'Background']


sc.pl.stacked_violin(adata, show_gene, groupby='selection', dendrogram=False, return_fig=False)
#vp.add_totals().style(ylim=(0,5)).show()


palette = sns.color_palette("muted",4)
sc.pl.stacked_violin(adata, show_gene, groupby='human_anno_region', dendrogram=True,row_palette=palette )


#####Stroma Cells investigation
# Another example: color cells based on multiple conditions
adata.obs['stroma_groups'] = 'Back'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (-)') & (adata.obs['human_anno_region']== 'Stroma'), 'stroma_groups'] = 'Stroma(-)'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (+)') & (adata.obs['human_anno_region']== 'Stroma'), 'stroma_groups'] = 'Stroma(+)'
custom_colors = {'Stroma(-)': 'red', 'Stroma(+)': 'blue', 'Back': 'gray'}

sc.pl.embedding(adata, basis="coord", color='stroma_groups', s=120, show=True, frameon=False, palette=custom_colors)


print(adata)
pd.DataFrame(adata.to_df()).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/Processed-A/Normalized_PDAC.csv",index=True)



# Load and process MYEOV data
MYEOV_selected = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-PDAC/MYEOV_selected_al5e-1.csv',
    index_col=0, header=0
)
MYEOV_selected['MYEOV_x'] = MYEOV_selected['x'].replace({
    0: 'Background',
    1: 'Graphist (+)',
    2: 'Graphist (-)'
})
MYEOV_selected = MYEOV_selected.drop('x', axis=1)  # Drop the original 'x' column

# Load and process DDX60L data
DDX60L_selected = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-PDAC/DDX60L_selected_al6e-1.csv',
    index_col=0, header=0
)
DDX60L_selected['DDX60L_x'] = DDX60L_selected['x'].replace({
    0: 'Background',
    1: 'Graphist (+)',
    2: 'Graphist (-)'
})
DDX60L_selected = DDX60L_selected.drop('x', axis=1)  # Drop the original 'x' column

# Join the data with adata.obs
adata.obs = adata.obs.join(MYEOV_selected)
adata.obs = adata.obs.join(DDX60L_selected)

# Create selection columns
adata.obs['MYEOV_selection'] = adata.obs.MYEOV_x
adata.obs['DDX60L_selection'] = adata.obs.DDX60L_x

# Create Phenotypes column
adata.obs['Phenotypes'] = 'Rest'
adata.obs.loc[(adata.obs['DDX60L_selection'] == 'Graphist (-)') & (adata.obs['MYEOV_selection'] == 'Graphist (-)'), 'Phenotypes'] = 'Common'
adata.obs.loc[(adata.obs['DDX60L_selection'] == 'Graphist (-)') & (adata.obs['MYEOV_selection'] != 'Graphist (-)'), 'Phenotypes'] = 'DDX60L Only'
adata.obs.loc[(adata.obs['DDX60L_selection'] != 'Graphist (-)') & (adata.obs['MYEOV_selection'] == 'Graphist (-)'), 'Phenotypes'] = 'MYEOV Only'

pd.DataFrame(adata.obs['Phenotypes']).to_csv("/Users/naminiyakan/Documents/BulkToST/Res-PDAC/Comparison-DDX60L-MYEOC.csv",index=True)

# Define custom colors
custom_colors = {'DDX60L Only': 'purple', 'MYEOV Only': 'darkorange', 'Rest': 'gray', 'Common': 'teal'}

# Plot the embedding
sc.pl.embedding(adata, basis="coord", color='Phenotypes', s=120, show=True, frameon=False, palette=custom_colors)


show_gene = ["TM4SF1","S100A4","S100A6"]



sc.pl.stacked_violin(adata, show_gene, groupby='Phenotypes', dendrogram=True, return_fig=False)


n_cancer = 0
n_stroma = 0 
n_pancreatic = 0
n_duct = 0
n_tot = 0
for i in range(adata.n_obs):
    if ((adata.obs['Phenotypes'].iloc[i] == 'Back') & (adata.obs['human_anno_region'].iloc[i]== 'Cancer region')):
        n_cancer = n_cancer + 1
    if ((adata.obs['Phenotypes'].iloc[i] == 'Back') & (adata.obs['human_anno_region'].iloc[i]== 'Stroma')):
        n_stroma = n_stroma + 1
    if ((adata.obs['Phenotypes'].iloc[i] == 'Back') & (adata.obs['human_anno_region'].iloc[i]== 'Pancreatic tissue')):
        n_pancreatic = n_pancreatic + 1
    if ((adata.obs['Phenotypes'].iloc[i] == 'Back') & (adata.obs['human_anno_region'].iloc[i]== 'Duct epithelium')):
        n_duct = n_duct + 1
    
    

print('n_cancer:', n_cancer) 
print('n_stroma:', n_stroma)  
print('n_duct:', n_duct)  
print('n_pancreatic', n_pancreatic)  
