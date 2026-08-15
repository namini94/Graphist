import scanpy as sc
import pandas as pd 
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


data_root = Path('/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/outs')
annot_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/Nami-Annotations.csv'

anno_df = pd.read_csv(annot_file, index_col=0)



adata = sc.read_visium(data_root)
adata.var_names_make_unique()

anno_df['V1'] = anno_df['V1'].replace({'fat': 'Fat'})
anno_df['V1'] = anno_df['V1'].replace({'necrosis': 'Necrosis'})
anno_df['V1'] = anno_df['V1'].replace({'tumor': 'Invasive Carcinoma'})
anno_df['V1'] = anno_df['V1'].replace({'immune': 'Immune Cell'})
anno_df['V1'] = anno_df['V1'].replace({'fibrous': 'Fibrous Tissue'})
anno_df['V1'] = anno_df['V1'].replace({'Not Annotated': 'na'})


adata.obs = adata.obs.join(anno_df)
adata.obs['annotations'] = adata.obs.V1




adata.layers['count'] = adata.X.toarray()
sc.pp.filter_genes(adata, min_cells=10)
print('After flitering: ', adata.shape)


# Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.set_figure_params(dpi=80, figsize=(4, 4))
sc.pl.spatial(adata, color='annotations', spot_size=275,frameon=False, alpha_img=1,
              palette = ['lightskyblue','darkorchid' , 'gold', 'seagreen','coral', 'gray']
, title = "Breast Cancer")

#'saddlebrown'
print(adata)
print(adata.obsm['spatial'])
#pd.DataFrame(adata.to_df()).to_csv("/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/processed/Normalized_BRCA_PACSI.csv",index=True)



selected = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/selected_al3e-2.csv',
    index_col=0, header=0
)

#selected = pd.read_csv('/Users/naminiyakan/Documents/Multi-Agent-Design/BioMaster/outputs/BRCA_output/data.csv',
#    index_col=0, header=0
#)


selected['x'] = selected['x'].replace({0: 'Background'})
selected['x'] = selected['x'].replace({1: 'Graphist (+)'})
selected['x'] = selected['x'].replace({2: 'Graphist (-)'})
adata.obs = adata.obs.join(selected)
adata.obs['selection'] = adata.obs.x

n_tot = 0
n_NA = 0
n_tumor = 0 
n_necrosis = 0
n_fibr = 0
n_immune = 0
n_fat = 0
for i in range(adata.n_obs):
    if((adata.obs['selection'].iloc[i] == 'Graphist (-)')):
        n_tot = n_tot + 1 
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'na')):
        n_NA = n_NA + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'Invasive Carcinoma')):
        n_tumor = n_tumor + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'Necrosis')):
        n_necrosis = n_necrosis + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'Fibrous Tissue')):
        n_fibr = n_fibr + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'Immune Cell')):
        n_immune = n_immune + 1
    if ((adata.obs['selection'].iloc[i] == 'Graphist (-)') & (adata.obs['annotations'].iloc[i]== 'Fat')):
        n_fat = n_fat + 1
    
    
    
print('n_tot',n_tot)
print('n_NA:', n_NA) 
print('n_tumor:', n_tumor)  
print('n_necrosis:', n_necrosis)  
print('n_fibr', n_fibr)  
print('n_immune', n_immune)  
print('n_fat', n_fat)  


sc.pl.spatial(adata, color='selection', spot_size=275,frameon=False,title='Phenotype: Survival' ,alpha_img=1, palette=["gray", "blue", "red"])

# Select marker genes
#show_gene=["FTH1","ARPC1B","PPDPF","KRT7","C1R","DCN","C1S"]
#show_gene = ["FOS", "FN1","EGR1","SERPINA3"]
show_gene = ["CXCL12", "CXCR4","EDN1", "EDNRA"]

sc.pl.spatial(adata, color=show_gene, spot_size=275,frameon=False)
#sc.pl.embedding(adata, basis="coord", color=show_gene, s=275, show=True, frameon=False)
#category_order = ['Graphist (+)','Graphist (-)', 'Background']


sc.pl.stacked_violin(adata, show_gene, groupby='selection', dendrogram=False, return_fig=False)
#vp.add_totals().style(ylim=(0,5)).show()

##### Immune cells investigation
adata.obs['immune_groups'] = 'Background'
adata.obs.loc[(adata.obs['annotations']== 'Invasive Carcinoma'), 'immune_groups'] = 'Invasive Carcinoma'
adata.obs.loc[(adata.obs['annotations']== 'Fibrous Tissue'), 'immune_groups'] = 'Fibrous Tissue'
adata.obs.loc[(adata.obs['selection'] == 'Background') & (adata.obs['annotations']== 'Immune Cell'), 'immune_groups'] = 'Immune(Background)'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (-)') & ((adata.obs['annotations']!= 'Invasive Carcinoma') & (adata.obs['annotations']!= 'Fibrous Tissue')), 'immune_groups'] = 'Graphist(-)'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (+)') & ((adata.obs['annotations']!= 'Invasive Carcinoma') & (adata.obs['annotations']!= 'Fibrous Tissue')), 'immune_groups'] = 'Graphist(+)'
custom_colors = {'Graphist(-)': 'red', 'Graphist(+)': 'blue', 'Immune(Background)': 'gold', 'Invasive Carcinoma': 'seagreen','Background': 'gray','Fibrous Tissue': 'darkorchid'}

sc.pl.spatial(adata, color='immune_groups', spot_size=275, show=True, frameon=False, palette=custom_colors)

pd.DataFrame(adata.obs['immune_groups']).to_csv("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/Comparison-immune(-)-immune(+).csv",index=True)


##### NEW Immune cells investigation (pathway DE)
adata.obs['immune_groups_DE'] = 'Background'
#adata.obs.loc[(adata.obs['annotations']== 'Invasive Carcinoma'), 'immune_groups_DE'] = 'Invasive Carcinoma'
#adata.obs.loc[(adata.obs['annotations']== 'Fibrous Tissue'), 'immune_groups_DE'] = 'Fibrous Tissue'
#adata.obs.loc[(adata.obs['selection'] == 'Background') & (adata.obs['annotations']== 'Immune Cell'), 'immune_groups_DE'] = 'Immune(Background)'
adata.obs.loc[(adata.obs['selection'] == 'Background') & (adata.obs['annotations']== 'Immune Cell'), 'immune_groups_DE'] = 'Immune(Background)'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (-)') & (adata.obs['annotations']== 'Immune Cell'), 'immune_groups_DE'] = 'Immune(-)'
adata.obs.loc[(adata.obs['selection'] == 'Graphist (+)') & (adata.obs['annotations']== 'Immune Cell'), 'immune_groups_DE'] = 'Immune(+)'
#adata.obs.loc[(adata.obs['selection'] == 'Graphist (-)') & ((adata.obs['annotations']!= 'Invasive Carcinoma') & (adata.obs['annotations']!= 'Fibrous Tissue')), 'immune_groups_DE'] = 'Graphist(-)'
#adata.obs.loc[(adata.obs['selection'] == 'Graphist (+)') & ((adata.obs['annotations']!= 'Invasive Carcinoma') & (adata.obs['annotations']!= 'Fibrous Tissue')), 'immune_groups_DE'] = 'Graphist(+)'
custom_colors = {'Immune(-)': 'red', 'Immune(+)': 'blue', 'Immune(Background)': 'gold','Background': 'gray'}

sc.pl.spatial(adata, color='immune_groups_DE', spot_size=275, show=True, frameon=False, palette=custom_colors)

pd.DataFrame(adata.obs['immune_groups_DE']).to_csv("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/DE_Comparison-immune(-)-immune(+).csv",index=True)



adata_temp = adata[adata.obs['annotations']=='Immune Cell']

adata_temp.obs['immune_gp']= 'Immune(Background)'
adata_temp.obs.loc[(adata_temp.obs['selection']== 'Graphist (-)'), 'immune_gp'] = 'Immune(Graphist(-))'
adata_temp.obs.loc[(adata_temp.obs['selection']== 'Graphist (+)'), 'immune_gp'] = 'Immune(Graphist(+))'


show_gene2 = ['S100A11','CCL19','FTH1','SLC11A1','TRBC2','FN1','C3']
#show_gene2 = ["CXCL12", "CXCR4","EDN1", "EDNRA"]
sc.pl.stacked_violin(adata_temp, show_gene2, groupby='immune_gp', dendrogram=False, return_fig=False)
sc.pl.spatial(adata_temp, color=show_gene2, spot_size=275,frameon=False)

# Fixing sparse matrix access
cxcl12_expr = adata[:, 'CXCL12'].X.toarray().flatten()
cxcr4_expr  = adata[:, 'CXCR4'].X.toarray().flatten()

# Create interaction score
adata.obs['CXCL12_CXCR4_score'] = cxcl12_expr * cxcr4_expr

cxcl12_expr = adata_temp[:, 'FN1'].X.toarray().flatten()
cxcr4_expr  = adata_temp[:, 'ITGA5'].X.toarray().flatten()

adata_temp.obs['FN1_ITGA5_score'] = cxcl12_expr * cxcr4_expr

# Plot
sc.pl.spatial(adata_temp, color='FN1_ITGA5_score', show=True,cmap='plasma', spot_size=275, frameon=False, title='FN1 × ITGA5 Co-expression')


pairs = [
    ("CXCL12", "CXCR4"),
    ("CCL19","CCR7"),
    ("MIF", "CD74"),
    ("CCL5", "CCR5"),
    ("IL7", "IL7R") , # for prostaglandin signaling (if relevant)
    ("FN1","ITGA5"),
    ("FN1","ITGB1")
]
sc.pl.dotplot(adata, [
                      'FN1', 'ITGA5', 'ITGB1', 'ITGAV', 'ITGB3'], 
              groupby='immune_groups_DE', 
              standard_scale='var')

for ligand, receptor in pairs:
    try:
        ligand_expr = adata[:, ligand].X.toarray().flatten()
        receptor_expr = adata[:, receptor].X.toarray().flatten()
        score = ligand_expr * receptor_expr
        colname = f"{ligand}_{receptor}_score"
        adata.obs[colname] = score
        
        sc.pl.spatial(
            adata,
            color=colname,
            cmap='plasma',
            spot_size=275,
            frameon=False,
            title=f"{ligand} × {receptor} Co-expression"
        )
    except KeyError:
        print(f"Skipping {ligand}-{receptor}: gene not found in data")
        
        

adata.obs['CXCL12_CXCR4_score'] = adata[:, 'CXCL12'].X.toarray().flatten() * adata[:, 'CXCR4'].X.toarray().flatten()

sc.pl.violin(
    adata,
    keys='CXCL12_CXCR4_score',
    groupby='immune_groups_DE',
    stripplot=True,
    jitter=0.4,
    rotation=45,
    show=True
)
plt.title("CXCL12 × CXCR4 Score by Tissue Type")


fig, ax = plt.subplots(figsize=(6,6))
sc.pl.spatial(
    adata,
    color='CXCL12_CXCR4_score',
    groups=None,
    spot_size=275,
    frameon=False,
    show=False,
    ax=ax,
    cmap='plasma'
)

# Overlay annotations
adata.obs['annot_overlay'] = adata.obs['annotations'].astype(str)
for category in adata.obs['annot_overlay'].unique():
    subset = adata[adata.obs['annot_overlay'] == category]
    plt.scatter(
        subset.obsm['spatial'][:, 0],
        subset.obsm['spatial'][:, 1],
        s=10,
        label=category,
        alpha=0.3
    )
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("CXCL12 × CXCR4 Score with Annotation Overlay")
plt.axis('off')
plt.show()