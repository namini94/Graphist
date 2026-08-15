import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import stlearn as st
from PIL import Image
import numpy as np
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


data_root = Path('/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/outs')
annot_file = '/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/Nami-Annotations.csv'
morph_file = '/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/Morph-ResNet50.csv'

anno_df = pd.read_csv(annot_file, index_col=0)
morph_df = pd.read_csv(morph_file, index_col=0)



adata = sc.read_visium(data_root)
adata.var_names_make_unique()

anno_df['V1'] = anno_df['V1'].replace({'fat': 'Fat'})
anno_df['V1'] = anno_df['V1'].replace({'necrosis': 'Necrosis'})
anno_df['V1'] = anno_df['V1'].replace({'tumor': 'Invasive Carcinoma'})
anno_df['V1'] = anno_df['V1'].replace({'immune': 'Immune Cell'})
anno_df['V1'] = anno_df['V1'].replace({'fibrous': 'Fibrous Tissue'})
anno_df['V1'] = anno_df['V1'].replace({'Not Annotated': 'Not Annotated'})


adata.obs = adata.obs.join(anno_df)
adata.obs['annotations'] = adata.obs.V1

adata.obsm['morph'] = morph_df

#print(adata)
#print(adata.obsm['morph'])

temp_adata = sc.AnnData(X=adata.obsm['morph'])
sc.pp.neighbors(temp_adata, n_neighbors=15,random_state=0)
sc.tl.louvain(temp_adata, resolution=1, random_state=0)
adata.obs['louvain'] = temp_adata.obs['louvain']
sc.set_figure_params(dpi=80, figsize=(4, 4))
sc.pl.spatial(adata, color='louvain', spot_size=275, show=True, frameon=False, title='Image Feature Clustering (ResNet50)')
    
def plot_clustered_obsm_heatmap(adata, obsm_key, obs_key, figsize=(30, 20), cmap="viridis"):
    # Extract the data from obsm
    data = adata.obsm[obsm_key]
    print(f"Shape of data: {data.shape}")
    print(f"Type of data: {type(data)}")
    
    # Get annotation
    annotation = adata.obs[obs_key]
    print(f"Shape of annotation: {annotation.shape}")
    
    # Handle MultiIndex
    if isinstance(annotation.index, pd.MultiIndex):
        print("MultiIndex detected. Using the first level for coloring.")
        annotation = annotation.index.get_level_values(0)
    else:
        annotation = annotation.astype(str)
    
    # Create color map for annotations
    unique_annotations = sorted(annotation.unique())
    print(f"Number of unique annotations: {len(unique_annotations)}")
    color_map = dict(zip(unique_annotations, range(len(unique_annotations))))
    
    # Create annotation color list
    row_colors = pd.Series(annotation).map(color_map)
    
    # Perform clustering
    print("Performing clustering...")
    linkage = hierarchy.linkage(pdist(data.values), method='average')
    
    # Set up the matplotlib figure
    fig = plt.figure(figsize=figsize)
    
    # Create a gridspec for the dendrogram, heatmap, and color bar
    gs = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[0.05, 1], height_ratios=[0.2, 1],
                          left=0.05, right=0.9, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Plot the dendrogram
    ax_dendrogram = fig.add_subplot(gs[0, 1])
    hierarchy.dendrogram(linkage, ax=ax_dendrogram, labels=None, leaf_rotation=90)
    ax_dendrogram.set_xticks([])
    ax_dendrogram.set_yticks([])
    
    # Plot the heatmap
    ax_heatmap = fig.add_subplot(gs[1, 1])
    im = ax_heatmap.imshow(data.values[hierarchy.leaves_list(linkage)], aspect='auto', cmap=cmap)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    
    # Add colorbar for heatmap
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=f'{obsm_key} values')
    
    # Plot the annotation colors
    ax_colors = fig.add_subplot(gs[1, 0])
    color_array = row_colors.iloc[hierarchy.leaves_list(linkage)].to_numpy().reshape(-1, 1)
    ax_colors.imshow(color_array, aspect='auto', cmap=ListedColormap(sns.color_palette("husl", len(unique_annotations))))
    ax_colors.set_xticks([])
    ax_colors.set_yticks([])
    
    # Add a legend for annotations
    handles = [plt.Rectangle((0,0),1,1, color=sns.color_palette("husl", len(unique_annotations))[i]) for i in range(len(unique_annotations))]
    fig.legend(handles, unique_annotations, title=obs_key, 
               bbox_to_anchor=(1, 0.9), loc='upper left')
    
    plt.suptitle(f"Clustered Heatmap of {obsm_key} values\nAnnotated with {obs_key}", fontsize=16)
    plt.savefig('clustered_obsm_heatmap_annotated.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'clustered_obsm_heatmap_annotated.png'")
    plt.close()  # Close the figure to free up memory



plot_clustered_obsm_heatmap(adata,obsm_key='morph',obs_key='annotations')