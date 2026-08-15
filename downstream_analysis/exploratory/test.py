import scanpy as sc
import pandas as pd 
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches


warnings.filterwarnings('ignore')

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

print(adata)
print(adata.obsm['spatial'])

selected = pd.read_csv('/Users/naminiyakan/Documents/Multi-Agent-Design/BioMaster/outputs/BRCA_output/data.csv',
    index_col=0, header=0
)

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
show_genes = ["FOS", "FN1","EGR1","SERPINA3"]

# ===============================================================================
# GENE-REGION CORRELATION ANALYSIS
# ===============================================================================

def calculate_gene_region_correlations(adata, genes, region_col='annotations'):
    """
    Calculate correlations between gene expression and cancer regions
    """
    # Filter genes that exist in the dataset
    available_genes = [g for g in genes if g in adata.var_names]
    missing_genes = [g for g in genes if g not in adata.var_names]
    
    if missing_genes:
        print(f"Warning: The following genes are not found in the dataset: {missing_genes}")
    
    if not available_genes:
        print("No genes found in the dataset!")
        return None, None
    
    print(f"Analyzing {len(available_genes)} genes: {available_genes}")
    
    # Get expression data for available genes
    gene_expr = pd.DataFrame(
        adata[:, available_genes].X.toarray(), 
        index=adata.obs_names, 
        columns=available_genes
    )
    
    # Get region annotations
    regions = adata.obs[region_col].copy()
    
    # Remove spots with missing annotations
    valid_mask = regions != 'na'
    gene_expr_clean = gene_expr[valid_mask]
    regions_clean = regions[valid_mask]
    
    # Get unique regions
    unique_regions = regions_clean.unique()
    print(f"Analyzing regions: {unique_regions}")
    
    # Calculate correlations for each gene-region pair
    correlation_results = []
    pvalue_results = []
    
    for gene in available_genes:
        gene_correlations = []
        gene_pvalues = []
        
        for region in unique_regions:
            # Create binary indicator for region
            region_indicator = (regions_clean == region).astype(int)
            
            # Calculate Pearson correlation
            corr, p_val = pearsonr(gene_expr_clean[gene], region_indicator)
            
            gene_correlations.append(corr)
            gene_pvalues.append(p_val)
        
        correlation_results.append(gene_correlations)
        pvalue_results.append(gene_pvalues)
    
    # Create DataFrames
    corr_df = pd.DataFrame(
        correlation_results, 
        index=available_genes, 
        columns=unique_regions
    )
    
    pval_df = pd.DataFrame(
        pvalue_results, 
        index=available_genes, 
        columns=unique_regions
    )
    
    return corr_df, pval_df

def plot_correlation_heatmap(corr_df, pval_df, title="Gene-Region Correlations"):
    """
    Plot heatmap of gene-region correlations with significance markers
    """
    # Create significance markers
    sig_markers = np.where(pval_df < 0.001, '***',
                  np.where(pval_df < 0.01, '**',
                  np.where(pval_df < 0.05, '*', '')))
    
    # Create the heatmap
    plt.figure(figsize=(12, 8))
    
    # Plot heatmap
    sns.heatmap(corr_df, 
                annot=True, 
                cmap='RdBu_r', 
                center=0,
                fmt='.3f',
                cbar_kws={'label': 'Pearson Correlation'},
                linewidths=0.5)
    
    # Add significance markers
    for i in range(len(corr_df.index)):
        for j in range(len(corr_df.columns)):
            if sig_markers[i, j]:
                plt.text(j + 0.7, i + 0.3, sig_markers[i, j], 
                        fontsize=12, fontweight='bold', color='black')
    
    plt.title(f'{title}\n* p<0.05, ** p<0.01, *** p<0.001', fontsize=14, fontweight='bold')
    plt.xlabel('Cancer Regions', fontsize=12)
    plt.ylabel('Genes', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()
    
    return plt.gcf()

def analyze_expression_patterns(adata, genes, region_col='annotations'):
    """
    Analyze expression patterns of genes across different regions
    """
    available_genes = [g for g in genes if g in adata.var_names]
    
    # Create expression dataframe
    expr_df = pd.DataFrame(
        adata[:, available_genes].X.toarray(),
        columns=available_genes
    )
    expr_df['region'] = adata.obs[region_col].values
    
    # Remove na regions
    expr_df = expr_df[expr_df['region'] != 'na']
    
    # Calculate mean expression per region
    mean_expr = expr_df.groupby('region')[available_genes].mean()
    
    # Calculate fold changes relative to overall mean
    overall_mean = expr_df[available_genes].mean()
    fold_changes = mean_expr.div(overall_mean, axis=1)
    
    print("\n" + "="*80)
    print("GENE EXPRESSION ANALYSIS RESULTS")
    print("="*80)
    
    print(f"\nMean Expression per Region:")
    print(mean_expr.round(3))
    
    print(f"\nFold Change relative to overall mean:")
    print(fold_changes.round(3))
    
    # Identify strongly upregulated genes (>1.5x) per region
    print(f"\nStrongly UPREGULATED genes per region (>1.5x fold change):")
    for region in mean_expr.index:
        upregulated = fold_changes.loc[region][fold_changes.loc[region] > 1.5]
        if len(upregulated) > 0:
            print(f"  {region}: {', '.join([f'{gene} ({fc:.2f}x)' for gene, fc in upregulated.items()])}")
        else:
            print(f"  {region}: None")
    
    # Identify strongly downregulated genes (<0.67x) per region
    print(f"\nStrongly DOWNREGULATED genes per region (<0.67x fold change):")
    for region in mean_expr.index:
        downregulated = fold_changes.loc[region][fold_changes.loc[region] < 0.67]
        if len(downregulated) > 0:
            print(f"  {region}: {', '.join([f'{gene} ({fc:.2f}x)' for gene, fc in downregulated.items()])}")
        else:
            print(f"  {region}: None")
    
    return mean_expr, fold_changes

def plot_gene_expression_boxplots(adata, genes, region_col='annotations'):
    """
    Create boxplots showing gene expression across different regions
    """
    available_genes = [g for g in genes if g in adata.var_names]
    
    # Create expression dataframe
    expr_df = pd.DataFrame(
        adata[:, available_genes].X.toarray(),
        columns=available_genes
    )
    expr_df['region'] = adata.obs[region_col].values
    
    # Remove na regions
    expr_df = expr_df[expr_df['region'] != 'na']
    
    # Melt for plotting
    expr_melted = expr_df.melt(id_vars=['region'], var_name='gene', value_name='expression')
    
    # Create subplot for each gene
    n_genes = len(available_genes)
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()
    
    for i, gene in enumerate(available_genes):
        if i < len(axes):
            gene_data = expr_melted[expr_melted['gene'] == gene]
            sns.boxplot(data=gene_data, x='region', y='expression', ax=axes[i])
            axes[i].set_title(f'{gene} Expression Across Regions', fontweight='bold')
            axes[i].set_xlabel('Region')
            axes[i].set_ylabel('Log Expression')
            axes[i].tick_params(axis='x', rotation=45)
    
    # Hide unused subplots
    for i in range(len(available_genes), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def plot_fancy_gene_expression_boxplots(adata, genes, region_col='annotations'):
    """
    Create fancy boxplots showing gene expression across different regions with tissue-specific colors
    """
    # Filter to only the specified genes
    selected_genes = ["EGR1"]
    available_genes = [g for g in selected_genes if g in adata.var_names]
    
    if not available_genes:
        print("None of the selected genes found in dataset!")
        return None
    
    print(f"Creating fancy boxplots for: {available_genes}")
    
    # Define tissue-specific colors based on your legend
    tissue_colors = {
        'Invasive Carcinoma': '#2E8B57',  # Dark green
        'Necrosis': '#FF6347',            # Orange-red  
        'Fat': '#87CEEB',                 # Light blue
        'Fibrous Tissue': '#9370DB',      # Purple
        'Immune Cell': '#FFD700',         # Gold/Yellow
        'na': '#D3D3D3'                   # Light gray for missing
    }
    
    # Create expression dataframe
    expr_df = pd.DataFrame(
        adata[:, available_genes].X.toarray(),
        columns=available_genes
    )
    expr_df['region'] = adata.obs[region_col].values
    
    # Remove na regions for cleaner visualization
    expr_df = expr_df[expr_df['region'] != 'na']
    
    # Get unique regions (excluding 'na')
    unique_regions = [r for r in expr_df['region'].unique() if r != 'na']
    
    # Set up the plot style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 7))
    fig.suptitle('Gene Expression Across Breast Cancer Tissue Types', 
                 fontsize=20, fontweight='bold', y=0.98)
    
    # Add subtle background
    fig.patch.set_facecolor('#FAFAFA')
    
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        
        # Prepare data for this gene
        gene_data = []
        region_labels = []
        colors_for_plot = []
        
        for region in unique_regions:
            region_expr = expr_df[expr_df['region'] == region][gene].values
            gene_data.append(region_expr)
            region_labels.append(region)
            colors_for_plot.append(tissue_colors.get(region, '#CCCCCC'))
        
        # Create the boxplot
        bp = ax.boxplot(gene_data, 
                       labels=region_labels,
                       patch_artist=True,
                       notch=True,  # Add notches for median confidence intervals
                       showmeans=True,  # Show mean markers
                       meanprops=dict(marker='D', markerfacecolor='white', 
                                    markeredgecolor='black', markersize=6),
                       whiskerprops=dict(linewidth=2),
                       capprops=dict(linewidth=2),
                       medianprops=dict(linewidth=2.5, color='white'),
                       flierprops=dict(marker='o', markerfacecolor='red', 
                                     markersize=4, alpha=0.6))
        
        # Color the boxes according to tissue types
        for patch, color in zip(bp['boxes'], colors_for_plot):
            patch.set_facecolor(color)
            patch.set_alpha(0.8)
            patch.set_linewidth(2)
            patch.set_edgecolor('black')
        
        # Customize the axes
        ax.set_title(f'{gene} Expression', fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel('Tissue Type', fontsize=12, fontweight='bold')
        ax.set_ylabel('Log₁₀(Expression + 1)', fontsize=12, fontweight='bold')
        
        # Rotate x-axis labels for better readability
        ax.tick_params(axis='x', rotation=45, labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        # Set background color
        ax.set_facecolor('#FFFFFF')
        
        # Add statistical annotations
        # Calculate and display mean expression for each tissue
        y_max = ax.get_ylim()[1]
        for j, region in enumerate(unique_regions):
            region_expr = expr_df[expr_df['region'] == region][gene].values
            mean_expr = np.mean(region_expr)
            n_spots = len(region_expr)
            
            # Add text annotation above each box
            ax.text(j + 1, y_max * 0.95, f'n={n_spots}', 
                   ha='center', va='top', fontsize=9, 
                   fontweight='bold', alpha=0.7)
        
        # Enhance the appearance
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('#333333')
    
    # Create a custom legend for tissue types
    legend_elements = []
    for tissue, color in tissue_colors.items():
        if tissue != 'na' and tissue in unique_regions:
            legend_elements.append(mpatches.Patch(color=color, label=tissue, alpha=0.8))
    
    # Add legend to the figure
    fig.legend(handles=legend_elements, 
              loc='center right', 
              bbox_to_anchor=(0.98, 0.5),
              fontsize=12,
              title='Tissue Types',
              title_fontsize=14,
              frameon=True,
              fancybox=True,
              shadow=True)
    
    # Adjust layout to accommodate legend
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    
    # Add some statistics text
    fig.text(0.02, 0.02, 
             'Boxes show quartiles, whiskers show 1.5×IQR, diamonds show means, notches show median CI',
             fontsize=10, style='italic', alpha=0.7)
    
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*80)
    print("EXPRESSION SUMMARY STATISTICS")
    print("="*80)
    
    for gene in available_genes:
        print(f"\n{gene} Expression Summary:")
        print("-" * 40)
        
        gene_summary = expr_df.groupby('region')[gene].agg([
            'count', 'mean', 'std', 'median', 'min', 'max'
        ]).round(3)
        
        print(gene_summary)
        
        # Calculate statistical significance between Invasive Carcinoma and other tissues
        if 'Invasive Carcinoma' in unique_regions:
            cancer_expr = expr_df[expr_df['region'] == 'Invasive Carcinoma'][gene].values
            
            print(f"\nStatistical tests comparing Invasive Carcinoma vs other tissues for {gene}:")
            print("-" * 60)
            
            from scipy.stats import ttest_ind, mannwhitneyu
            
            for region in unique_regions:
                if region != 'Invasive Carcinoma':
                    other_expr = expr_df[expr_df['region'] == region][gene].values
                    
                    # Perform t-test
                    t_stat, t_pval = ttest_ind(cancer_expr, other_expr)
                    
                    # Perform Mann-Whitney U test (non-parametric)
                    u_stat, u_pval = mannwhitneyu(cancer_expr, other_expr, alternative='two-sided')
                    
                    # Calculate fold change
                    cancer_mean = np.mean(cancer_expr)
                    other_mean = np.mean(other_expr)
                    fold_change = cancer_mean / other_mean if other_mean > 0 else np.inf
                    
                    print(f"  vs {region}:")
                    print(f"    Fold Change: {fold_change:.2f}")
                    print(f"    T-test p-value: {t_pval:.2e}")
                    print(f"    Mann-Whitney p-value: {u_pval:.2e}")
                    
                    # Interpretation
                    if t_pval < 0.001:
                        significance = "highly significant (***)"
                    elif t_pval < 0.01:
                        significance = "very significant (**)"
                    elif t_pval < 0.05:
                        significance = "significant (*)"
                    else:
                        significance = "not significant"
                    
                    direction = "upregulated" if fold_change > 1 else "downregulated"
                    print(f"    Result: {gene} is {direction} in cancer ({significance})")
                    print()
    
    return fig


# ===============================================================================
# RUN THE ANALYSIS
# ===============================================================================

print("\n" + "="*80)
print("STARTING GENE-REGION CORRELATION ANALYSIS")
print("="*80)

# Calculate correlations
corr_df, pval_df = calculate_gene_region_correlations(adata, show_genes)

if corr_df is not None:
    # Plot correlation heatmap
    plot_correlation_heatmap(corr_df, pval_df, "Gene-Cancer Region Correlations")
    
    # Analyze expression patterns
    mean_expr, fold_changes = analyze_expression_patterns(adata, show_genes)
    
    # Create boxplots
    #plot_gene_expression_boxplots(adata, show_genes)
    # Usage - replace the original boxplot function call with this:
    plot_fancy_gene_expression_boxplots(adata, show_genes)    
    # Print detailed correlation results
    print(f"\nDETAILED CORRELATION RESULTS:")
    print("="*50)
    for gene in corr_df.index:
        print(f"\n{gene}:")
        for region in corr_df.columns:
            corr = corr_df.loc[gene, region]
            pval = pval_df.loc[gene, region]
            significance = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
            direction = "↑" if corr > 0 else "↓"
            print(f"  {region}: r={corr:.3f} (p={pval:.3e}) {significance} {direction}")
    
    # Save results to CSV
    corr_df.to_csv('gene_region_correlations.csv')
    pval_df.to_csv('gene_region_pvalues.csv')
    fold_changes.to_csv('gene_region_fold_changes.csv')
    
    print(f"\nResults saved to:")
    print(f"  - gene_region_correlations.csv")
    print(f"  - gene_region_pvalues.csv") 
    print(f"  - gene_region_fold_changes.csv")

# Plot spatial expression for each gene
print(f"\nPlotting spatial expression for marker genes...")
sc.pl.spatial(adata, color=show_genes, spot_size=275, frameon=False, ncols=2, alpha_img=0.7)

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)