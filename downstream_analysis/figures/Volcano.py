import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from adjustText import adjust_text

# Set random seed for reproducibility
np.random.seed(42)

# Load the data
#DDX60L_DE = pd.read_csv("/Users/naminiyakan/Documents/BulkToST/Res-PDAC/MYEOV_DE_Res.csv", header=0, index_col=None)
DDX60L_DE = pd.read_csv("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/BRCA_PACSI_DE_Res.csv", header=0, index_col=None)


DDX60L_DE['p_val_adj'] = -np.log10(DDX60L_DE['p_val_adj'])

DDX60L_DE = DDX60L_DE.loc[np.abs(DDX60L_DE['avg_log2FC']) <= 4]

# Define thresholds for significant genes
lfc_threshold =0.5849
#lfc_threshold = 0.3219
fdr_threshold = 1.3

# Categorize genes
DDX60L_DE['category'] = 'Not Significant'
DDX60L_DE.loc[(DDX60L_DE['avg_log2FC'] > lfc_threshold) & (DDX60L_DE['p_val_adj'] > fdr_threshold), 'category'] = 'Up-regulated'
DDX60L_DE.loc[(DDX60L_DE['avg_log2FC'] < -lfc_threshold) & (DDX60L_DE['p_val_adj'] > fdr_threshold), 'category'] = 'Down-regulated'

# Set up the plot
plt.figure(figsize=(6, 4))  # Increased size to accommodate labels
sns.set_style("white")

# Create the scatter plot
scatter = sns.scatterplot(
    data=DDX60L_DE,
    x='avg_log2FC',
    y='p_val_adj',
    hue='category',
    palette={'Up-regulated': 'red', 'Down-regulated': 'blue', 'Not Significant': 'gray'},
    alpha=0.7,
    s=5,
    legend=False
)

# Customize the plot
plt.title('Volcano Plot of Differential Gene Expression', fontsize=16)
plt.xlabel('Log2 Fold Change', fontsize=12)
plt.ylabel('-Log10 FDR', fontsize=12)

# Add threshold lines
plt.axvline(x=lfc_threshold, color='black', linestyle='--', alpha=0.9)
plt.axvline(x=-lfc_threshold, color='black', linestyle='--', alpha=0.9)
plt.axhline(y=fdr_threshold, color='black', linestyle='--', alpha=0.9)

# Create custom legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Up-regulated',
           markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Down-regulated',
           markerfacecolor='blue', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Not Significant',
           markerfacecolor='gray', markersize=8)
]

# Add custom legend
plt.legend(handles=legend_elements, title='Gene Regulation', 
           title_fontsize='12', fontsize='10', 
           bbox_to_anchor=(1.05, 1), loc='upper left')

# Add frame
for spine in plt.gca().spines.values():
    spine.set_visible(True)

# Get top 3 up-regulated and down-regulated genes
top_up = DDX60L_DE[DDX60L_DE['category'] == 'Up-regulated'].nlargest(5, 'p_val_adj')
top_down = DDX60L_DE[DDX60L_DE['category'] == 'Down-regulated'].nlargest(5, 'p_val_adj')



# Prepare texts for adjustment
texts = []

# Add labels for top up-regulated genes
for _, gene in top_up.iterrows():
    texts.append(plt.text(gene['avg_log2FC'], gene['p_val_adj'], gene['gene'],
                 color='red', fontsize=7, fontweight='bold'))

# Add labels for top down-regulated genes
for _, gene in top_down.iterrows():
    texts.append(plt.text(gene['avg_log2FC'], gene['p_val_adj'], gene['gene'],
                 color='blue', fontsize=7, fontweight='bold'))

# Adjust text positions to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            expand_points=(1.5, 1.5), force_points=0.1)

# Adjust layout to prevent legend from being cut off
plt.tight_layout()

# Show the plot
plt.show()