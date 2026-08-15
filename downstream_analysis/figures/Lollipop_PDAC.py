import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("🚀 Creating Median-Only Lollipop Plot...")

# Load data
print("📂 Loading data...")
pathway_df = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-PDAC/Graphencode/latent/pathway_encoded_df_2024.csv', index_col=0)
annotation_df = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Dataset/Raw-PDAC/PDAC/PDAC-A ST1/STAGE-PDAC-A-ST1-meta.txt', sep='\t', index_col=0)

print(f"✅ Data loaded: {pathway_df.shape[0]} spots, {pathway_df.shape[1]} pathways")

# Choose pathway to plot (you can change this)
#pathway_name = pathway_df.columns[0]  # First pathway
pathway_name = 'REACTOME_IMMUNE_SYSTEM'  # Or specify a particular one

print(f"📊 Plotting pathway: {pathway_name}")

# Get pathway scores and annotations
pathway_scores = pathway_df[pathway_name]
annotations = annotation_df['human_anno_region']

# Create aggregated data by annotation group
aggregated_data = []
for spot_id, score in pathway_scores.items():
    annotation = annotations[spot_id]
    aggregated_data.append({
        'score': score,
        'annotation': annotation
    })

data_df = pd.DataFrame(aggregated_data)

# Calculate median and mean for each annotation group
summary_stats = data_df.groupby('annotation')['score'].agg(['median', 'mean', 'std', 'count']).reset_index()
summary_stats = summary_stats.sort_values('annotation')  # Sort alphabetically

print(f"✅ Created aggregated data for {len(summary_stats)} annotation groups")

# Define colors for each tissue type
color_map = {
    'Cancer region': '#1f77b4',      # Blue
    'Stroma': '#d62728',             # Red  
    'Duct epithelium': '#ff7f0e',    # Orange
    'Pancreatic tissue': '#2ca02c'   # Green
}

print("🎨 Creating median-only lollipop plot...")

# Create the plot with median only
fig, ax = plt.subplots(figsize=(12, 8))

# Set up positions for bars
n_groups = len(summary_stats)
x_positions = np.arange(n_groups)

# Plot median bars
for i, row in summary_stats.iterrows():
    color = color_map[row['annotation']]
    
    # Plot median stem and head
    ax.plot([x_positions[i], x_positions[i]], 
           [0, row['median']], 
           color=color, 
           linewidth=2.5, 
           alpha=0.8,
           label='_nolegend_')
    
    ax.scatter(x_positions[i], row['median'], 
              color=color, 
              s=100,
              alpha=0.9,
              edgecolors='white',
              linewidths=1.5,
              marker='o',
              label=row['annotation'] if i == 0 else '_nolegend_')

# Customize the plot
ax.axhline(y=0, color='black', linewidth=0.8, alpha=0.6)
ax.set_xlabel('Tissue Regions', fontsize=14, fontweight='bold')
ax.set_ylabel('Pathway Activity Score', fontsize=14, fontweight='bold')
ax.set_title(f'Median Pathway Activity by Tissue Region\n{pathway_name}', 
             fontsize=16, fontweight='bold', pad=20)

# Set x-axis
ax.set_xticks(x_positions)
ax.set_xticklabels(summary_stats['annotation'], rotation=45, ha='right')

# Create custom legend
legend_elements = []
# Add tissue type legend
for annotation in summary_stats['annotation']:
    legend_elements.append(plt.scatter([], [], color=color_map[annotation], s=60, 
                                     label=f'{annotation}', edgecolors='white', linewidths=1))

ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1), loc='upper left', 
          frameon=True, fancybox=True, shadow=True)

# Customize grid and appearance
ax.grid(True, alpha=0.3, axis='y', linestyle='-', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_linewidth(0.8)

# Set background color
ax.set_facecolor('#fafafa')
fig.patch.set_facecolor('white')

plt.tight_layout()

# Save the plot
fig.savefig('median_lollipop_pathway_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
print("✅ Plot saved as 'median_lollipop_pathway_plot.png'")

# Show the plot
plt.show(block=False)
plt.draw()
print("✅ Plot displayed!")

# Print summary
print(f"\n📋 Plot Summary:")
print(f"   Pathway: {pathway_name}")
print(f"   Total annotation groups: {len(summary_stats)}")
print(f"\n   Summary Statistics:")
for _, row in summary_stats.iterrows():
    print(f"   {row['annotation']}:")
    print(f"     Count: {row['count']} spots")
    print(f"     Median: {row['median']:.3f}")
    print(f"     Mean: {row['mean']:.3f}")
    print(f"     Std: {row['std']:.3f}")
    print()

print("\n" + "="*50)
print("CREATING MULTIPLE PATHWAY SUBPLOT (MEDIAN ONLY)")
print("="*50)

# Create multiple pathway plot with aggregation
#CANCER:
#pathways_to_plot = ['REACTOME_CELL_JUNCTION_ORGANIZATION','REACTOME_IMMUNE_SYSTEM']
#DUCT
#pathways_to_plot = ['REACTOME_ER_PHAGOSOME_PATHWAY','REACTOME_INNATE_IMMUNE_SYSTEM']
#PANC
#pathways_to_plot = ['REACTOME_SYNTHESIS_OF_PA','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION']
#Stroma
pathways_to_plot = ['REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR','REACTOME_CELL_JUNCTION_ORGANIZATION']

n_pathways = len(pathways_to_plot)

fig, axes = plt.subplots(1, 2, figsize=(5, 3))  # Single row, 6 columns

for idx, current_pathway in enumerate(pathways_to_plot):
    ax = axes[idx]
    
    # Get scores for this pathway
    scores = pathway_df[current_pathway]
    
    # Create aggregated data for this pathway
    pathway_data = []
    for spot_id, score in scores.items():
        annotation = annotations[spot_id]
        pathway_data.append({
            'score': score,
            'annotation': annotation
        })
    
    pathway_df_subset = pd.DataFrame(pathway_data)
    pathway_summary = pathway_df_subset.groupby('annotation')['score'].agg(['median']).reset_index()
    pathway_summary = pathway_summary.sort_values('annotation')
    
    # Plot aggregated data
    n_groups_sub = len(pathway_summary)
    x_pos_sub = np.arange(n_groups_sub)
    
    # Plot median only
    for i, row in pathway_summary.iterrows():
        color = color_map[row['annotation']]
        
        # Median
        ax.plot([x_pos_sub[i], x_pos_sub[i]], [0, row['median']], 
               color=color, linewidth=1.5, alpha=0.8)
        ax.scatter(x_pos_sub[i], row['median'], color=color, s=40, 
                  alpha=0.9, edgecolors='white', linewidths=0.5, marker='o')
    
    # Customize subplot
    ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.5)
    ax.set_title(current_pathway.replace('REACTOME_', '').replace('_', ' '), 
                fontsize=7, fontweight='bold')
    
    # Set classic/clean style - white background, no grid
    ax.grid(False)  # Remove grid lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.set_facecolor('white')  # White background
    
    # Set x-ticks with numbers only (like in the reference image)
    ax.set_xticks(x_pos_sub)
    ax.set_xticklabels(range(len(pathway_summary)), fontsize=8)
    
    # Make subplots more compact
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_ylim(bottom=min(0, pathway_summary['median'].min() * 1.1),
               top=pathway_summary['median'].max() * 1.1)

# Add common labels
fig.text(0.5, 0.02, 'Tissue Regions', ha='center', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(bottom=0.15, top=0.85)

# Save multiple pathway plot
fig.savefig('multiple_pathway_median_lollipop.png', dpi=300, bbox_inches='tight', facecolor='white')
print("✅ Multiple pathway plot saved as 'multiple_pathway_median_lollipop.png'")

plt.show(block=False)
plt.draw()

print(f"\n🎉 All median-only plots created successfully!")
print(f"Files saved:")
print(f"  - median_lollipop_pathway_plot.png (single pathway)")
print(f"  - multiple_pathway_median_lollipop.png (6 pathways)")

# Keep script running to see plots
input("\nPress Enter to exit...")