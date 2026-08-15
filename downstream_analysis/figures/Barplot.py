import matplotlib.pyplot as plt
import seaborn as sns

# Define the data
cluster_labels = ['Cancer region', 'Stroma','Duct epithelium', 'Pancreatic tissue' ]
num_cells = [102, 54, 5, 34]
percentages = [52.3, 27.7, 2.6, 17.4]

# Create the figure and axes
fig, ax = plt.subplots(figsize=(3, 3.5))

# Set the color palette
palette = sns.color_palette("muted", len(cluster_labels))
palette = [palette[0], palette[3], palette[1], palette[2]]  # Custom order

# Plot the bar chart with reduced spacing between bars
ax.bar(cluster_labels, num_cells, edgecolor='black', color=palette, width=0.8)  # Thinner bars, less space

# Add data labels to the bars with both numbers and percentages above the bars
for i, cell_count in enumerate(num_cells):
    ax.text(i, cell_count + 6, f"{cell_count}", ha='center', va='bottom', fontsize=7)  # Numbers above
    ax.text(i, cell_count + 1, f"({percentages[i]}%)", ha='center', va='bottom', fontsize=7)  # Percentages further above

# Remove the top, right, and bottom spines (but keep the left/y-axis spine)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Make the y-axis thicker
ax.spines['left'].set_linewidth(2.0)  # Thicker y-axis

# Rotate the x-axis labels by 45 degrees
plt.xticks(rotation=45, ha='right', fontsize = 8)

# Set the x-axis label
#ax.set_xlabel('Cluster')

# Set the y-axis label
ax.set_ylabel('Number of spots')

# Set the title
#ax.set_title('Cell Count per Cluster')

# Adjust the spacing
plt.tight_layout()

# Show the plot
plt.show()
