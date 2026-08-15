import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Data for the stacked bars
data = np.array([
    [2.8, 72.3, 24.8],   # A
    [24.0, 32.3, 43.7],   # B
    [54.2, 8.5, 37.3],  # C
    [11.5, 55.7, 32.8]  # D
])

# Labels for x-axis
labels = ['Cancer region', 'Stroma','Duct epithelium', 'Pancreatic tissue']

# Create the stacked bar plot
fig, ax = plt.subplots(figsize=(3, 3.5))

# Define specific colors
colors = ['blue', 'red', 'gray']  # Light green, light orange, light blue

bottom = np.zeros(4)

for i, color in enumerate(colors):
    values = data[:, i]
    ax.bar(labels, values, bottom=bottom, color=color)
    
    # Add percentage labels inside the bars
    for j, value in enumerate(values):
        if value > 0:  # Only add label if the value is greater than 0
            ax.text(j, bottom[j] + value/2, f'{value:.1f}%', 
                    ha='center', va='center', fontsize=6, color='white',
                    fontweight='bold')
    
    bottom += values


# Customize the plot
ax.set_ylim(0, 100)
ax.set_ylabel('Percentage')

# Remove the top, right, and bottom spines (but keep the left/y-axis spine)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Make the y-axis thicker
ax.spines['left'].set_linewidth(2.0)  # Thicker y-axis

# Remove x-axis ticks
ax.tick_params(axis='x', which='both', bottom=False, top=False)

# Adjust x-axis labels
plt.xticks(rotation=45, ha='right', fontsize=8)

# Remove y-axis ticks but keep labels
ax.tick_params(axis='y', which='both', left=False, right=False)
ax.yaxis.set_ticks_position('none')

# Display the plot
plt.tight_layout()
plt.show()