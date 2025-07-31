import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import sys
from sklearn.preprocessing import LabelEncoder

# Read command-line argument
if len(sys.argv) < 2:
    print("Usage: python script.py output_filename.pdf")
    sys.exit(1)

output_file = sys.argv[1]

# Load segment annotation
df = pd.read_csv("seg_anno", sep="\t")

# Expand segments into bins
bin_size = 500000
expanded_rows = []
for _, row in df.iterrows():
    start, end = int(row['start']), int(row['end'])
    for b_start in range(start, end, bin_size):
        b_end = min(b_start + bin_size, end)
        if b_end > b_start:
            new_row = row.copy()
            new_row['start'] = b_start
            new_row['end'] = b_end
            new_row['segment'] = f"{row['chr']}:{b_start}-{b_end}"
            expanded_rows.append(new_row)
df_exp = pd.DataFrame(expanded_rows)

# Pivot to matrix
matrix = df_exp.pivot_table(index='sample', columns='segment', values='state', aggfunc='first').fillna(3)
sample_labels = df[['sample', 'label']].drop_duplicates().set_index('sample').loc[matrix.index]

# Sort segments
segment_info = matrix.columns.str.extract(r'(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+)')
segment_info['start'] = segment_info['start'].astype(int)
segment_info['chr'] = pd.Categorical(segment_info['chr'], categories=[str(i) for i in range(1, 23)] + ['X', 'Y'], ordered=True)
segment_info['original'] = matrix.columns
segment_info = segment_info.sort_values(['chr', 'start'])
matrix_sorted = matrix[segment_info['original'].tolist()]

# Set sample labels and colors
label_encoder = LabelEncoder()
sample_labels['label_code'] = label_encoder.fit_transform(sample_labels['label'])
unique_labels = label_encoder.classes_
palette = sns.color_palette("Set2", len(unique_labels))
label_to_color = dict(zip(unique_labels, palette))
row_colors = sample_labels['label'].map(label_to_color)
row_colors.index.name = None     # Remove index name
row_colors.name = None           # Remove Series name
matrix_sorted.index.name = None

# Generate heatmap
sns.set(style="white")
plt.rcParams.update({'font.family': 'DejaVu Sans'})

g = sns.clustermap(
    matrix_sorted,
    cmap=["white"],
    row_cluster=True,
    col_cluster=False,
    row_colors=row_colors,
    xticklabels=False,
    yticklabels=False,
    figsize=(18, 10),
    dendrogram_ratio=(0.1, 0.01),
    cbar_pos=None
)

# Apply clustering order
clustered_order = g.dendrogram_row.reordered_ind
clustered_samples = matrix_sorted.index[clustered_order]
df_exp = df_exp[df_exp['sample'].isin(clustered_samples)].copy()

# Center sample labels on y-axis
g.ax_heatmap.set_yticks(np.arange(len(clustered_samples)) + 0.5)
g.ax_heatmap.set_yticklabels(clustered_samples, fontsize=8)

# Colorbar setup
state_labels = ['homozygous deletion', 'heterozygous deletion', 'neutral', 'gain', 'amplification', 'high-level amplification']
state_colors = ['#77dd77', '#c2f0c2', 'white', '#f7c2c2', '#f77e7e', '#ff0000']
cmap_custom = mcolors.ListedColormap(state_colors)
norm_custom = mcolors.BoundaryNorm(boundaries=[1, 2, 3, 4, 5, 6, 7], ncolors=6)


cb_ax = g.fig.add_axes([1.07, 0.60, 0.2, 0.02])  # move color palette on top of sample legend
cb = plt.colorbar(
    plt.cm.ScalarMappable(cmap=cmap_custom, norm=norm_custom),
    cax=cb_ax,
    orientation='horizontal',
)
cb.ax.tick_params(axis='x', length=0)

for i, label in enumerate(state_labels):
    xpos = (i + 0.5) / 6
    if (label=="homozygous deletion"):
        xpos=xpos+0.11
    elif label == "heterozygous deletion":
        xpos=xpos+0.11
    elif label == "neutral":
        xpos=xpos+0.02
    elif label == "amplification":
        xpos=xpos+0.08
    elif label == "high-level amplification":
        xpos=xpos+0.12
    cb_ax.text(xpos, 1.2, label, ha='center', va='bottom', fontsize=8, rotation=45, transform=cb_ax.transAxes)
handles = [mpatches.Patch(color=palette[i], label=label) for i, label in enumerate(unique_labels)]
g.ax_heatmap.legend(
    handles=handles,
    loc='upper left',
    bbox_to_anchor=(1.10, 0.5),  # move sample legend slightly left
    fontsize=12,
    frameon=False
)

# Patch segments to heatmap
sample_order = {sample: i for i, sample in enumerate(clustered_samples)}
df_exp['sample_index'] = df_exp['sample'].map(sample_order)
segment_locs = dict(zip(segment_info['original'], range(len(segment_info))))
df_exp['segment_id'] = df_exp['chr'].astype(str) + ":" + df_exp['start'].astype(str) + "-" + df_exp['end'].astype(str)
df_exp.drop_duplicates(subset=['sample', 'segment_id'], inplace=True)

for _, row in df_exp.iterrows():
    if row['chr'] == 'Y':
        continue
    segment = f"{row['segment']}" 
    seg_idx = segment_locs.get(segment)


    row_idx = row['sample_index']
    if seg_idx is not None and row_idx is not None:
        try:
            color = cmap_custom.colors[int(row['state']) - 1]
            g.ax_heatmap.add_patch(Rectangle((seg_idx, row_idx), 1, 1, color=color, linewidth=0))
        except:
            pass

# Draw chromosome boundaries and labels
chr_boundaries = segment_info.groupby('chr').size().cumsum()
chr_starts = chr_boundaries.shift(fill_value=0)
chr_midpoints = (chr_starts + chr_boundaries) / 2
for chr_name, pos in chr_boundaries.items():
    if chr_name != 'Y':
        g.ax_heatmap.axvline(x=pos, color='black', linewidth=0.3, zorder=100)

# Label chromosome axis
g.ax_heatmap.set_xlabel('chromosome', labelpad=25)
g.ax_heatmap.set_ylabel('')
for pos, label in zip(chr_midpoints, chr_midpoints.index):
    g.ax_heatmap.text(pos, matrix_sorted.shape[0]+0.3, label, ha='center', va='top', fontsize=8, rotation=0)

# Add legend
handles = [mpatches.Patch(color=label_to_color[label], label=label) for label in unique_labels]
g.ax_heatmap.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.10, 0.5), fontsize=12, frameon=False)

# Save output
plt.savefig(output_file, bbox_inches='tight')
