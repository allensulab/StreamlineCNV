import pandas as pd
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(description="aneuploidy chromosome segment lengths")
parser.add_argument("input_file", help="Path to the input .tsv file")
parser.add_argument("output_file", help="Path to the output summary .tsv file")
parser.add_argument("sample_name", help="Sample name to include in the output")
args = parser.parse_args()

# Load the data
df = pd.read_csv(args.input_file, sep="\t")

# Filter rows where state is not 3
df_filtered = df[df['state'] != 3].copy()

# Calculate length of each segment
df_filtered['length'] = df_filtered['end'] - df_filtered['start'] + 1

# Summarize total length per chromosome
summary = df_filtered.groupby('chr')['length'].sum().reset_index()
summary.columns = ['chromosome', 'aneuploidy_length_bp']

# Add total row
total = pd.DataFrame([['genome', summary['aneuploidy_length_bp'].sum()]], columns=summary.columns)
summary_with_total = pd.concat([summary, total], ignore_index=True)

# Add sample name as the first column
summary_with_total.insert(0, 'sample', args.sample_name)

# Save to a file
summary_with_total.to_csv(args.output_file, sep="\t", index=False)
