import argparse
import csv

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Aneuploidy portion calculation.")
parser.add_argument("infile", help="merge.segments.tsv)")
parser.add_argument("infile1", help="chr length file")
parser.add_argument("outfile", help="Output TSV file path")
args = parser.parse_args()

# Step 1: Read segment length file into a dictionary
cnv = {}
with open(args.infile1, 'r') as f2:
    for line in f2:
        line = line.strip().rstrip("\r")
        partscnv = line.split("\t")
        if len(partscnv) >= 3:
            pos = partscnv[1]
            value = partscnv[2]
            cnv[pos] = value

# Step 2: Read chr length file and annotate with segment length data, and calculate cnv segment portion
with open(args.infile, 'r') as f1, open(args.outfile, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    n = 0
    for line in f1:
        line = line.strip().rstrip("\r")
        parts = line.split("\t")
        index = parts[0]

        if index in cnv:
            if n == 0:
                writer.writerow(["sampleID"] + parts + [cnv[index], "CNVscore_aneuploidy_portion"])
            else:
                value = float(cnv[index])
                denom = float(parts[1]) if len(parts) > 1 and parts[1].replace('.', '', 1).isdigit() else 1
                portion = value / denom if denom != 0 else 0
                writer.writerow([partscnv[0]] + parts + [cnv[index], portion])
        else:
            writer.writerow([partscnv[0]] + parts + ["0", "0"])
        n += 1

