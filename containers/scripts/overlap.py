#!/usr/bin/env python3

import sys

if len(sys.argv) != 4:
    sys.exit("\nUsage: python overlaphash.py <infile> <infile1> <outfile>\n")

infile = sys.argv[1]
infile1 = sys.argv[2]
outfile = sys.argv[3]

# Step 1: Read infile1 into a dictionary
clusteringLabel = {}
with open(infile1, 'r') as f1:
    for line in f1:
        fields = line.strip().split('\t')
        if len(fields) >= 2:
            key = fields[0]
            value = fields[1]
            clusteringLabel[key] = value

# Step 2: Read infile and match keys
with open(infile, 'r') as f, open(outfile, 'w') as out:
    # Write header
    out.write("sample\tchr\tstart\tend\tstate\tlabel\n")
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 1:
            key = fields[0]
            if key in clusteringLabel:
                out.write(f"{line.strip()}\t{clusteringLabel[key]}\n")
