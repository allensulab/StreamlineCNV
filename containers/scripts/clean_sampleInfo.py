import sys

if len(sys.argv) != 5:
    sys.exit("\nUsage: <infile> <infile1> <outfile> <outfile1>\n")

infile = sys.argv[1]
infile1 = sys.argv[2]
outfile = sys.argv[3]
outfile1 = sys.argv[4]

wigfiles = {}
with open(infile1, 'r') as f:
    for line in f:
        line = line.strip().rstrip('\r')
        parts = line.split('/')
        s2 = parts[-1].split('.')
        pos = s2[0]
        wigfiles[pos] = line

with open(infile, 'r') as fin, open(outfile, 'w') as fout, open(outfile1, 'w') as fout1:
    fout1.write("Warning: The following samples listed in the sampleInfo files do not have wig files\n")
    for line in fin:
        line = line.strip().rstrip('\r')
        s1 = line.split('\t')
        index = s1[1]
        if index in wigfiles:
            fout.write(line + '\n')
        else:
            fout1.write(s1[1] + '\n')

