import sys

if len(sys.argv) != 4:
    sys.exit("\nUsage: <infile> <infile1> <outfile>\n\n")

infile = sys.argv[1]
infile1 = sys.argv[2]
outfile = sys.argv[3]

# Read the second input file into a dictionary
gene = {}
with open(infile1, 'r') as f:
    for line in f:
        line = line.strip()
        gene[line] = line

# Process the first input file
with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    n = 0
    for current_line in fin:
        current_line = current_line.strip()
        cnvgene = []
        brokengene = []

        if n == 0:
            fout.write(f"SampleName\tIntervalSequentialID\tchr\tstart\tend\tstate\tmedian\tGenes\tBroken_Genes\n")
        else:
            s = current_line.split("\t")

            if s[5] == '"3"':
                fout.write(f"{current_line}\teuploid\teuploid\n")
            else:
                for key in gene:
                    s1 = gene[key].split("\t")
                    s[2] = s[2].replace('"', '')
                    if s[2] == s1[0] and float(s[3]) < float(s1[1]) and float(s[4]) > float(s1[2]):
                        cnvgene.append(f"{s1[3]}:{s1[4]}")
                    elif s[2] == s1[0] and float(s[3]) > float(s1[1]) and float(s[3]) < float(s1[2]):
                        cnvgene.append(f"{s1[3]}:{s1[4]}")
                        brokengene.append(f"{s1[3]}:{s1[4]}")
                    elif s[2] == s1[0] and float(s[4]) > float(s1[1]) and float(s[4]) < float(s1[2]):
                        cnvgene.append(f"{s1[3]}:{s1[4]}")
                        brokengene.append(f"{s1[3]}:{s1[4]}")

                fout.write(f"{current_line}\t")
                fout.write(";".join(cnvgene[:-1]) + (";" if len(cnvgene) > 1 else ""))
                if cnvgene:
                    fout.write(cnvgene[-1])
                fout.write("\t")
                fout.write(";".join(brokengene[:-1]) + (";" if len(brokengene) > 1 else ""))
                if brokengene:
                    fout.write(brokengene[-1])
                else:
                    print("NA")
                    fout.write("NA")
                fout.write("\n")

        n += 1

