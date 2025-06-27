import sys

if len(sys.argv) != 4:
    print("\nUsage: <infile> <infile1> <outfile>\n")
    sys.exit(1)

infile = sys.argv[1]
infile1 = sys.argv[2]
outfile = sys.argv[3]

try:
    with open(infile, 'r') as instream, open(infile1, 'r') as instream1, open(outfile, 'w') as out:
        out.write("chr\tstart\tend\tgene\n")

        s = {}
        for line1 in instream1:
            line1 = line1.strip().rstrip('\r')
            parts1 = line1.split('\t')
            pos = parts1[4].upper()
            s[pos] = f"{parts1[0]}\t{parts1[1]}\t{parts1[2]}"

        for line in instream:
            line = line.strip().rstrip('\r')
            parts = line.split('\t')
            index = parts[0].upper()
            if index in s:
                out.write(f"{s[index]}\t{parts[0]}\n")

except FileNotFoundError as e:
    print(f"Couldn't open file: {e.filename}")
    sys.exit(1)
