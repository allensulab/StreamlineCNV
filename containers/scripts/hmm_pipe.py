import os
import sys
import argparse
import subprocess
import shutil
import statistics
from pathlib import Path

def check_programs():
    programs = ['R', 'samtools', 'dirname', 'readCounter']
    for program in programs:
        if shutil.which(program) is None:
            raise RuntimeError(f"CHECK: {program} not found")

def average(values):
    return sum(values) / len(values)

def get_segstr(inputfile):
    with open(inputfile) as fi:
        next(fi)  # Skip header
        segments = ""
        for line in fi:
            segments += f"{inputfile}\t{line}"
        return segments

def get_vs(inputfile):
    winsize = 30
    values = {}
    with open(inputfile) as fi:
        for line in fi:
            line = line.replace('"', '').strip()
            chr_, v = line.split("\t")
            if v != "NA":
                values.setdefault(chr_, []).append(float(v))

    averages = []
    for chr_, vals in values.items():
        if chr_ in ["chrX", "chrY"]:
            continue

        sds = []
        for i in range(len(vals) - winsize):
            window = vals[i:i + winsize]
            mean = average(window)
            variance = sum((x - mean) ** 2 for x in window) / winsize
            sds.append(variance ** 0.5)
        if sds:
            averages.append(average(sds))

    return average(sorted(averages, reverse=True)[:3]) if averages else 0

def main():
    parser = argparse.ArgumentParser(description="Python script to run HMMpipe")
    parser.add_argument("input_bam_folder", help="Folder containing BAM files")
    parser.add_argument("ref_genome", help="Reference genome (e.g., GRCh38)")
    parser.add_argument("--downsample", type=int, default=0, help="Downsample to this number of reads (default: 0)")
    parser.add_argument("--dosomeclean", action="store_true", help="Perform additional cleaning (remove chrM, sorting)")
    parser.add_argument("--evalue", type=float, default=0.9999999, help="E-value for HMMcopy (default: 0.9999999)")
    parser.add_argument("--keeptmp", action="store_true", help="Keep temporary files for debugging")
    args = parser.parse_args()

    check_programs()

    script_dir = Path(__file__).parent
    input_folder = Path(args.input_bam_folder)
    ref_genome = args.ref_genome

    if not input_folder.is_dir():
        raise FileNotFoundError(f"Input folder {input_folder} not found")

    input_files = sorted([file for file in input_folder.iterdir() if file.suffix == ".bam"], key=lambda x: x.name)
    total_reads = {}
    mapped_reads = {}
    vss = {}
    all_segments = ""

    for bam_file in input_files:
        print(f"Processing {bam_file}", file=sys.stderr)

        sampnum = int(subprocess.getoutput(f"samtools view -F {0x100 |0x800} -c {bam_file}"))
        mappednum = int(subprocess.getoutput(f"samtools view -F {0x0004 | 0x100 | 0x800} -c {bam_file}"))

        total_reads[bam_file.name] = sampnum
        mapped_reads[bam_file.name] = mappednum

        tmpfolder = bam_file.stem + ".tmp"
        if Path(tmpfolder).exists():
            shutil.rmtree(tmpfolder)
        Path(tmpfolder).mkdir()

        if args.downsample > 0:
            frac = args.downsample / sampnum if sampnum > 0 else 0
            if sampnum == 0 or frac > 1:
                print(f"WARNING: {frac} greater than 1: no downsampling performed", file=sys.stderr)
            else:
                downsampled_bam = Path(tmpfolder) / f"{bam_file.stem}.down.bam"
                subprocess.run(f"samtools view -F {0x0004 | 0x100 | 0x800} -b -s {frac} {bam_file} > {downsampled_bam}", shell=True, check=True)
                bam_file = downsampled_bam

        if args.dosomeclean:
            clean_sam = Path(tmpfolder) / f"{bam_file.stem}.tmp.sam"
            clean_bam = Path(tmpfolder) / f"{bam_file.stem}.tmp.bam"
            subprocess.run(f"samtools view -h {bam_file} | grep -v -w MT | grep -v KI|grep -v GL|grep -v JH|grep -v MU > {clean_sam}", shell=True, check=True)
            subprocess.run(f"samtools view -Shu {clean_sam} | samtools sort -T {tmpfolder}/{bam_file.stem} -o {clean_bam}", shell=True, check=True)
            bam_file = clean_bam

        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
        wig_file = Path(tmpfolder) / "input.wig"
        subprocess.run(f"readCounter -w 500000 {bam_file} > {wig_file}", shell=True, check=True)
        r_script = script_dir / f"run_hmmcopy.{ref_genome}.r"
        subprocess.run(f"Rscript {r_script} {args.evalue}", cwd=tmpfolder, shell=True, check=True)

        # Move the single.fixed.y.pdf to the parent directory
        pdf_file = Path(tmpfolder) / "single.fixed.y.pdf"
        if pdf_file.exists():
            dest_pdf = Path.cwd() / f"{bam_file.stem}.pdf"
            shutil.move(str(pdf_file), str(dest_pdf))
            print(f"PDF file moved to: {dest_pdf}")
        else:
            print(f"Error: PDF file not found in {tmpfolder}", file=sys.stderr)

        vss[bam_file.name] = get_vs(Path(tmpfolder) / "data.txt")
        all_segments += get_segstr(Path(tmpfolder) / "segments.txt")

        if not args.keeptmp:
            shutil.rmtree(tmpfolder)

    print("\t".join(["sample", "totalreads", "mappedreads", "mapped%", "VS"]))
    for sample, total in total_reads.items():
        mapped = mapped_reads[sample]
        print("\t".join(map(str, [sample, total, mapped, mapped / total if total > 0 else 0, vss[sample]])))

    with open("merge.segments.txt", "w") as seg_file:
        seg_file.write(all_segments)

if __name__ == "__main__":
    main()
