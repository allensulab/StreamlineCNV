#!/bin/bash
# 
# Change these options depending on
# your HPC environment's policies
# 
#SBATCH -N 1
#SBATCH -n 4

REF_FILES_DIR="$PWD/data/ref_files"
FASTQ_DATA="$PWD/data/fastq_files"
OUTDIR="$PWD/results"
GENE_LIST="$PWD/data/GeneListExample"
SAMPLE_INFO="$PWD/data/SampleInfo"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
QUEUE="normal"
EXTRA_MOUNT="/net:/net"

# Make sure to include any steps necessary
# to load nextflow and singularity into
# your HPC environment. For example:
# 
# module load singularity
# module load nextflow

nextflow clean -f
nextflow run StreamlineCNV.nf \
            -profile slurm,singularity \
            --ref_files_dir $REF_FILES_DIR \
            --fastq "$FASTQ_DATA/*fq" \
            --species 'Homo_sapiens' \
            --assembly 'GRCh38' \
            --outdir $OUTDIR \
            --plotGeneDensity true \
            --geneList $GENE_LIST \
            --recolor true \
            --sampleInfo $SAMPLE_INFO \
            --clustering true \
            --clusteringLabel $CLUSTERING_LABEL \
            --queue $QUEUE \
            --extra_mount $EXTRA_MOUNT \
            -resume
