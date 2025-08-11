#!/bin/bash
# 
# Change these options depending on
# your HPC environment's policies
# 
#SBATCH -N 1
#SBATCH -n 4

DATA_DIR="$PWD/data/"
OUTDIR="results"
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

nextflow run StreamlineCNV.nf \
            -profile slurm,singularity \
            --data_dir $DATA_DIR \
            --fastq "$DATA_DIR/*fq" \
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
