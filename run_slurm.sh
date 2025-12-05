#!/bin/bash
# 
# Change these options depending on
# your HPC environment's policies
# 
#SBATCH -N 1
#SBATCH -n 4

DATA_DIR="$PWD/data/"
FASTQ_DATA="$PWD/data/fastq_files"
OUTDIR="results"
GENE_LIST="$PWD/data/GeneListExample"
SAMPLE_INFO="$PWD/data/SampleInfo"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
QUEUE="normal"
CONTAINER_PROFILE="singularity"
EXTRA_MOUNT="/net:/net"

# Make sure to include any steps necessary
# to load nextflow and your containerization software
# into your HPC environment. For example:
# 
# module load singularity
# module load nextflow

nextflow run StreamlineCNV.nf \
            -profile slurm,$CONTAINER_PROFILE \
            --data_dir $DATA_DIR \
            --fastq "$FASTQ_DATA/*.fastq" \
            --species 'Homo_sapiens' \
            --assembly 'GRCh38' \
            --outdir $OUTDIR \
            --plotGeneDensity false \
            --geneList $GENE_LIST \
            --recolor true \
            --sampleInfo $SAMPLE_INFO \
            --clustering true \
            --clusteringLabel $CLUSTERING_LABEL \
            --clustering_memory '16.GB' \
            --clustering_cpus 4 \
            --clustering_time '20.h' \
            --dropChr '' \
            --tissue_color_palette 'Set3' \
            --no_sampleLabel false \
            --queue $QUEUE \
            --extra_mount $EXTRA_MOUNT \
            -resume
