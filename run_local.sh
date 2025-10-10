#!/bin/bash

DATA_DIR="$PWD/data/"
FASTQ_DATA="$PWD/data/fastq_files"
OUTDIR="results"
GENE_LIST="$PWD/data/GeneListExample"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
SAMPLE_INFO="$PWD/data/SampleInfo"
CONTAINER_PROFILE="singularity"

nextflow run StreamlineCNV.nf \
            -profile local,$CONTAINER_PROFILE \
            --data_dir $DATA_DIR \
            --fastq "$FASTQ_DATA/*.fastq" \
            --species 'Homo_sapiens' \
            --assembly 'GRCh38' \
            --outdir $OUTDIR \
            --plotGeneDensity true \
            --geneList $GENE_LIST \
            --recolor true \
            --sampleInfo $SAMPLE_INFO \
            --clustering true \
            --clusteringLabel $CLUSTERING_LABEL \
            -resume
