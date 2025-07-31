#!/bin/bash

DATA_DIR="$PWD/data/"
OUTDIR="$PWD/results"
GENE_LIST="$PWD/data/geneListExample"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
SAMPLE_INFO="$PWD/data/SampleInfo"

nextflow clean -f
nextflow run StreamlineCNV.nf \
            -profile local,singularity \
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
            -resume
