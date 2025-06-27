#!/bin/bash

REF_FILES_DIR="$PWD/data/ref_files"
FASTQ_DATA="$PWD/data/fastq_files"
OUTDIR="$PWD/results"
GENE_LIST="$PWD/data/geneListExample"
SAMPLE_INFO="$PWD/data/SampleInfo"

nextflow clean -f
nextflow run MIT_IGB_CNV.nf \
            -profile local,singularity \
            --ref_files_dir $REF_FILES_DIR \
            --fastq "$FASTQ_DATA/*fq" \
            --species 'Homo_sapiens' \
            --assembly 'GRCh38' \
            --outdir $OUTDIR \
            --plotGeneDensity true \
            --geneList $GENE_LIST \
            --recolor true \
            --sampleInfo $SAMPLE_INFO \
            -resume
