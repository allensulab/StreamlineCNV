#!/bin/bash

REF_FILES_DIR="$PWD/data/ref_files"
FASTQ_DATA="$PWD/data/fastq_files"
OUTDIR="$PWD/results"
GENE_LIST="$PWD/data/GeneListExample"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
SAMPLE_INFO="$PWD/data/SampleInfo"

singularity exec -B "$REF_FILES_DIR:/data" docker://ghcr.io/biomicrocenter/streamlinecnv:main nextflow clean -f
singularity exec -B "$REF_FILES_DIR:/data" docker://ghcr.io/biomicrocenter/streamlinecnv:main nextflow run \
                                                                                                StreamlineCNV.nf \
                                                                                                -profile singularity_no_nextflow \
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
                                                                                                -resume
