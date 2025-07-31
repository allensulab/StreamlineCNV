#!/bin/bash

DATA_DIR="$PWD/data"
OUTDIR="$PWD/results"
GENE_LIST="$PWD/data/GeneListExample"
CLUSTERING_LABEL="$PWD/data/ClusteringLabel"
SAMPLE_INFO="$PWD/data/SampleInfo"

singularity exec -B "$DATA_DIR:/data" docker://ghcr.io/biomicrocenter/streamlinecnv:main nextflow clean -f
singularity exec -B "$DATA_DIR:/data" docker://ghcr.io/biomicrocenter/streamlinecnv:main nextflow run \
                                                                                                StreamlineCNV.nf \
                                                                                                -profile singularity_no_nextflow \
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
