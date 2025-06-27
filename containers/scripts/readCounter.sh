#!/bin/sh
#SBATCH -N 1 #request Bourne shell as shell for job
#SBATCH -n 8
#SBATCH --exclude=b[14-15]
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate cnv
readCounter -w 500000 $1 > input.wig
