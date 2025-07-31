#!/bin/sh
#SBATCH -N 1 #request Bourne shell as shell for job
#SBATCH -n 8
#SBATCH --exclude=b[14-15]
readCounter -w 500000 $1 > input.wig
