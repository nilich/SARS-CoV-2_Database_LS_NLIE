#!/bin/bash
#SBATCH --mem=80g
#SBATCH --cpus-per-task=32
#SBATCH --output=snakemake-%A.out
#SBATCH --error=snakemake-%A.err
#SBATCH --job-name=BioHub

conda activate BioHub

snakemake -j 32 --use-conda
