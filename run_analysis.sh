#!/bin/bash
#SBATCH --mem=50g
#SBATCH --cpus-per-task=24
#SBATCH --output=snakemake-%A.out
#SBATCH --error=snakemake-%A.err
#SBATCH --job-name=BioHub

conda activate BioHub_test

snakemake -j 24 --use-conda
