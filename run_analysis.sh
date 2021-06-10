#!/bin/bash
#SBATCH --mem=50g
#SBATCH --cpus-per-task=24
#SBATCH --output=snakemake-%A.out
#SBATCH --error=snakemake-%A.err
#SBATCH --job-name=Ampli_snakemake

conda activate SARS-Airolo

snakemake -j 24 --use-conda
