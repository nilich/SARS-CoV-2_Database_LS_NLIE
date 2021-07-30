#!/bin/bash
#SBATCH --mem=80g
#SBATCH --cpus-per-task=32
#SBATCH --output=snakemake-%A.out
#SBATCH --error=snakemake-%A.err
#SBATCH --job-name=BioHub

conda activate BioHub

# rename files to remove barcode:
rename 's/(.*)IonCode_\d+.(bam)/\1\2/' Reads/*.bam
rename 's/(.*)IonCode_\d+.(vcf)/\1\2/' Reads/*.vcf

snakemake -j 32 --use-conda
