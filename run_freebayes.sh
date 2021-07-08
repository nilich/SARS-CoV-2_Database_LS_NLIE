#!/bin/bash
#SBATCH --array=1-1
#SBATCH --mem=50g
#SBATCH --cpus-per-task=1
#SBATCH --output=snakemake-%A.out
#SBATCH --error=snakemake-%A.err
#SBATCH --job-name=Ampli_snakemake


array=( mock "3015_P1")
export i=${array["$SLURM_ARRAY_TASK_ID"]}

conda activate SARS-Airolo

ALIGNMENT=/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results_old/Alignment
REF=/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/References/NC_045512.2.fasta


#freebayes -f $REF "$ALIGNMENT"/"$i".sorted.bam -m 20 -q 20 -F 0.01 --min-coverage 10 -V > test_freebayes/"$i".vcf

#freebayes -f $REF "$ALIGNMENT"/"$i".sorted.bam -F 0.01 -C 1 > test_freebayes/"$i".F003.vcf

samtools mpileup -f $REF "$ALIGNMENT"/"$i".sorted.bam | varscan pileup2cns --variants --output-vcf 1 > test_freebayes/"$i".varscan.vcf 

#samtools mpileup -f $REF "$ALIGNMENT"/"$i".sorted.bam | varscan mpileup2cns > test_freebayes/"$i".varscan.vcf
#samtools mpileup -aa -A -d 0 -B -f $REF "$ALIGNMENT"/"$i".sorted.bam | ivar variants -m 100 0.1 -p test_freebayes/"$i"_ivar_01 -r $REF 

#test snippy parameters
#freebayes -p 2 -P 0 -C 5 --min-repeat-entropy 1.5 --strict-vcf -q 13 -m 60 --min-coverage 5 -F 0.05 -f $REF "$ALIGNMENT"/"$i".sorted.bam > test_freebayes/"$i".snippyfreebayes.vcf
