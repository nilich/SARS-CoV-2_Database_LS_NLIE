### Pipeline for the analysis of SARS-CoV-2 AmpliSeq Data

## TODO: Documentation

import os
import glob
configfile: "config.yaml"

# 1. Initialize Result Directory
RESULTS_DIR = config["general"]["RESULTS_DIR"]
RAW_READS = config["general"]["RAW_READS"]

SAMPLES, = glob_wildcards((RAW_READS + "/{sample}.bam"))
CWD=os.path.abspath(RESULTS_DIR)

## write functions for rule all, see https://github.com/cbg-ethz/V-pipe/blob/sars-cov2/rules/common.smk
## integrate OpengenomeBrowser: see OGB folder

## set variant Calling and Consensus in config file
if config["general"]["VariantCalling"] == "S5":
    include: "Scripts/variantsS5.snake"

if config["general"]["VariantCalling"] == "Custom":
    include: "Scripts/variantCalling.snake"


if config["general"]["Consensus"] == "Custom":
    include: "Scripts/consensusIvar.snake"

if config["general"]["Consensus"] == "S5":
    include: "Scripts/consensusS5.snake"

rule all:
    input:
        expand(RESULTS_DIR + "/VariantAnnotation/{sample}.variants.annot.tab", sample=SAMPLES),
        expand(RESULTS_DIR + "/VariantAnnotation/{sample}.Sprotein.annot.tab", sample=SAMPLES),
        RESULTS_DIR + "/Consensus/pangolin_output.csv",
        RESULTS_DIR + "/Consensus/allSequences.fasta",
        expand(RESULTS_DIR + "/VADR/{sample}/{sample}.vadr.alt.list", sample=SAMPLES),
        CWD + "/Summary.html",

### bam to fastq
rule bamTofastq:
    input:
        RAW_READS + "/{sample}.bam",
    output:
        RAW_READS + "/{sample}.fastq",
    log: "log/bamTofastq_{sample}.log"
    shell:
        """
        samtools bam2fq {input} > {output}
        """

### trimming and Quality control
rule trimmomatic:
    input:
        RAW_READS + "/{sample}.fastq",
    output:
        RESULTS_DIR + "/Trimming/{sample}_trimmed.fq.gz"
    params:
        adapters = config["References"]["adapters"]
    threads: config["general"]["threads"]
    log: "log/trimmomatic_{sample}.log"
    shell:
        """
        trimmomatic SE -threads {threads} -phred33 {input} {output} ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -trimlog {log}
        """

rule fastqc:
    input:
         RESULTS_DIR + "/Trimming/{sample}_trimmed.fq.gz"
    output:
        html = RESULTS_DIR + "/Trimming/{sample}_fastqc.html"
    threads: config["general"]["threads"]
    log: "log/fastqc_{sample}.log"
    shell:
        """
        fastqc {input} -t {threads} 2> {log}
        """

### Alignment
## bwa index is in reference folder, add indexing step?
rule bwa_map:
    input:
        RESULTS_DIR + "/Trimming/{sample}_trimmed.fq.gz"
    output:
        aligned=temp(RESULTS_DIR + "/Alignment/{sample}.sam"),
        sorted=RESULTS_DIR + "/Alignment/{sample}.sorted.bam",
        bai=RESULTS_DIR + "/Alignment/{sample}.sorted.bam.bai"
    params:
        genome=config["References"]["Genome"]
    threads: config["general"]["threads"]
    shell:
        """
        bwa mem -t {threads} {params.genome} {input} > {output.aligned}
        samtools sort -@{threads} -O bam -o {output.sorted} {output.aligned}
        samtools index {output.sorted}
        """

### Mapping statistics and coverage analysis
# Try using pysam: mapped: print reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(filename) ])
#                  unmapped: print reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[3]) for l in pysam.idxstats(filename) ])
### Number of mapped reads
rule samtools_readcount:
    input:
        RESULTS_DIR + "/Alignment/{sample}.sorted.bam"
    output:
        temp(RESULTS_DIR + "/Alignment/{sample}.readcount.txt")
    shell:
        """
        MAP=$(samtools idxstats {input} | head -n1 | awk -F '\t' '{{print $3}}')
        UMAP=$(samtools idxstats {input} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}')
        SAMPLE={wildcards.sample}
        combined="$SAMPLE,$MAP,$UMAP"
        echo $combined >> {output}
        """

rule cat_readcount:
    input:
        expand(RESULTS_DIR + "/Alignment/{sample}.readcount.txt", sample=SAMPLES)
    output:
        temp(RESULTS_DIR + "/Alignment/readcount_all.txt")
    shell:
        """
        cat {input} > {output}
        """

### Coverage
rule samtools_coverage:
    input:
        RESULTS_DIR + "/Alignment/{sample}.sorted.bam"
    output:
        RESULTS_DIR + "/Coverage/{sample}.cov"
    shell:
        """
        samtools depth -d 0 -aa {input} > {output}
        """

rule mean_cov:
    input:
        RESULTS_DIR + "/Coverage/{sample}.cov"
    output:
        temp(RESULTS_DIR + "/Coverage/{sample}.covStats.txt"),
    shell:
        """
        SAMPLE={wildcards.sample}
        MEAN=$(awk '{{ total += $3 }} END {{ print total/NR }}' {input})
        B10=$(awk '$3>=10' {input} | wc -l)
        B10frac=$(awk -v var1=$B10 -v var2=29903 'BEGIN {{ print  ( var1 / var2 *100 ) }}')
        combined="$SAMPLE,$MEAN,$B10,$B10frac"
        echo $combined >> {output}
        """

rule cat_cov:
    input:
        expand(RESULTS_DIR + "/Coverage/{sample}.covStats.txt", sample=SAMPLES),
    output:
        temp(RESULTS_DIR + "/Coverage/covStats_all.txt"),
    shell:
        """
        cat {input} > {output}
        """

## get % Ns per Seq
rule get_N:
    input:
        RESULTS_DIR + "/Consensus/{sample}.fa",
    output:
        temp(RESULTS_DIR + "/Consensus/{sample}.N.txt")
    threads:
        config["general"]["threads"]
    shell:
        """
        python Scripts/get_Ns.py -i {input} -s {wildcards.sample} -o {output}
        """

rule cat_N:
    input:
        expand(RESULTS_DIR + "/Consensus/{sample}.N.txt", sample=SAMPLES),
    output:
        temp(RESULTS_DIR + "/Consensus/Ns_all.txt"),
    shell:
        """
        cat {input} > {output}
        """

### run vadr for Quality Check of assembly
### VADR is installed on the nfs, use this version or install it locally and
### set the path to the vadr installation in the config file
rule vadr:
    input:
        RESULTS_DIR + "/Consensus/{sample}.fa"
    output:
        trimmed=temp(RESULTS_DIR + "/VADR/{sample}.tr.fa"),
        alt=RESULTS_DIR + "/VADR/{sample}/{sample}.vadr.alt.list",
        dir=directory(RESULTS_DIR + "/VADR/{sample}"),
        sum=temp(RESULTS_DIR + "/VADR/{sample}.summary.csv")
    params:
        vadrdir=config["general"]["vadrdir"],
        vadrmdir=config["general"]["vadrmdir"]
    conda:
        "envs/vadr.yaml"
    shell:
        """
        export VADRINSTALLDIR={params.vadrdir}
        export VADRSCRIPTSDIR="$VADRINSTALLDIR/vadr"
        export VADRMODELDIR="$VADRINSTALLDIR/vadr-models-calici"
        export VADRINFERNALDIR="$VADRINSTALLDIR/infernal/binaries"
        export VADREASELDIR="$VADRINSTALLDIR/infernal/binaries"
        export VADRHMMERDIR="$VADRINSTALLDIR/hmmer/binaries"
        export VADRBIOEASELDIR="$VADRINSTALLDIR/Bio-Easel"
        export VADRSEQUIPDIR="$VADRINSTALLDIR/sequip"
        export VADRBLASTDIR="$VADRINSTALLDIR/ncbi-blast/bin"
        export VADRFASTADIR="$VADRINSTALLDIR/fasta/bin"
        export PERL5LIB="$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch"
        export PATH="$VADRSCRIPTSDIR":"$PATH"

        perl $VADRSCRIPTSDIR/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 {input} > {output.trimmed}
        v-annotate.pl -f --split --cpu 8 --glsearch -s --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 \
        --alt_fail lowscore,fstukcnf,insertnn,deletinn --mdir {params.vadrmdir} \
        {output.trimmed} {output.dir}

        if [[ $(grep -vc '#' {output.alt}) > 1 ]]; then echo "{wildcards.sample},FAILED" >> {output.sum} ; else echo "{wildcards.sample},PASSED" >> {output.sum}; fi
        """

## get vadr summary Pass, failed
rule cat_vadr:
    input:
        expand(RESULTS_DIR + "/VADR/{sample}.summary.csv", sample=SAMPLES),
    output:
        summary=temp(RESULTS_DIR + "/VADR/vadr_summary.csv"),
    shell:
        """
        cat {input} > {output}
        """

## cat consensus for pangolin
rule cat_cns:
    input:
        expand(RESULTS_DIR + "/Consensus/{sample}.fa", sample=SAMPLES),
    output:
        RESULTS_DIR + "/Consensus/allSequences.fasta",
    shell:
        """
        cat {input} > {output}
        """

##### run Pangolin for Linage identification
### pangolin ist installed in the BioHub envrionment
### --> update regularly according to https://cov-lineages.org/pangolin_docs/updating.html
### --> for major releases --> update BioHub environment?
rule pangolin:
    input:
        RESULTS_DIR + "/Consensus/allSequences.fasta",
    output:
        full=RESULTS_DIR + "/Consensus/pangolin_output.csv",
        reduced = temp(RESULTS_DIR + "/Consensus/pangolin_output.red.csv"),
    shell:
        """
        pangolin {input} --outfile {output.full}
        sed '1d' {output.full} | cut -f 1,2,5 -d ',' > {output.reduced}
        """

#### Count variants
rule count_var:
    input:
        RESULTS_DIR + "/VariantAnnotation/{sample}.variants.annot.tab",
    output:
        temp(RESULTS_DIR + "/VariantAnnotation/{sample}.countVar.txt"),
    shell:
        """
        SAMPLE={wildcards.sample}
        HOMO=$(awk '$5>=0.5{{c++}} END{{print c+0}}' {input})
        HET=$(awk '$5<0.5{{c++}} END{{print c+0}}' {input})
        FUR=$(cut -f 12,12 {input} | grep -c '68[1-6]')
        combined="$SAMPLE,$HOMO,$HET,$FUR"
        echo $combined >> {output}
        """

rule cat_count:
    input:
        expand(RESULTS_DIR + "/VariantAnnotation/{sample}.countVar.txt", sample=SAMPLES),
    output:
        temp(RESULTS_DIR + "/VariantAnnotation/counts_Var.txt"),
    shell:
        """
        cat {input} > {output}
        """

### create Summary
rule cat_Stats:
    input:
        cov = RESULTS_DIR + "/Coverage/covStats_all.txt",
        rc = RESULTS_DIR + "/Alignment/readcount_all.txt",
        pg = RESULTS_DIR + "/Consensus/pangolin_output.red.csv",
        vadr = RESULTS_DIR + "/VADR/vadr_summary.csv",
        n = RESULTS_DIR + "/Consensus/Ns_all.txt",
        var = RESULTS_DIR + "/VariantAnnotation/counts_Var.txt",
    output:
        out=temp(CWD + "/Summary_Results.csv")
    run:
        import pandas as pd
        import numpy as np
        from functools import reduce
        cov = pd.read_csv(input.cov, names=["Name", "Mean_Coverage", "Coverage_minDP10", "Perc_Coverage_minDP10"])
        rc = pd.read_csv(input.rc, names=["Name", "MappedReads", "TotalReads"])
        pg = pd.read_csv(input.pg, names=["Name", "PangoLinage", "WHOLinage"])
        vadr = pd.read_csv(input.vadr, names=["Name", "Vadr_QC"])
        n = pd.read_csv(input.n, names=["Name", "Perc_N"])
        var = pd.read_csv(input.var, names=["Name", "Variants Homozygot", "Variants Heterozygot", "Mutation_PRRARS"])
        data_frames = [pg, cov, rc, n, vadr, var]
        sum = reduce(lambda  left,right: pd.merge(left,right,on=['Name'], how='outer'), data_frames)
        sum.to_csv(output.out, sep=',', encoding='utf-8', index = False, header=True, na_rep='NA')

## add link to variant table in Summary?
rule summary:
    input:
        CWD + "/Summary_Results.csv"
    output:
        CWD + "/Summary.html",
    params:
        wd=CWD
    shell:
        """
        Rscript -e "rmarkdown::render('Scripts/summary.Rmd', params=list(input = '{input}', wd='{params.wd}'),  output_file = '{output}')"
        """
