### Pipeline for the analysis of SARS-CoV-2 AmpliSeq Data

## TODO: Documentation
## TODO: "transportable using conda envs"
## TODO: separate different analysis steps
## TODO: test if triming is working with adapters as params
configfile: "config.yaml"

# 1. Initialize Result Directory
RESULTS_DIR = config["general"]["RESULTS_DIR"]
RAW_READS = config["general"]["RAW_READS"]

SAMPLES, = glob_wildcards((RAW_READS + "/{sample}.bam"))
CWD=os.getcwd()

##Modularisation of workflow: ogb + variant calling
## write functions for rule all, see https://github.com/cbg-ethz/V-pipe/blob/sars-cov2/rules/common.smk
include: "/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/ogb_snake"

## defines which files to output...
## TODO: add outfiles from ogb
rule all:
    input:
        #RESULTS_DIR + "/Alignment/readcount_all.txt",
        #RESULTS_DIR + "/Coverage/covStats_all.txt",
        #expand(RESULTS_DIR + "/Coverage/{sample}.cov", sample=SAMPLES),
        #expand(RESULTS_DIR + "/LoFreq/{sample}.filtered.Sprot.vcf", sample=SAMPLES), # do we need this file? pobably not
        expand(RESULTS_DIR + "/Var_annot/{sample}.lofreq.DP400.AF003.annot.vcf", sample=SAMPLES), # DP and AF filtered vcf, keep
        expand(RESULTS_DIR + "/Var_annot/{sample}.variants.DP400.AF003.annot.tab", sample=SAMPLES),
        expand(RESULTS_DIR + "/Var_annot/{sample}.Sprot.DP400.AF003.annot.tab", sample=SAMPLES),
        RESULTS_DIR + "/Consensus/pangolin_output.csv",
        RESULTS_DIR + "/Consensus/all_samples.fa",
        RESULTS_DIR + "/Summary_Results.csv",
        #expand(RESULTS_DIR + "/Prokka/{sample}/{sample}.txt", sample=SAMPLES),
        expand(RESULTS_DIR + "/Var_annot/{sample}.md", sample=SAMPLES),
        expand(RESULTS_DIR + "/QC/VADR/{sample}/{sample}.vadr.alt.list", sample=SAMPLES),
        RESULTS_DIR + "/QC/VADR/vadr_summary.csv",
        CWD + "/" + RESULTS_DIR + "/summary.html",
        expand(RESULTS_DIR + "/ogb/{sample}.json", sample=SAMPLES),

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
## TODO: test if triming is working with adapters as params
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
rule bwa_map:
    input:
        RESULTS_DIR + "/Trimming/{sample}_trimmed.fq.gz"
    output:
        aligned=temp(RESULTS_DIR + "/Alignment/{sample}.sam"),
        sorted=RESULTS_DIR + "/Alignment/{sample}.sorted.bam"
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
        combined="$SAMPLE,$MEAN,$B10,$B10frac
        echo $combined >> {output}
        """
        #B400=$(awk '$3>=400' {input} | wc -l)
        #B400frac=$(awk -v var1=$B400 -v var2=29903 'BEGIN {{ print  ( var1 / var2 *100 ) }}')
        #combined="$SAMPLE,$MEAN,$B10,$B10frac,$B400,$B400frac"

rule cat_cov:
    input:
        expand(RESULTS_DIR + "/Coverage/{sample}.covStats.txt", sample=SAMPLES),
    output:
        temp(RESULTS_DIR + "/Coverage/covStats_all.txt"),
    shell:
        """
        cat {input} > {output}
        """

### create consensus --> use samtools as long as ivar did not update Manual (https://andersen-lab.github.io/ivar/html/manualpage.html)
rule consensus:
    input:
        RESULTS_DIR + "/Alignment/{sample}.sorted.bam"
    output:
        RESULTS_DIR + "/Consensus/{sample}.fa"
    threads:
        config["general"]["threads"]
    shell:
        """
        samtools mpileup -d 0 -A {input} | ivar consensus -p {output} -i {wildcards.sample} -q 20 -m 10 -n N
        """
### run vadr for Quality Check of assembly
### you need to install the VADR locally and set the path to the vadr installation in the config file
rule vadr:
    input:
        RESULTS_DIR + "/Consensus/{sample}.fa"
    output:
        trimmed=temp(RESULTS_DIR + "/QC/VADR/{sample}.tr.fa"),
        alt=RESULTS_DIR + "/QC/VADR/{sample}/{sample}.vadr.alt.list",
        dir=directory(RESULTS_DIR + "/QC/VADR/{sample}"),
        sum=RESULTS_DIR + "/QC/VADR/{sample}.summary.csv"
    params:
        vadrdir=config["general"]["vadrdir"]
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
        v-annotate.pl -f --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 \
        --alt_fail lowscore,fstukcnf,insertnn,deletinn --mdir /mnt/nfs/bio/software/QC/VADR/vadr-models-sarscov2-1.2-2 \
        {output.trimmed} {output.dir}

        if [[ $(grep -vc '#' {output.alt}) > 1 ]]; then echo "{wildcards.sample},FAILED" >> {output.sum} ; else echo "{wildcards.sample},PASSED" >> {output.sum}; fi
        """
## get vadr summary Pass, failed
rule cat_vadr:
    input:
        expand(RESULTS_DIR + "/QC/VADR/{sample}.summary.csv", sample=SAMPLES),
    output:
        summary=temp(RESULTS_DIR + "/QC/VADR/vadr_summary.csv"),
    shell:
        """
        cat {input} > {output}
        """

## cat consensus for pangolin
rule cat_cns:
    input:
        expand(RESULTS_DIR + "/Consensus/{sample}.fa", sample=SAMPLES),
    output:
        temp(RESULTS_DIR + "/Consensus/all_samples.fa"),
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
        RESULTS_DIR + "/Consensus/all_samples.fa",
    output:
        full=RESULTS_DIR + "/Consensus/pangolin_output.csv",
        reduced = RESULTS_DIR + "/Consensus/pangolin_output.red.csv",
    shell:
        """
        pangolin {input} --outfile {output.full}
        sed '1d' {output.full} | cut -f 1,2,5 -d ',' > {output.reduced}
        """

### create Summary
rule cat_Stats:
    input:
        cov = RESULTS_DIR + "/Coverage/covStats_all.txt",
        rc = RESULTS_DIR + "/Alignment/readcount_all.txt",
        pg = RESULTS_DIR + "/Consensus/pangolin_output.red.csv",
        vadr = RESULTS_DIR + "/QC/VADR/vadr_summary.csv"
    output:
        out=RESULTS_DIR + "/Summary_Results.csv"
    run:
        import pandas as pd
        import numpy as np
        from functools import reduce
        cov = pd.read_csv(input.cov, names=["Name", "Mean_Coverage", "Breath_10", "Perc_Breath_10"])
        rc = pd.read_csv(input.rc, names=["Name", "MappedReads", "TotalReads"])
        pg = pd.read_csv(input.pg, names=["Name", "Linage", "scorpio_call"])
        vadr = pd.read_csv(input.vadr, names=["Name", "Vadr_QC"])
        data_frames = [pg, cov, rc, vadr]
        sum = reduce(lambda  left,right: pd.merge(left,right,on=['Name'], how='outer'), data_frames)
        sum.to_csv(output.out, sep=',', encoding='utf-8', index = False, header=True, na_rep='NA')

## TODO: test if setting of the WorkingDir is working!!
rule summary:
    input:
        RESULTS_DIR + "/Summary_Results.csv",
    output:
        CWD + "/" + RESULTS_DIR + "/summary.html",
    params:
        wd=CWD
    shell:
        """
        Rscript -e "rmarkdown::render('Scripts/summary.Rmd', params=list(input = '{input}', dir='{params.wd}'),  output_file = '{output}')"
        """

################## Variant Calling
### TODO: this part is optional, user can also provide Variantcalling of S5!
## Generalize filenames
rule indelqual:
    input:
        RESULTS_DIR + "/Alignment/{sample}.sorted.bam"
    output:
        temp(RESULTS_DIR + "/LoFreq/{sample}.sorted.indelqual.bam"),
        temp(RESULTS_DIR + "/LoFreq/{sample}.sorted.indelqual.bam.bai")
    threads:
        config["general"]["threads"]
    shell:
        """
        lofreq indelqual -u 4 -o {output} {input}
        samtools index {output}
        """

rule lofreq:
    input:
        RESULTS_DIR + "/LoFreq/{sample}.sorted.indelqual.bam"
    output:
        RESULTS_DIR + "/LoFreq/{sample}.lofreq.vcf"
    threads:
        config["general"]["threads"]
    params:
        ref=config["References"]["Genome"]
    shell:
        """
        lofreq call-parallel --pp-threads {threads} --call-indels -f {params.ref} -o {output} {input}
        """

rule bgzip_tabix:
    input:
        RESULTS_DIR + "/LoFreq/{sample}.lofreq.vcf"
    output:
        temp(RESULTS_DIR + "/LoFreq/{sample}.lofreq.vcf.gz"),
    shell:
        """
        bgzip {input}
        tabix {output}
        """

rule filter_vcf:
    input:
        RESULTS_DIR + "/LoFreq/{sample}.lofreq.vcf.gz"
    output:
        all=RESULTS_DIR + "/LoFreq/{sample}.lofreq.DP400.AF003.vcf",
        #sprot=RESULTS_DIR + "/LoFreq/{sample}.filtered.Sprot.vcf"
    shell:
        """
        bcftools filter -i "DP>400 & AF>0.03" {input} -o {output.all}
        """
        #grep -v '#' {output.all} | awk '$2>21562 && $2<25385' > {output.sprot}

rule annot_vcf:
    input:
        RESULTS_DIR + "/LoFreq/{sample}.lofreq.DP400.AF003.vcf"
    output:
        all = temp(RESULTS_DIR + "/Var_annot/{sample}.variants.annot.vcf"),
        table = RESULTS_DIR + "/Var_annot/{sample}.variants.annot.tab",
        sprot = RESULTS_DIR + "/Var_annot/{sample}.Sprotein.annot.tab"
    shell:
        """
        snpEff ann NC_045512.2 -noStats {input} > {output.all}
        cat {output.all} | Scripts/vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT AF DP "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_P" | sort -nu -k 2,2 > {output.table}
        sed -i 's/CHROM.*/CHROM\tPOS\tREF\tALT\tAF\tDP\tEFFECT\tGENE\tAA_CHANGE/' {output.table}
        cat {output.table} | awk '$2>21562 && $2<25385' > {output.sprot}
        """
        #snpEff ann NC_045512.2 {input} > {output.all}
        #header="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tVariant\tEffect\tProtein\tReplacement"
        #echo -e $header >> {output.table}
        #grep -v '#' {output.all} | cut -f 1,2,3,4,11 -d '|' | sed 's/ANN=[A-Z]//' | \
        #sed 's/p\.//' | sed 's/|/\t/g' >> {output.table}
        #cat {output.table} | awk '$2>21562 && $2<25385' > {output.sprot}
