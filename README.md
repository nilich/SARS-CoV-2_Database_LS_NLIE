# SARS-CoV-2_Database_LS_NLIE

Pipeline to summarize the Results of the IonTorrent S5 SARS-CoV-2 AmpliSeq Panel

## Installation

1. Clone the repository:

  `git clone https://github.com/nilich/SARS-CoV-2_Database_LS_NLIE.git`

2. create a conda environment using the BioHub.yaml files:

  `conda env create --file BioHub.yaml`

3. Install VADR (https://github.com/ncbi/vadr) locally and set the paths to your installation and the model in the config.yaml:
    ```
    vadrdir: "/mnt/nfs/bio/software/QC/VADR/" # path to your local vadr installation`
    vadrmdir: "/mnt/nfs/bio/software/QC/VADR/vadr-models-sarscov2-1.2-2" # vadr model dir'
    ```

4. Configure your pipeline:
  * set the path to the directory containing the raw reads (RAW_READS) and the consensus sequences (CNS), the output directory (RESULTS_DIR) in the config.yaml
    ```
    RAW_READS: "/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Reads"
    CNS_DIR: "/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Consensus"
    RESULTS_DIR: "Results_S5"
    ```
    Note: the directory containing consensus sequences is only needed when working with S5 plugins

  * set the parameter VariantCalling in the config.yaml to "S5" for using S5 computed vcf files or to "Custom" for variant calling using LoFreq
    `VariantCalling: "S5" # "S5 or Custom"`

  * set the parameter Consensus in the config.yaml to "S5" for using S5 computed consensus sequences or to "Custom" for building consensus sequences using Ivar
    `Consensus: "S5" # "S5 or Custom"`

5. Run the Pipeline:
  * to run the pipeline on the cluster submit the run_analysis.sh script using:
  `sbatch run_analysis.sh`

  * to run the pipeline on your local computer use:
  `snakemake -j {threads}`

    Note: remove the barcode in the read file name before starting the pipeline (already integrated in the run_analysis.sh script) using:
    ```
    rename -n 's/(.*)IonCode_\d+.(bam)/\1\2/' *.bam
    rename -n 's/(.*)IonCode_\d+.(vcf)/\1\2/' *.vcf
    ```
6. Update Pangolin within the BioHub environemnt regularly in order to analyse newest linages:
   ```
   pangolin --update
   ```
   For major releases see: https://cov-lineages.org/resources/pangolin/updating.html

## Output
The output of the pipeline is structured as followed:
```
Results
|-- Alignment
|   |-- sample.readcount.txt
|   |-- sample.sorted.bam --> aligned reads
|   |-- sample.sorted.bam.bai
|-- Consensus
|   |-- allSequences.fasta --> Fasta with all sequences
|   |-- sample.fa --> sample consensus fasta
|   `-- pangolin_output.csv --> full pangolin output
|-- Coverage
|   `-- sample.cov --> Coverage tabke
|-- LoFreq
|   |-- sample.lofreq.DP400.AF003.vcf --> Variant file, only computed if VariantCalling is set to "Custom"
|   |-- sample.lofreq.vcf.gz.tbi
|   |-- sample.sorted.indelqual.bam
|   `-- sample.sorted.indelqual.bam.bai
|-- Trimming
|   `-- sample_trimmed.fq.gz
|-- VADR
|   |-- sample
|   |   |-- sample.vadr.alc --> Summary of the VADR quality check, consult if the quality check fails
|   |   |-- sample.vadr.alt
|   |   |-- sample.vadr.alt.list
|   |    `-- ...
|-- VariantAnnotation
|    |-- sample.Sprotein.annot.tab
|    `-- sample.variants.annot.tab --> List of variants and their effect
`-- Summary.html --> Summary of the results
```
Note: If the VADR Quality check fails, consult the VADR summary and visualize aligned reads (Alignment/*.bam) and vcf files in IGV (https://software.broadinstitute.org/software/igv/) to verify coverage and alignment.
