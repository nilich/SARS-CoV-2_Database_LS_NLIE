# SARS-CoV-2_Database_LS_NLIE

Pipeline to summarize the Results of the IonTorrent S5 SARS-CoV-2 AmpliSeq Panel

# Installation

1. Clone the repository:

`git clone https://github.com/nilich/SARS-CoV-2_Database_LS_NLIE.git`

2. create a conda environment using the BioHub.yaml files:

`conda env create --file BioHub.yaml`

3. Install VADR (https://github.com/ncbi/vadr) locally and set the paths to your installation and the model in the config.yaml
  * `vadrdir: "/mnt/nfs/bio/software/QC/VADR/" # path to your local vadr installation
    vadrmdir: "/mnt/nfs/bio/software/QC/VADR/vadr-models-sarscov2-1.2-2" # vadr model dir'`

4. Configure your pipeline
  * set the path to the directory containing the raw reads (RAW_READS) and the consensus sequences(CNS), the output directory (RESULTS_DIR) in the config.yaml
  `general:
    RAW_READS: "/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Reads"
    CNS_DIR: "/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Consensus"
    RESULTS_DIR: "Results_S5"`

  * set the parameters VariantCalling and Consensus to "S5" if you want to use the results of the S5 plugins or to "Custom" if you prefer to compute the variant calling using LoFreq and the consensus sequence using ivar.
