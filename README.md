# SARS-CoV-2_Database_LS_NLIE

Pipeline to summarize the Results of the IonTorrent S5 SARS-CoV-2 AmpliSeq Panel

# Installation

*  clone the repository:

`git clone https://github.com/nilich/SARS-CoV-2_Database_LS_NLIE.git`

* create a conda environment using the BioHub.yaml files:

`conda env create --file BioHub.yaml`

* Install VADR (https://github.com/ncbi/vadr) locally and set the paths to your installation and the model in the config.yaml

* In the config.yaml:
* set the path to the directory containing the raw reads (RAW_READS) and the consensus sequences(CNS), the output directory (RESULTS_DIR)
* set the parameters VariantCalling and Consensus to "S5" if you want to use the results of the S5 plugins or to "Custom" if you prefer to compute the variant calling using LoFreq and the consensus sequence using ivar.
