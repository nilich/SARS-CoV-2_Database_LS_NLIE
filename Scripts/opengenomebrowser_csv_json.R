## RScript to automatically generate organism.json and genome.json files for ogb ;)
## libraries müssen in conda umgebung installiert sein!
library(jsonlite)
library(readr)

args = commandArgs(trailingOnly=TRUE)
#args[1] usw.
## read files from cmd line inputfile, sample, pangofile, xslx file, output
sample <- args[1]
print(args[1])
# read default genome.json
genome <- fromJSON(args[2],
                   simplifyVector = FALSE)

genome <- jsonlite:::null_to_na(genome)
# user defined csv file -> read from snakemake input
t <- data.frame(read_csv(args[3], col_names = T))


# read pangolin linage
#"/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/Consensus/pangolin_output.csv"
p <- data.frame(read_csv(args[4]))
##
# do not use for loop, instead read sample from snakemake input

#define output filename
filename=args[5]
#dir="/home/ubuntusnpi/Documents/TestOpenGenomeBrowser"
#filename=paste(dir, "/", sample, ".json", sep ="")

#set json fields
genome$identifier <- as.character(sample)
genome$isolation_date <- t[1,sample]
genome$growth_condition <- t[2,sample]
genome$geographical_name <- t[3,sample]
genome$library_preparation <- t[4,sample]
genome$sequencing_tech <- t[5,sample]
genome$assembly_tool <- t[6,sample]
genome$assembly_date <- format(Sys.time(), "%Y-%m-%d ") # -> get current date
genome$assembly_fasta_file <- paste(sample, "fna", sep=".")
# get filenames
genome$cds_tool_faa_file <- paste(sample, "faa", sep=".")
genome$cds_tool_ffn_file <- paste(sample, "ffn", sep=".")
genome$cds_tool_gbk_file <- paste(sample, "gbk", sep=".")
genome$cds_tool_gff_file <- paste(sample, "gff", sep=".")
#pangolinage as tag
l <- (p[p$taxon == sample, ]$lineage)
genome$tags <- list(l)
# convert to JSON
j <- toJSON(genome, pretty=T, auto_unbox = TRUE, na="null")
# write json file
write(j, filename)
