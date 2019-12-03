############################################
# Compute pileup frequencies for multiple samples, creating fingerprinting DB
# mcadosch 10/17
############################################

rm(list=ls())
suppressMessages(library(dplyr))
suppressPackageStartupMessages(library("argparse"))

##################################################################
#### Helper functions
##################################################################
write.delim <- function(data, file='', coln=T, rown=F, append=F, ...){
    write.table(data, file, col.names=coln, sep="\t", row.names=rown, quote=F, append=append)
}

##################################################################
#### MAIN
##################################################################
parser <- ArgumentParser(description='Create fingerprinting DB by compiling pileup frequencies for all samples')
parser$add_argument('--pileup_pct_count_paths', type="character", nargs="+",
                   help='pileup_pct_count files', required=T)
parser$add_argument('--samples_tsca_ids', type="character", nargs="+", 
					help='sample tsca_id', required=T)

args <- parser$parse_args()
files <- args$pileup_pct_count_paths
samples_tsca_ids <- args$samples_tsca_ids

# Number of samples
N = length(files)
samples_data <- vector("list", N)

# Compile all samples
for (i in seq_along(files)) {
	sample_tsca_id <- samples_tsca_ids[i]
	sample_data <- read.delim(files[i]) %>% mutate(batch=sample_tsca_id)
	samples_data[[i]] <- sample_data
}

my_db <- samples_data %>% bind_rows()
output = paste0("fingerprinting_db.txt")
write.delim(my_db, output)