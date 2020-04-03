############################################
# Compute pileup frequencies for multiple samples, creating fingerprinting DB
# mcadosch 10/17
# Gwen Miller edited 12/2019
############################################

rm(list=ls())
suppressMessages(library(dplyr))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(methods))

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
parser$add_argument('--prev_fng_db_path', type="character", nargs="+", 
                    help='previous FNG database to merge new FNG data into', required=T)

args <- parser$parse_args()
files <- args$pileup_pct_count_paths
samples_tsca_ids <- args$samples_tsca_ids
prev_fng_db <- args$prev_fng_db_path

##########
## Compile table for all new samples
##########
# Number of samples
N = length(files)
samples_data <- vector("list", N)

# Compile all samples
for (i in seq_along(files)) {
	sample_tsca_id <- samples_tsca_ids[i]
	sample_data <- read.delim(files[i]) %>% mutate(batch=sample_tsca_id)
	samples_data[[i]] <- sample_data
}

new_fngs <- samples_data %>% bind_rows()

##########
## Merge new FNG db with the previously existing FNG db
##########
# read in previous FNG db
prev_fng_db <- read.delim(prev_fng_db)

# stop running if they do no not share the same column names (order doesn't matter)
if (setequal(colnames(prev_fng_db), colnames(new_fngs))){
  stop("Error: The old and new fingerprinting databases do not have the same columns")
}

# merge and delete duplicate rows
no_dups_merged <- bind_rows(new_fngs, prev_fng_db) %>% distinct(merged)

##########
## write out the final compiled FNG db
##########
output = paste0("fingerprinting_db.txt")
write.delim(no_dups_merged, output)