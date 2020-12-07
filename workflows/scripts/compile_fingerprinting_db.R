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
print("Compiling all the new samples into one table...")

# Number of samples
N = length(files)
samples_data <- vector("list", N)

# Compile all samples
for (i in seq_along(files)) {
	sample_tsca_id <- samples_tsca_ids[i]
	sample_data <- read.delim(files[i], colClasses = "character") %>% mutate(batch=sample_tsca_id)
	samples_data[[i]] <- sample_data
}

new_fngs <- samples_data %>% bind_rows()

##########
## Merge new FNG db with the previously existing FNG db
##########
# read in previous FNG dbs
prev_fng_db <- read.delim(prev_fng_db, colClasses = "character")

# stop running if not all of the columns in the old FNG database appear in the new database (order doesn't matter)
cols_of_interest <- colnames(prev_fng_db)[!(colnames(prev_fng_db) %in% c('X.', 'NA.'))]
if (!all(cols_of_interest %in% colnames(new_fngs))){
  cols_missing_in_new = setdiff(cols_of_interest, colnames(new_fngs))
  print("The new fingerprinting databases is missing some columns.")
  print("The missing columns are:")
  print(cols_missing_in_new)
  stop("The new fingerprinting databases is missing some columns")
}


# merge and delete duplicate rows
print("Merging the old and newly created FNG databases, keeping only columns that existed in old database...")
restricted_new_fng_db <- new_fngs[, colnames(new_fngs) %in% cols_of_interest]
restricted_old_fng_db  <- prev_fng_db[, colnames(prev_fng_db) %in% cols_of_interest]
no_dups_merged <- bind_rows(restricted_new_fng_db, restricted_old_fng_db) %>% distinct()

##########
## write out the final compiled FNG db
##########
output = paste0("fingerprinting_db.txt")
write.delim(no_dups_merged, output)
print("Done merging and writing output fingerprinting table.")

