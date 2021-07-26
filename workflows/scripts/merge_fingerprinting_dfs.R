############################################
# Merge fingerprinting databases together
# Goal: maintain a database of all FNG data in one table.
# Gwen Miller May 2021
############################################
suppressMessages(library(readr))
suppressMessages(library(argparse))
# install.packages("argparse")


##################################################################
#### Function to merge FNG data into one table
##################################################################
merge_fng_dbs <- function(path_old, path_new){
  twist_old <- readr::read_tsv(file = path_old)
  twist_new <- readr::read_tsv(file = path_new)
  
  cols_of_interest <- colnames(twist_new)[!(colnames(twist_new) %in% c('X.', 'NA.'))]
  twist_old <- twist_old[, colnames(twist_old) %in% cols_of_interest]
  twist_new <- twist_new[, colnames(twist_new) %in% cols_of_interest]
  
  # do they share the same columns?
  print("Checking if the new table contains all the columns found in the old table...")
  if (length(setdiff(colnames(twist_old), colnames(twist_new))) != 0) {
    stop("The columns don't match between the old and new fingerprinting file. Stopping.")
  }
  
  # dim
  print("Checking the table dimensions...")
  print(paste0("dim(twist_old):",dim(twist_old)))
  print(paste0("dim(twist_new):",dim(twist_new)))
  
  # merge and delete duplicate rows
  print("Merge FNG tables and delete duplicate rows...")
  merged = bind_rows(twist_new, twist_old)
  no_dups_merged = distinct(merged)
  print(paste0("dim(merged): ", dim(merged)))
  print(paste0("dim(no_dups_merged): ", dim(no_dups_merged)))
  
  return(no_dups_merged)
}



##################################################################
#### MAIN
##################################################################
parser <- ArgumentParser(description='Merge the fingerprinting databases all the way through a given batch (sample set)')
parser$add_argument('--sample_set_id', metavar='sample_set_id', type="character",
                    help='sample_set_id', required=T)
parser$add_argument('--old_fng_db', metavar='bam_path', type="character",
                    help='Path to old FNG db', required=T)
parser$add_argument('--new_fng_db', type="character",
                    help='Path to new FNG db', required=T)

args <- parser$parse_args()
sample_set_id <- args$sample_set_id
old_fng_db <- args$old_fng_db
new_fng_db <- args$new_fng_dbs

merged_fng_df <- merge_fng_dbs(old_fng_db, new_fng_db)
write.table(merged_fng_df, 
          file = paste0('fingerprinting_db_through_', sample_set_id, '.txt'), 
          sep='\t', row.names = FALSE, quote=FALSE)
