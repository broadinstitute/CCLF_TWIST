suppressMessages(library(plyr))
suppressMessages(library(dplyr))

##################################################################
#### Helper functions
##################################################################
read.maf <- function(maf_path, ...){
    if(is.null(maf_path)){
        cat("No maf found.\n")
        return(NA)
    } else if (maf_path=='' | is.na(maf_path)){
        cat("No maf found.\n")
        return(NA)
    }
    if(grepl('/xchip/cga/gdac-prod/', maf_path)){
        cat("maf file not in FUSE format.\n")
        return('maf_not_in_FUSE')
    } else if(!file.exists(maf_path)){
        cat("maf path outdated or file not found:",maf_path, "\n")
        return(NA)
    }
    myMAF <- read.delim(maf_path, stringsAsFactors=F, comment.char="#")
    return(myMAF)
}

unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

select_minor_allele <- function(nucleotide, pct_count, major_allele){
    ## set 40% as approximation for 50%, but please EDIT this in the future, this is a hard cutoff, therefore not the most accurate. Better way might be to model distribution, LOD scores, find the likelihood of a value being within the 50% range......
    cutoff_for_minor_allele <- 0.4  


    if(length(unique(nucleotide))==1){
        minor_allele <- nucleotide
    } else {
        second_highest_pct_count <- sort(pct_count, decreasing=T)[2]
        if(second_highest_pct_count > cutoff_for_minor_allele){
            minor_allele <- nucleotide[which(pct_count==second_highest_pct_count & nucleotide!=unique(major_allele))] ## added nucleotide!=major_allele in the case of 2 alleles having equal pct_count (otherwise will return 2 nucleotides as minor allele, instead of 1 nucleotide)
        } else {
            minor_allele <- nucleotide[which(pct_count==max(pct_count))]
        }
    }

    return(minor_allele)
}

write.delim <- function(data, file='', coln=T, rown=F, append=F, ...){
    write.table(data, file, col.names=coln, sep="\t", row.names=rown, quote=F, append=append)
}