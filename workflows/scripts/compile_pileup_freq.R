############################################
# Compute pileup frequencies for a given sample
# mcadosch 10/17
############################################

rm(list=ls())
suppressMessages(library(Rsamtools))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressPackageStartupMessages(library("argparse"))

# source('/TSCA/helpers.R')

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

##################################################################
#### PERFORM PILEUP 
##################################################################
read_bam_and_run_pileup <- function(bam_path, sample_id, interval_list.granges, max_depth=12000, distinguish_strand=F) {
	# read in bam file with .bai index file instead of .bam.bai (Rsamtools default), because picard indexes into .bai
	bam_file <- BamFile(bam_path, index=gsub('\\.bam$', '.bai', bam_path))
	## Set up pileup parameters
	sbp_param <- ScanBamParam(which=interval_list.granges)
	pu_param <- PileupParam(max_depth=max_depth, distinguish_strand=distinguish_strand)
	## Run pileup
	interval_list.granges.df <- as.data.frame(interval_list.granges) %>% unfactorize() %>% mutate(pos=start)
	res_pileup <- 
		pileup(bam_file, scanBamParam=sbp_param, pileupParam=pu_param) %>%
		merge(interval_list.granges.df, all.x=T, all.y=T) %>%
		mutate(
			count = ifelse(is.na(count), 0, count),
			sample = sample_id,
			snp_id=paste0('chr', seqnames, ':', pos)
			)
	return(res_pileup)
}

##################################################################
#### COMPILE PILEUP FREQUENCY
##################################################################
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

##################################################################
#### MAIN
##################################################################
parser <- ArgumentParser(description='Compile pileup frequencies for a given sample')
parser$add_argument('--sample_id', metavar='sample_id', type="character",
                   help='sample_id', required=T)
parser$add_argument('--bam_path', metavar='bam_path', type="character",
                   help='BAM PATH', required=T)
parser$add_argument('--fluidgm_snp_path', type="character",
                   help='Path to fluidgm_snp database', required=T)

args <- parser$parse_args()
sample_id <- args$sample_id
bam_path <- args$bam_path
fluidgm_snp_path <- args$fluidgm_snp_path

fluidgm_snp.raw <- read.maf(fluidgm_snp_path)
fluidgm_snp <- fluidgm_snp.raw[, c('CHROM', 'POSITION', 'SNP', 'SNP_ALLELES')] %>% unique()
fluidgm_snp.granges <- makeGRangesFromDataFrame(fluidgm_snp, seqnames.field='CHROM',  start.field='POSITION', end.field='POSITION', keep.extra.columns=T)

cat('Performing pileup for sample', sample_id, '...')
res_pileup <- read_bam_and_run_pileup(bam_path, sample_id, fluidgm_snp.granges)

res_pileup_with_allele <- 	
	res_pileup %>% 
	group_by(sample, seqnames, pos) %>%
	mutate(
		pct_count=count/sum(count), 
		pct_count=ifelse(is.na(pct_count), 0, pct_count),
		major_allele=sort(nucleotide[which(pct_count==max(pct_count))])[1],
		minor_allele=select_minor_allele(nucleotide, pct_count, major_allele),
		total_cvg = sum(count),
		genotype=paste0(major_allele, minor_allele) #Edit: Sahar make space between alleles
		) %>% 
	ungroup() %>%
	mutate(genotype=ifelse(genotype=='NANA', NA, genotype))

write.delim(res_pileup_with_allele, paste0(sample_id, '_fng.pileup.txt'))

res_pileup_freq <- dcast(res_pileup_with_allele, sample + snp_id + seqnames + pos + major_allele + minor_allele + total_cvg + genotype ~ nucleotide, value.var='count') %>% unfactorize()
res_pileup_pctCount <- dcast(res_pileup_with_allele, sample + snp_id + seqnames + pos + major_allele + minor_allele + total_cvg +genotype ~ nucleotide, value.var='pct_count') %>% unfactorize()

write.delim(res_pileup_pctCount, paste0(sample_id, '_fng.pileup_pct_count.txt'))

