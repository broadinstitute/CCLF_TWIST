rm(list=ls())
suppressMessages(library(Rsamtools))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressPackageStartupMessages(library("argparse"))

##################################################################
#### HELPER FUNCTIONS
##################################################################

read.maf <- function(maf_path, ...){
    if(is.null(maf_path)){
        cat("No maf found.\n")
        return(NA)
    } else if (maf_path=='' | is.na(maf_path)){
        cat("No maf found.\n")
        return(NA)
    }
    myMAF <- read.delim(maf_path, stringsAsFactors=F, comment.char="#")
    return(myMAF)
}

unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

read_bam_and_run_pileup <- function(bam_path, sample_id, interval_list.granges, max_depth=12000, distinguish_strand=F) {
	
	cat('Performing pileup for sample', sample_id, '\n')
	# read in bam file with .bai index file instead of .bam.bai (Rsamtools default), because picard indexes into .bai
	bam_file <- BamFile(bam_path, index=gsub('\\.bam$', '.bai', bam_path))

	##################################################################
	#### 2. PERFORM PILEUP 
	##################################################################
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
			snp_id=paste0('chr', seqnames, ':', pos, Reference_Allele, '>', Tumor_Seq_Allele2)		
			)
	return(res_pileup)
}

select_minor_allele <- function(nucleotide, pct_count, major_allele){
	cutoff_for_minor_allele <- 0.4 	## set 40% as approximation for 50%, but please EDIT this in the future, this is a hard cutoff, therefore not the most accurate. Better way might be to model distribution, LOD scores, find the likelihood of a value being within the 50% range......

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

assign_mouse_calls <- function(info_row){
	## if total number of reads does not meet minimum coverage of 20 reads, then assign indeterminate (not enough power to determine)
	if(info_row[['total_cvg']]<20){
		prelim_call <- 'indeterminate'
		mouse_AF <- NA
	} else {
		if(info_row[['Tumor_Seq_Allele2']] %in% c(info_row[['major_allele']], info_row[['minor_allele']])){
			prelim_call <- 'mouse_contamination'
			mouse_AF <- info_row[[info_row[['Tumor_Seq_Allele2']]]]
		} else if (info_row[['Reference_Allele']] %in% c(info_row[['major_allele']], info_row[['minor_allele']])){
			prelim_call <- 'human'
			mouse_AF <- 0
		} else {
			prelim_call <- 'indeterminate'
			mouse_AF <- NA
		}
	}
	
	return(data.frame(prelim_call=prelim_call, mouse_AF=mouse_AF, stringsAsFactors=F))
}

##################################################################
#### 1. INPUT FILES 
#### a. Reading interval_list file denoting positions for pileups
#### b. Reading BAM files for samples
##################################################################

parser <- ArgumentParser(description='Run Mouse QC for a given sample')
parser$add_argument('--sample_id', metavar='sample_id', type="character",
                   help='sample_id', required=T)
parser$add_argument('--bam_path', metavar='bam_path', type="character",
                   help='BAM PATH', required=T)
parser$add_argument('--mouse_reference_snps', type="character",
                   help='Path to mouse_reference_snps', required=T)

args <- parser$parse_args()
sample_id <- args$sample_id
bam_path <- args$bam_path
mouse_reference_snps <- args$mouse_reference_snps

## Interval list file
interval_list_path <- mouse_reference_snps
interval_list <- read.maf(interval_list_path)
interval_list.granges <- makeGRangesFromDataFrame(interval_list, seqnames.field='Chromosome',  start.field='Start_Position', end.field='End_Position', keep.extra.columns=T)


res_pileup <- read_bam_and_run_pileup(bam_path, sample_id, interval_list.granges)

res_pileup_with_allele <- 	
	res_pileup %>% 
	group_by(sample, seqnames, pos) %>%
	mutate(
		pct_count=count/sum(count), 
		pct_count=ifelse(is.na(pct_count), 0, pct_count),
		major_allele=sort(nucleotide[which(pct_count==max(pct_count))])[1],
		minor_allele=select_minor_allele(nucleotide, pct_count, major_allele),
		total_cvg = sum(count)
		) %>% 
	ungroup()

res_pileup_freq <- dcast(res_pileup_with_allele, sample + snp_id + seqnames + pos + Reference_Allele + Tumor_Seq_Allele2 + major_allele + minor_allele + total_cvg ~ nucleotide, value.var='count') %>% unfactorize()
res_pileup_pctCount <- dcast(res_pileup_with_allele, sample + snp_id + seqnames + pos + Reference_Allele + Tumor_Seq_Allele2 + major_allele + minor_allele + total_cvg ~ nucleotide, value.var='pct_count') %>% unfactorize()
called_pileup_pctCount <- adply(res_pileup_pctCount, 1, function(x) assign_mouse_calls(x))

write.table(called_pileup_pctCount, paste0(sample_id, ".mouse_qc.txt"), sep="\t", quote=FALSE)

