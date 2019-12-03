############################################
# Query samples from batch in fingerprinting db
# mcadosch 10/17
############################################

rm(list=ls())
suppressMessages(library(Rsamtools))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressPackageStartupMessages(library("argparse"))

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

write.delim <- function(data, file='', coln=T, rown=F, append=F, ...){
    write.table(data, file, col.names=coln, sep="\t", row.names=rown, quote=F, append=append)
}

query_fp_database <- function(mySample, myDatabase.fp, cutoff_similarity=0.70){
	## 70% cutoff might be better than 80%
	
	cat('Querying', mySample, '\n')
	
	myQuery.fp <- subset(myDatabase.fp, sample==mySample, select=c('sample', 'snp_id', 'genotype'))

	#### Selecting most likely observed genotype for each position per sample based on LOD score - already done so with major/minor allele 

	#### Do comparison between query & fp_database:
	database_plus_query <- 
		merge(myQuery.fp, myDatabase.fp, all.x=T, all.y=T, suffixes=c('.query', '.database'), by='snp_id') %>%
		mutate(
			does_query_match_database=(genotype.query==genotype.database)
			)
	
	scores_per_sample_in_database <- 
		database_plus_query %>% 
		group_by(sample.query, sample.database) %>% 
		summarise(
			n_powered_snps= sum(!is.na(does_query_match_database)),
			pct_matching=sum(does_query_match_database, na.rm=T)/n_powered_snps
			) %>% 
		ungroup() %>%
		arrange(desc(pct_matching), desc(n_powered_snps))
		
		
	scores_per_sample_in_database.85pt <- 
		subset(scores_per_sample_in_database, pct_matching >=cutoff_similarity) %>% 
		arrange(desc(n_powered_snps), desc(pct_matching)) 

	if(dim(scores_per_sample_in_database.85pt)[1]==0){
		cat('...fingerprinting query completed. No likely matches (>', cutoff_similarity, ') found.\n')
	} else {
		cat('...fingerprinting query completed. Mostly likely matches (>', cutoff_similarity, ') are:\n')
		print(data.frame(scores_per_sample_in_database.85pt))
	}

	return(scores_per_sample_in_database)
}

getIndivID = function(sample_id){
	sampleName1 = gsub("CCLF_", "", sample_id);
	sampleName2 = gsub("CCLFTSCA-", "", sampleName1);
	splitData = strsplit(sampleName2, "-")
	return(splitData[[1]][1])
}


##################################################################
#### MAIN
##################################################################
parser <- ArgumentParser(description='Query samples from batch in fingerprinting db')
parser$add_argument('--sample_ids', type="character", nargs="+",
                   help='sample ids', required=T)
parser$add_argument('--external_ids', type="character", nargs="+",
                   help='sample external id validations', required=T)
parser$add_argument('--tsca_id', type="character", 
					help='tsca_id', required=T)
parser$add_argument('--fingerprinting_db', type="character", 
					help='fingerprinting_db file', required=T)
parser$add_argument('--fluidgm_snp_path', type="character",
                   help='Path to fluidgm_snp database', required=T)

args <- parser$parse_args()
sample_ids <- args$sample_ids
external_ids <- args$external_ids
tsca_id <- args$tsca_id
fingerprinting_db <- args$fingerprinting_db
fluidgm_snp_path <- args$fluidgm_snp_path


### Read resources
fluidgm_snp.raw <- read.maf(fluidgm_snp_path)
fluidgm_snp <- fluidgm_snp.raw[, c('CHROM', 'POSITION', 'SNP', 'SNP_ALLELES')] %>% unique()
fluidgm_snp.granges <- makeGRangesFromDataFrame(fluidgm_snp, seqnames.field='CHROM', start.field='POSITION', end.field='POSITION', keep.extra.columns=T)

### Create FNG DB
myDatabase = read.delim(fingerprinting_db, stringsAsFactors=F)
myDatabase.fp = myDatabase %>% subset(select=c('sample', 'snp_id', 'genotype'))

print('myDatabase.fp')
print(head(myDatabase.fp))

### Create sample_info DF
print(sample_ids)
print(external_ids)
sample_info <- data.frame("sample_id"=sample_ids, "external_id"=external_ids)
print(sample_info)
res_fp <- ddply(sample_info, .(sample_id), function(x) query_fp_database(x[['external_id']], myDatabase.fp))

print('res_fp')
print(head(res_fp))

write.delim(res_fp, paste0(tsca_id,  '.res_fp.txt'))

### Fix some IDs
res_fp$sampleName = sapply(res_fp$sample_id, getIndivID)

res_fp.detailed = res_fp %>% 
    mutate(
        hasSameIndivID = mapply(grepl, sampleName,sample.database),
        sample.database.batch = NA,
        sample.query.batch  = NA) %>% 
        subset((sample.query != sample.database)& !(sample.database %in% c('3T3J2','PEDS002TW','MP_pos'))) %>% 
        subset(hasSameIndivID == TRUE | pct_matching > 0.7
    )

print('res_fp.detailed')
print(head(res_fp.detailed))

if (nrow(res_fp.detailed) >0){
	for (ind in c(1:nrow(res_fp.detailed ))) {
		sample.database = res_fp.detailed$sample.database[ind]
		matchInd = which(myDatabase$sample == sample.database)[1]
		res_fp.detailed$sample.database.batch[ind] = myDatabase$batch[matchInd]
		
		sample.query = res_fp.detailed$sample.query[ind]
		matchIndquery = which(myDatabase$sample == sample.query)[1]
		res_fp.detailed$sample.query.batch[ind] = myDatabase$batch[matchIndquery]
	} 
}

### Check if samples from the same individual match (should)
res_fp_mstLikelyMatches_summary.sameIndivs = res_fp.detailed %>%  subset(hasSameIndivID == TRUE)

print('res_fp_mstLikelyMatches_summary.sameIndivs')
print(head(res_fp_mstLikelyMatches_summary.sameIndivs))

res_fp_mstLikelyMatches_summary.notMatchingIndivs = res_fp_mstLikelyMatches_summary.sameIndivs %>% 
	group_by(sample_id) %>%
	summarise(num_indivs = n(),
		num_matching = length(which(pct_matching > 0.70)),
		frac_matching = num_matching/num_indivs,
		failed_FNG_IndivsMatch = ifelse(frac_matching < 0.75, 1,0)
	) 

failedSamples = res_fp_mstLikelyMatches_summary.notMatchingIndivs %>% subset(failed_FNG_IndivsMatch == 1, select = sample_id)
failedSamples.NotMatchingIndivs = res_fp_mstLikelyMatches_summary.sameIndivs %>% subset(sample_id %in% failedSamples$sample_id, select = c(-sampleName,-hasSameIndivID) )

print('failedSamples.NotMatchingIndivs')
print(head(failedSamples.NotMatchingIndivs))

write.delim(failedSamples.NotMatchingIndivs, paste0(tsca_id,  '.failedSamples.NotMatchingIndivs.txt'))

### Check between samples that are not from the same individual
res_fp_mstLikelyMatches_summary.final = res_fp.detailed %>%  subset(hasSameIndivID == FALSE & pct_matching > 0.7) 

print('res_fp_mstLikelyMatches_summary.final')
print(head(res_fp_mstLikelyMatches_summary.final))

write.delim(res_fp_mstLikelyMatches_summary.final, paste0(tsca_id,  '.res_fp_mstLikelyMatches_summary.final.txt'))

