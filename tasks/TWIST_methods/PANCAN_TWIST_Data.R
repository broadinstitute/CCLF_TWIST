gcp_Data=read.table("~/Downloads/sample.tsv",header = T,stringsAsFactors = F)
gcp_Data_sub=gcp_Data[,c("entity.sample_id","crai_or_bai_path","cram_or_bam_path")]
gcp_Data_sub$cram_or_bam_path=gsub("fc-51648132-739d-4b71-9800-03fb9f9990c6","fc-26a55b2a-ed2f-4c75-92cc-c4f9b943d92d",gcp_Data_sub$cram_or_bam_path)
gcp_Data_sub$crai_or_bai_path=gsub("fc-51648132-739d-4b71-9800-03fb9f9990c6","fc-26a55b2a-ed2f-4c75-92cc-c4f9b943d92d",gcp_Data_sub$crai_or_bai_path)
gcp_Data_sub$bsp_sample_id_validation=do.call(rbind, strsplit(gcp_Data_sub$entity.sample_id, '\\_'))[,3]
pancan_external=data.frame(openxlsx::read.xlsx("~/Downloads/Untitled spreadsheet.xlsx",sheet=1),stringsAsFactors = F)
pancan_external_sub=pancan_external[,c("Sample.ID","Stock.Sample","Collection","Original.Material.Type","Material.Type","Collaborator.Participant.ID","Sample.Type","Tumor.Type","Tissue.Site")]

pancan_external2=data.frame(openxlsx::read.xlsx("~/Downloads/PanCan1 External Sample ID.xlsx",sheet=1),stringsAsFactors = F)
pancan_external2_sub=pancan_external2[,c("Exported.DNA.SM.ID","External.ID")]

pancan_external_combined=merge(pancan_external2_sub,pancan_external_sub,by.x=c("Exported.DNA.SM.ID"),by.y=c("Sample.ID"),all.x=T)
pancan_external_combined$sample_id=paste(pancan_external_combined$Collaborator.Participant.ID,pancan_external_combined$Sample.Type,pancan_external_combined$Exported.DNA.SM.ID,sep="-")
colnames(pancan_external_combined)=c("bsp_sample_id_validation","external_id_validation","stock_sample_id_validation","Collection","source_subtype_validation","processed_subtype_validation","participant","sample_type","tumor_subtype","tissue_site","sample_id")
pancan_external_combined$tsca_id="PANCAN1"

pancan_final=merge(pancan_external_combined,gcp_Data_sub,by=c("bsp_sample_id_validation"))
colnames(pancan_final)[13]="reference_id"
pancan_final$participant_id=pancan_final$participant


#### CHROMOSOME PLOT##################

common_twist_default=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_pancan_default.bed",header = F,stringsAsFactors = F)
common_tsca_default=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_tsca_default.bed",header = F,stringsAsFactors = F)
common_twist_10percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_pancan_10percent.bed",header = F,stringsAsFactors = F)
common_tsca_10percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_tsca_10percent.bed",header = F,stringsAsFactors = F)
common_twist_25percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_pancan_25percent.bed",header = F,stringsAsFactors = F)
common_tsca_25percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_tsca_25percent.bed",header = F,stringsAsFactors = F)
common_twist_50percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_pancan_50percent.bed",header = F,stringsAsFactors = F)
common_tsca_50percent=read.table("/Users/aniketshetty/pancancer_TWIST/bedtool_tsca_50percent.bed",header = F,stringsAsFactors = F)

tsca_all_regions=read.table("~/pancancer_TWIST/plots/rcd2_complete.tsv",header = F,stringsAsFactors = F)
twist_all_regions=read.table("~/pancancer_TWIST/plots/PANCAN_complete.tsv",header = F,stringsAsFactors = F)

tsca_all_regions$design_function="UNIQUE_TO_TSCA"
twist_all_regions$design_function="UNIQUE_TO_TWIST"
colnames(tsca_all_regions)=c("chr","pos1","pos2","id","design_function")
colnames(twist_all_regions)=c("chr","pos1","pos2","id","design_function")

tsca_all_regions_default=tsca_all_regions
twist_all_regions_default=twist_all_regions
tsca_all_regions_10percent=tsca_all_regions
twist_all_regions_10percent=twist_all_regions
tsca_all_regions_25percent=tsca_all_regions
twist_all_regions_25percent=twist_all_regions
tsca_all_regions_50percent=tsca_all_regions
twist_all_regions_50percent=twist_all_regions

#tsca_all_regions$design_function[tsca_all_regions$id %in% common_tsca$V4]="COMMON_TO_BOTH"
#twist_all_regions$design_function[twist_all_regions$id %in% common_twist$V4]="COMMON_TO_BOTH"

tsca_all_regions_default$design_function[tsca_all_regions_default$id %in% common_tsca_default$V4]="COMMON_TO_BOTH"
twist_all_regions_default$design_function[twist_all_regions_default$id %in% common_twist_default$V4]="COMMON_TO_BOTH"

tsca_all_regions_10percent$design_function[tsca_all_regions_10percent$id %in% common_tsca_10percent$V4]="COMMON_TO_BOTH"
twist_all_regions_10percent$design_function[twist_all_regions_10percent$id %in% common_twist_10percent$V4]="COMMON_TO_BOTH"

tsca_all_regions_25percent$design_function[tsca_all_regions_25percent$id %in% common_tsca_25percent$V4]="COMMON_TO_BOTH"
twist_all_regions_25percent$design_function[twist_all_regions_25percent$id %in% common_twist_25percent$V4]="COMMON_TO_BOTH"

tsca_all_regions_50percent$design_function[tsca_all_regions_50percent$id %in% common_tsca_50percent$V4]="COMMON_TO_BOTH"
twist_all_regions_50percent$design_function[twist_all_regions_50percent$id %in% common_twist_50percent$V4]="COMMON_TO_BOTH"

both_combined=rbind(twist_all_regions,tsca_all_regions)
both_combined$id=NULL
colnames(both_combined)=c("chromosome","start","end","design_function")
chr_length_guide <- read.delim('Desktop/ucsc_chromosome_length_guide_hg19_20150218.txt', stringsAsFactors=F)


p_total <- ggplot() + 
  geom_bar(data=chr_length_guide, aes(x=Chrom, y=size), position="dodge", stat="identity", fill="grey80") +
  # geom_bar(data=chr_length_guide, aes(x=Chrom, y=size), position="dodge", stat="identity", fill="white", color='black') +
  geom_errorbar(data=both_combined, aes(x=chromosome, y=start, ymin=start, ymax=end, colour = design_function), width=0.8) +
  scale_color_discrete(name = "Category") +
  coord_flip() + 
  ggtitle(paste0('Rapid Cancer Detection Panel Design plot')) +
  scale_x_discrete(limits=rev(c(1:22,"X","Y"))) + 
  theme_bw() + 
  theme(legend.key=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  xlab("Chromosome")


ggsave(filename = paste0('~/pancancer_TWIST/plots/DesignPlot.png'), plot=p_total, width=15, height=8)
############

firecloud_twist=read.table("~/pancancer_TWIST/plots/copy_number/PANCAN_TWIST_list.txt",header = T,stringsAsFactors = F)
firecloud_twist$tsca_id=NULL
firecloud_twist$cnv_calls=NULL
colnames(firecloud_twist)[3]="tumor_seg_twist"
firecloud_tsca=read.table("~/pancancer_TWIST/plots/copy_number/TSCA_list.txt",header = T,stringsAsFactors = F)
firecloud_tsca$tsca_id=NULL
firecloud_tsca$cnv_calls=NULL
colnames(firecloud_tsca)[3]="tumor_seg_tsca"
firecloud_combined=merge(firecloud_tsca,firecloud_twist,by=c("external_id_validation"))
write.table(firecloud_combined,file="~/pancancer_TWIST/plots/copy_number/Seg_files_download.txt",sep = "\t",quote = F,row.names = F,col.names = T)
firecloud_files <- list.files("~/pancancer_TWIST/plots/copy_number/tsca", pattern="*.seg", full.names=TRUE, recursive=FALSE)
f=list()
for(j in 1:length(firecloud_files)){
  t=read.table(firecloud_files[j],header=T,stringsAsFactors=F)
  t$Sample=paste(do.call(rbind, strsplit(tools::file_path_sans_ext(tools::file_path_sans_ext(basename(firecloud_files[j]))), '\\.'))[,1],"TSCA",sep="-")
  #  t$ID=tools::file_path_sans_ext(tools::file_path_sans_ext(basename(firecloud_files[j])))
  #colnames(t)=c("Sample","chrom","loc.start","loc.end","num.mark","seg.mean","ID")
  #colnames(t)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  #t=t[,c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")]
  colnames(t)=c("ID","chromosome","start","end","probes","log2")
  t$log2=log2(t$log2)
  t$gene=""
  write.table(t,file=paste("~/pancancer_TWIST/plots/copy_number/edited_tsca/",unique(t$ID),".tumor.seg",sep=""),sep = "\t",quote = F,row.names = F,col.names = T)
  #f[[j]]=t
}
firecloud_cnv_seg_twist=do.call("rbind",f)
write.table(firecloud_cnv_seg_tsca,file="~/pancancer_TWIST/plots/copy_number/TWIST_seg_File_for_plotting.seg",sep = "\t",quote = F,row.names = F,col.names = T)
final_cnv_seg=rbind(firecloud_cnv_seg_twist,firecloud_cnv_seg_tsca)
colnames(final_cnv_seg)=c("ID","chromosome","start","end","probes","log2")
final_cnv_seg$gene=""
write.table(final_cnv_seg,file="~/pancancer_TWIST/plots/copy_number/BOTH_seg_File_for_plotting.seg",sep = "\t",quote = F,row.names = F,col.names = T)



################################# TSCA ######################
data_for_import=read.delim("~/Downloads/TSCA_mapping_for_delete.tsv",header = T,stringsAsFactors = F)
file_names_deleted=read.table("~/TSCA_deleted_fileNames.txt",header = T,stringsAsFactors = F)
colnames_import=c("individual_id","aggregation_product_name_validation","external_id_validation","bsp_sample_id_validation","stock_sample_id_validation","sample_type","picard_aggregation_type_validation","sample_id","tumor_subtype","squid_sample_id_validation","source_subtype_validation","processed_subtype_validation","primary_disease","media","Collection","tissue_site","clean_bam_file_capture")
setdiff(colnames_import,colnames(data_for_import))
data_for_import_ed=data_for_import
colnames(data_for_import_ed)[1]="sample_id"
colnames(data_for_import_ed)[19]="individual_id"
for(i in 22:41){
  j=i+1
  tsca_id=paste("TSCA",j,sep="")
  data_for_import_sub=data_for_import_ed[which(data_for_import_ed$tsca_id==tsca_id),]
  write.table(data_for_import_sub[,c(colnames_import)],file=file_names_deleted$File_name[file_names_deleted$tsca_id==tsca_id],sep = "\t",col.names = T,row.names = F,quote = F)
}
