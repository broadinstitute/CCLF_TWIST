####################################
## Chromosome plot of target / bait distribution for TWIST panel ##
####################################
library(ggplot2)

# read in data
twist_all_regions <- read.table("/Users/gmiller/Documents/Work/GitHub/CCLF_TWIST/tasks/TWIST_methods/BroadPanCancer2019.Homo_sapiens_assembly19.targets.txt",header = T,stringsAsFactors = F)
# test_twist <- read.table("/Users/gmiller/Downloads/BroadPanCancer2019.Homo_sapiens_assembly19.targets.interval_list.txt", header = F, stringsAsFactors = F)
# colnames(test_twist)=c("chr","start","end","direction", "chr:start_end")

# colnames(twist_all_regions)=c("chr","pos1","pos2","id","design_function")
# colnames(twist_all_regions)=c("chromosome","start","end","id","design_function")

colnames(twist_all_regions)=c("chr","start","end")
twist_all_regions$id=NULL

chr_length_guide <- read.delim('/Users/gmiller/Documents/Work/GitHub/CCLF_TWIST/tasks/TWIST_methods/ucsc_hg19.chrom.sizes.txt', header = F, stringsAsFactors=F)
colnames(chr_length_guide) <- c("Chrom","size")
chr_length_guide$Chrom <- gsub("chr", "", chr_length_guide$Chrom)

# create plot
p_total <- 
  ggplot() + 
  geom_bar(data=chr_length_guide, aes(x=Chrom, y=size), position="dodge", stat="identity", fill="grey80") +
  # geom_bar(data=chr_length_guide, aes(x=chr, y=size), position="dodge", stat="identity", fill="white", color='black') +
  geom_errorbar(data=twist_all_regions, aes(x=chr, ymin=start, ymax=end, colour = 'salmon'), width=0.6) +
  # geom_errorbar(data=twist_all_regions, aes(x=chr, y=start, ymin=start, ymax=end, colour = design_function), width=0.8) +
  scale_color_discrete(name = "Category") +
  coord_flip() + 
  ggtitle(paste0('TWIST Panel Design')) +
  scale_x_discrete(limits=rev(c(1:22,"X","Y"))) + 
  theme_bw() + 
  theme(legend.key=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "none") +
  xlab("Chromosome") + 
  ylab("Position on chromosome")

# save plot
ggsave(filename = paste0('/Users/gmiller/Documents/Work/GitHub/CCLF_TWIST/tasks/TWIST_methods/DesignPlot.png'), 
       plot=p_total, 
       width=8, height=6)
       # width=15, height=8)

