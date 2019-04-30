#!/bin/bash
## Usage: Plot a bunch o' diagnostic plots to figure out what the issue is.
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 28th January 2019


# library(reshape)
library(ggplot2)
library(sinaplot)
# library(plyr)
# library(lemon)
setwd("/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug")

## Read in genotype information per sample/variant
vcf <- read.table("BATCHmerge_whole.nosq.txt", header = TRUE, stringsAsFactors = F)

# Get the column numbers which contain the samples
sample.cols <- which(grepl("EGAN", colnames(vcf), fixed=TRUE))



# Iterate over the samples, and sum how many of each genotype exists
comp=data.frame()
 # 
for (i in sample.cols) {
  
  raw=data.frame()
  x=vcf[,i]
  
  DEL=x[which(vcf[,2]=="DEL")]
  DEL.REF=table(DEL)["0/0"]
  DEL.HET=table(DEL)["0/1"]
  DEL.HOM=table(DEL)["1/1"]
  del.raw=cbind(DEL.REF,DEL.HET,DEL.HOM)
  raw=rbind(raw,del.raw[1,])
  
  DUP=x[which(vcf[,2]=="DUP")]
  DUP.REF=table(DUP)["0/0"]
  DUP.HET=table(DUP)["0/1"]
  DUP.HOM=table(DUP)["1/1"]
  dup.raw=cbind(DUP.REF,DUP.HET,DUP.HOM)
  raw=rbind(raw,dup.raw[1,])
  
  INV=x[which(vcf[,2]=="INV")]
  INV.REF=table(INV)["0/0"]
  INV.HET=table(INV)["0/1"]
  INV.HOM=table(INV)["1/1"]
  inv.raw=cbind(INV.REF,INV.HET,INV.HOM)
  raw=rbind(raw,inv.raw[1,])    
  
  raw$SVTYPE=c("DEL", "DUP", "INV")
  raw$sample=names(vcf[i])
  colnames(raw)=c("REF","HET","HOM", "SVTYPE", "SAMPLE")
  
  comp=rbind(comp,raw)
  
}
# Remove the "gt" from the samplenames
comp$SAMPLE <- gsub(".GT", "", comp$SAMPLE)

# add batch number
batch.nums <- read.table("sampleList_batchname.txt", header = F, stringsAsFactors = F)
colnames(batch.nums) <- c("sample", "batch")
comp.with.batch <- merge(comp,batch.nums, by.y = "sample", by.x = "SAMPLE")

# add a column for homs/hets
comp.with.batch$prophet <- comp.with.batch$HET/comp.with.batch$HOM
comp.with.batch$prophom <- comp.with.batch$HOM/comp.with.batch$HET
comp.with.batch$propNonRef <- (2*(comp.with.batch$HOM)+comp.with.batch$HET)/comp.with.batch$REF

# Save data
write.table(comp.with.batch, "frequency_genotypes.txt", quote=F, row.names=F, sep="\t")


svtype.cols <- c("svpipeline" = "#D95F02", "1000g" = "#7570B3", "genomestrip" = "#E7298A", "BND" = "#66A61E")
zyg.cols <- c("Het" = "#D95F02", "Hom" = "#7570B3", "Ref" = "#E7298A")


## Make some violins describing number of each type per sample
pdf(width = 10, height = 3, file = "per_batch_het_per_hom.pdf")
ggplot(comp.with.batch) + #[which(comp.with.batch$SVTYPE=="DUP"),]
  geom_violin(aes(x=batch, y=propNonRef), col = "#D95F02", fill = "#D95F02", alpha = 0.5) + 
  # geom_violin(aes(x=batch, y=HOM), col = "#7570B3", fill = "#7570B3", alpha = 0.5) + 
  # geom_violin(aes(x=batch, y=REF), fill = "#E7298A") + 
  # scale_fill_manual(values = zyg.cols, 
                      # name="Type",
                      # labels=c("Het", "Hom", "Ref")) + 
  ylab("Number of Hets per Hom") +
  xlab("") +
  facet_grid(.~comp.with.batch$SVTYPE) +
  theme_bw() 
graphics.off()




  
  
ggplot(comp.with.batch[which(comp.with.batch$SVTYPE=="DEL"),]) + 
  geom_violin(aes(x=batch, y=HOM, fill = batch)) + 
  
  # facet_wrap(comp.with.batch$SVTYPE)
  theme_bw() 


  

# ## MSQ
#     # Basic density plot per type of MSQ
#     ggplot(vcf.no.bend[which(vcf.no.bend$MSQ<1000 & vcf.no.bend$MSQ>100),]) +
#       geom_density(aes(
#         fill=SVTYPE,
#         x=MSQ
#       ),
#       alpha = 0.2)
    
#     ## three violin plots, all values
#     pMSQ.all <- ggplot(vcf.no.bend[which(vcf.no.bend$MSQ<10000 ),]) +
#                   theme_bw() + 
#                   scale_fill_brewer(palette="Dark2") +
#                   geom_violin(aes(
#                     x=SVTYPE,
#                     y=MSQ,
#                     fill=SVTYPE
#                   ),
#                   alpha=0.7) +
#                   theme(legend.position="none") +
#                   scale_y_continuous(
#                     # breaks=seq(0,1000,100), 
#                     name="Mean sample quality") +
#                   xlab("SVTYPE\nMSQ < 1000\nn=53718") 
    
    
#     ## Three violins, only quality less than 1000
#     pMSQ.1000 <- ggplot(vcf.no.bend[which(vcf.no.bend$MSQ<1000 ),]) +
#                   theme_bw() + 
#                   scale_fill_brewer(palette="Dark2") +
#                   geom_violin(aes(
#                     x=SVTYPE,
#                     y=MSQ,
#                     fill=SVTYPE
#                   ),
#                   alpha=0.7) +
#                   theme(legend.position="none") +
#                   scale_y_continuous(
#                     breaks=seq(0,1000,100), 
#                     name="Mean sample quality") +
#                   xlab("SVTYPE\nMSQ < 1000\nn=53121") 
    
#     ## Three violins, less than 1000, more than 100
#     ggplot(vcf.no.bend[which(vcf.no.bend$MSQ<1000 & vcf.no.bend$MSQ>100 ),]) +
#       theme_bw() + 
#       scale_fill_brewer(palette="Dark2") +
#       geom_violin(aes(
#         x=SVTYPE,
#         y=MSQ,
#         fill=SVTYPE
#       ),
#       alpha=0.7) +
#       theme(legend.position="none") +
#       scale_y_continuous(
#         breaks=seq(0,1000,100), 
#         name="Mean sample quality") +
#       xlab("SVTYPE\nMSQ > 100, MSQ < 1000\nn=24584")
    
    
#     ggsave(filename = "MSQ-all.pdf", plot = pMSQ.all,  width = 7, height = 8, units = 'cm' )
#     ggsave(filename = "MSQ-1000.pdf", plot = pMSQ.1000,  width = 7, height = 8, units = 'cm' )

# ##
# # Just use PLINK - HW
# # /lustre/scratch118/infgen/team133/db22/software/plink1.9/plink --vcf /lustre/scratch119/realdata/mdt2/projects/cnv_15x/svtools/Phase-I/BATCHmerge/WDL_scripts/cromwell-executions/Post_Merge_SV/48e5155e-4aeb-4718-b698-c86a256855ac/call-Sort_Index_VCF/execution/BATCHmerge.vcf.gz --hardy --out HW.txt


# ## Make vcf minimal and remove unecessary columns
#       vcf.reduce <- subset(vcf.no.bend, select=c("ID", "MSQ", "AF", "SVTYPE", "NSAMP", "SVLEN"))

# ## Read in HW table and remove unecessary columns
#       HW.raw <- read.table("BATCHmerge.hw", header = T)
#       HW.reduce <- subset(HW.raw,select=c ("SNP" , "P"  ))

# ## Merge VCF with HW data
#       vcf.with.hw <- merge(vcf.reduce,HW.reduce, by.y = "SNP", by.x = "ID")
#       vcf.with.hw$SVLEN <- as.numeric(vcf.with.hw$SVLEN)

# ## Add a size column, with delimiters 500, 10000, and larger
#       vcf.with.hw$size <-  sapply(vcf.with.hw$SVLEN, function(x){
#         if (abs(x) > 0 & abs (x) < 500 ) {
#           return ("small") 
#         } else if (abs(x) > 499 & abs (x) < 10000 ) {
#           return ("medium") 
#         } else if (abs(x) > 9999 ) {
#           return ("large") 
#         }
        
#       })
      
      
# ## Plot HW against AF, in bins      
      
#       hw <- ggplot(data=vcf.with.hw[which(vcf.with.hw$AF<0.5),]) +
#                 geom_boxplot(
#                   aes(
#                     x=factor(as.character(round_any(AF, 0.05)), 
#                              level = seq(0,1,0.05)
#                     ),
#                     y=log(P)
#                   ),
#                   outlier.shape = NA,
#                   fill="468",
#                   alpha=0.3
#                 ) + 
#                 scale_x_discrete(breaks=seq(0,0.5,0.01), name="minor allele frequency") +
#           ylab("log10 p-value") +
#               # scale_y_continuous(trans='log10') +
#                 # facet_grid( size~ .) +
#                 theme_bw() +
#                 theme(axis.text.x = element_text(angle = 40, hjust = 1), 
#                       legend.position = "none") 
      
      
#       ggsave(filename = "~/Documents/reports/first-year-report/figIdeas/HW-by-AF.pdf", plot = hw,  width = 14, height = 7, units = 'cm' )
  