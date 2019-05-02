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
vcf <- read.table("BATCHmerge.AF.CN.GT.nosq.txt", header = TRUE, stringsAsFactors = F)
# vcf <- vcf[which((abs(vcf$SVLEN)>100 & abs(vcf$SVLEN)<1000000)),]
# vcf <- vcf[which(vcf$SU>5),]
# vcf <- vcf[which(vcf$AF>0.1),]

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
# comp.with.batch$propNonRef <- (2*(comp.with.batch$HOM)+comp.with.batch$HET)/comp.with.batch$REF

# Save data
write.table(comp.with.batch, "frequency_genotypes_PoM_all.txt", quote=F, row.names=F, sep="\t")

# Read data
setwd("Documents/data/Phase-I/debug/")
comp.with.batch <- read.table("frequency_genotypes_PoM_100_1mb.txt", stringsAsFactors = F, header = T)

# svtype.cols <- c("svpipeline" = "#D95F02", "1000g" = "#7570B3", "genomestrip" = "#E7298A", "BND" = "#66A61E")
# zyg.cols <- c("Het" = "#D95F02", "Hom" = "#7570B3", "Ref" = "#E7298A")


## Make some violins describing number of each type per sample
pdf(width = 10, height = 5, file = "PoM_100_1mb.pdf")
ggplot(comp.with.batch) + #[which(comp.with.batch$SVTYPE=="DUP"),]
  geom_violin(aes(x=batch, y=HET), col = "#D95F02", fill = "#D95F02", alpha = 0.5) + 
  geom_violin(aes(x=batch, y=HOM), col = "#7570B3", fill = "#7570B3", alpha = 0.5) +
  # geom_violin(aes(x=batch, y=REF), fill = "#E7298A") + 
  # scale_fill_manual(values = zyg.cols, 
                      # name="Type",
                      # labels=c("Het", "Hom", "Ref")) + 
  ylab("Number of Hets per Hom") +
  xlab("") +
  facet_grid(.~comp.with.batch$SVTYPE) +
  theme_bw() 
graphics.off()

###### Read in PreMerge VCFs

setwd("/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/")
for (num in 1:5){
  print(num)
d <- dir(paste("/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/BATCH",num,sep=""), full.names=TRUE)
fn <- dir(paste("/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/BATCH",num,sep=""))
comp=data.frame()

for (i in 1:length(d)) {
    
  sample.vars <- read.table(d[i], stringsAsFactors=F)
  colnames(sample.vars) <- c("chr", "pos", "SVTYPE", "QUAL", "gt", "len")
  
  # sample.vars <- sample.vars[which(sample.vars$QUAL>100),]
  # sample.vars <- sample.vars[which(abs(sample.vars$len)>100),]
  sample.vars <- sample.vars[which(abs(sample.vars$len)>100 & abs(sample.vars$len)<1000000),]

  x <- sample.vars$gt
  length(x)

    raw=data.frame()
    
    DEL=x[which(sample.vars[,3]=="DEL")]
    DEL.qual=mean(sample.vars[which(sample.vars[,3]=="DEL"),]$QUAL)
    DEL.REF=table(DEL)["0/0"]
    DEL.HET=table(DEL)["0/1"]
    DEL.HOM=table(DEL)["1/1"]
    del.raw=cbind(DEL.REF,DEL.HET,DEL.HOM,DEL.qual)
    raw=rbind(raw,del.raw[1,])
    
    DUP=x[which(sample.vars[,3]=="DUP")]
    DUP.qual=mean(sample.vars[which(sample.vars[,3]=="DUP"),]$QUAL)
    DUP.REF=table(DUP)["0/0"]
    DUP.HET=table(DUP)["0/1"]
    DUP.HOM=table(DUP)["1/1"]
    dup.raw=cbind(DUP.REF,DUP.HET,DUP.HOM,DUP.qual)
    raw=rbind(raw,dup.raw[1,])
    
    INV=x[which(sample.vars[,3]=="INV")]
    INV.qual=mean(sample.vars[which(sample.vars[,3]=="INV"),]$QUAL)
    INV.REF=table(INV)["0/0"]
    INV.HET=table(INV)["0/1"]
    INV.HOM=table(INV)["1/1"]
    inv.raw=cbind(INV.REF,INV.HET,INV.HOM,INV.qual)
    raw=rbind(raw,inv.raw[1,])    
    

    raw$SVTYPE=c("DEL", "DUP", "INV")
    raw$sample=gsub("_vars.txt","",fn[i])
    colnames(raw)=c("REF","HET","HOM", "meanQUAL" ,"SVTYPE", "SAMPLE")
    
    comp=rbind(comp,raw)
}

batch.nums <- read.table("/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/sampleList_batchname.txt", header = F, stringsAsFactors = F)
colnames(batch.nums) <- c("sample", "batch")
comp.with.batch <- merge(comp,batch.nums, by.y = "sample", by.x = "SAMPLE")

write.table(comp.with.batch, paste("frequency_genotypes_batch",num, "_len_100_1mb.txt", sep = ""), quote=F, row.names=F, sep="\t")
}




# read into local


listSuff=c("_len_100.txt", "_len_100_1mb.txt", "_qual_100.txt","_all.txt")

  suffix=listSuff[4]
  suffix
  pm.batchAll=data.frame()
  
for (num in 1:5) {
pm.batch <- read.table(paste("frequency_genotypes_batch",num,suffix, sep=""), stringsAsFactors = F, header = T)
pm.batchAll <- rbind(pm.batchAll, pm.batch)
}




## Make some violins describing number of each type per sample
pdf(width = 10, height = 5, file = paste("preM_per_batch",gsub(".txt","",suffix),".pdf", sep = ""))
ggplot() + #[which(comp.with.batch$SVTYPE=="DUP"),]
  geom_violin(data = pm.batchAll, aes(x=batch, y=HOM), col = "#D95F02", fill = "#D95F02", alpha = 0.5) +
  geom_violin(data = pm.batchAll, aes(x=batch, y=HET), col = "#7570B3", fill = "#7570B3", alpha = 0.5) +
  geom_violin(data = pm.batchAll, aes(x=batch, y=meanQUAL), col = "red", fill = "red", alpha = 0.1) +
  # geom_violin(data = pm.batchAll[which(pm.batchAll$batch=="BATCH1"),], aes(x=batch, y=HET, fill = SVTYPE), alpha = 0.5) +
  # geom_violin(data = comp.with.batch[which(comp.with.batch$batch=="BATCH1"),], aes(x=batch, y=HET, fill = SVTYPE), alpha = 0.5) + # "#E7298A"
  # scale_fill_manual(values = zyg.cols, 
  # name="Type",
  # labels=c("Het", "Hom", "Ref")) + 
    ylab("Number of hets/homs per sample") +
 xlab("") +
  facet_grid(.~pm.batchAll$SVTYPE) +
  theme_bw() 
graphics.off()




### BATCH4


pdf(width = 10, height = 5, file = "pre_merge_batch4.pdf")
ggplot() + #[which(comp.with.batch$SVTYPE=="DUP"),]
  geom_jitter(data = pm.batchAll[which(pm.batchAll$batch=="BATCH4"),], aes(x=batch, y=HET), col = "#D95F02", fill = "#D95F02", alpha = 0.1) +
  geom_jitter(data = pm.batchAll[which(pm.batchAll$batch=="BATCH4"),], aes(x=batch, y=HOM), col = "#7570B3", fill = "#7570B3", alpha = 0.1) +
  # geom_violin(data = pm.batchAll[which(pm.batchAll$batch=="BATCH1"),], aes(x=batch, y=HET, fill = SVTYPE), alpha = 0.5) +
  # geom_violin(data = comp.with.batch[which(comp.with.batch$batch=="BATCH1"),], aes(x=batch, y=HET, fill = SVTYPE), alpha = 0.5) + # "#E7298A"
  # scale_fill_manual(values = zyg.cols, 
  # name="Type",
  # labels=c("Het", "Hom", "Ref")) + 
  ylab("Number of hets/homs per sample") +
  xlab("") +
  facet_grid(.~pm.batchAll[which(pm.batchAll$batch=="BATCH4"),]$SVTYPE) +
  theme_bw() 
graphics.off()

high.samples <- pm.batchAll$SAMPLE[which(pm.batchAll$HET>2100 & pm.batchAll$batch=="BATCH4" &  pm.batchAll$SVTYPE=="DEL")]





