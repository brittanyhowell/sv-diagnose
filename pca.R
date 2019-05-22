library(ggfortify)

convertSVTYPE=function(x) {
  if (x == "DUP"){
    val <- 1
  } else if (x == "DEL"){
    val <- 2
  } else if (x == "INV"){
    val <- 3
  }
  return(val)
}

setwd("~/Documents/data/Phase-I/debug/pca/")
rand.vars <- read.table("BATCHall_rand.txt", stringsAsFactors = F)
colnames(rand.vars) <- c("chr", "pos", "ID", "QUAL", "SVTYPE", "AF", "NSAMP", "MSQ", "batch")




rand.vars=subset(rand.vars,select = -c(batch, chr, pos))

# vcf.reduce <- vcf[,c(7,11,12,13,17,18,19)]

rand.vars$SVTYPE=sapply(rand.vars$SVTYPE, convertSVTYPE)

# PCA with function PCA
library(FactoMineR)

# apply PCA
pca3 = PCA(rand.vars)
           


rand.vars=as.data.frame(apply(rand.vars,2,function(x){as.numeric(x)}))
# vcf.reduce$SVLEN=abs(vcf.reduce$SVLEN)




autoplot(prcomp(rand.vars),
         colour = 'SVTYPE',
         loadings = TRUE,
         loadings.label = TRUE, loadings.label.size = 3)


## USING PLINKEROO


setwd("/Users/bh10/Documents/data/Phase-I/debug/plinkPca/bcftoolsExclude///")

# options(scipen=100, digits=3)

#Read in the eigenvectors
eigenvec <- data.frame(read.table("SU_pca.eigenvec", header=FALSE), stringsAsFactors = F)
eigenvec$V2 <- as.character(eigenvec$V2)
eigenval <- data.frame(read.table("SU_pca.eigenval", header=FALSE, skip=0, sep=" "))
rownames(eigenvec) <- eigenvec[,2]



key <- read.table("../../sampleList_batchname.txt", header = F, stringsAsFactors = F)

eigen.merge <- (merge(key, eigenvec, all=TRUE, by.y = "V2", by.x='V1'))
eigen.merge <- eigen.merge[,2:ncol(eigen.merge)]
colnames(eigen.merge)[1]="BATCH"

sum_eigs<-sum(eigenval$V1)
sum_eigs<-lapply(eigenval$V1,function(x){
  rt<-(x/sum_eigs)*100
  rt<-round(rt)
  return(rt)
})


pdf("SU_5.pdf", width = 6, height = 4)
ggplot(eigen.merge)+
  geom_point(aes(eigen.merge[,3], eigen.merge[,4],
                 col = BATCH) , alpha=0.9) +
  xlab(paste0("PC1 (",sum_eigs[[1]],"% variance)")) +
  ylab(paste0("PC2 (",sum_eigs[[2]],"% variance)")) +
  theme_bw() +
  scale_colour_brewer(palette="Dark2")  + 
  theme(panel.grid=element_blank(), legend.title=element_blank()) + 
  theme(legend.position=c(.1, .2))
graphics.off()
  
