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
