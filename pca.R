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

df <- iris[c(1, 2, 3, 4)]
vcf.reduce <- vcf[,c(7,11,12,13,17,18,19)]

vcf.reduce$SVTYPE=sapply(vcf.reduce$SVTYPE, convertSVTYPE)

vcf.reduce=as.data.frame(apply(vcf.reduce,2,function(x){as.numeric(x)}))
# vcf.reduce$SVLEN=abs(vcf.reduce$SVLEN)




autoplot(prcomp(vcf.reduce),
         colour = 'SVTYPE',
         loadings = TRUE,
         loadings.label = TRUE, loadings.label.size = 3)
