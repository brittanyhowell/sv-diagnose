
# Usage qq.plink.R gwas.txt "A plot title"
# Note: Required that columns name with P values is labeled as 'P'.
setwd("~/Documents/data/Phase-I/debug/caseControl/")


library(ggplot2)

# args <- commandArgs(trailingOnly = TRUE)

# plink.f = "BATCHmerge.1.out.assoc"
title = "Batches 1 & 2"
# qq.f = paste(plink.f, ".qq.png", sep="")

# message("Input file = ", plink.f)
# message("Title = ", title)
# message("Output file = ", qq.f)

# message("Reading data ...")

d <- read.table("BATCHmerge.1.2.A.assoc", header=TRUE)
d <- read.table("BATCHmerge.1.out.assoc", header=TRUE)

# <- tryCatch({ readRDS(plink.f) },
#               warning = function(war) { read.table(plink.f, header=TRUE) },
#               error = function(err) { read.table(plink.f, header=TRUE) },
#               finally = {
#                 message("")
#               }
# )

do <- d

obs <- sort(d$P)
obs <- obs[!is.na(obs)]
obs <- obs[is.finite(-log10(obs))]

lobs <- -(log10(obs))

exp <- c(1:length(obs)) 
lexp <- -(log10(exp / (length(exp)+1)))

chisq <- qchisq(1-obs,1)

lambda <- round(median(chisq)/qchisq(0.5,1),2)




labelA<-as.expression(bquote(lambda == .(lambda)))
labelB<- as.expression(bquote(N == .(length(lobs))))	

# message("Drawing Plot")
png(file="BATCH1.png", res=300, height = 1500, width = 2000)
ggplot()+
  geom_point(aes(lexp, lobs), alpha = 0.1) +
  geom_abline(intercept=0,slope=1) +
  xlab(expression("Expected (-log"[10]*"(P))")) +
  ylab(expression("Observed (-log"[10]*"(P))")) +
  annotate(geom="text", x=4, y=max(lobs)/10*2, label=labelA,
           color="black") +
  annotate(geom="text", x=4, y=max(lobs)/10*1, label=labelB,
           color="black") +
  theme_bw()
dev.off()

# Plot p 
# for (chr in 1:5) {
png(file="chr1-5_batch1.png", res=300, height = 2000, width = 2000)

ggplot(d[which(d$CHR==1 | d$CHR==2 | d$CHR==3 | d$CHR==4 ),])+
        geom_point(aes(BP, -log10(P)), alpha = 0.1) +
        facet_grid(CHR~.) +
        theme_bw() +
        geom_hline(aes(yintercept=0.05))
  # ggprint(p)
  dev.off()
# }
