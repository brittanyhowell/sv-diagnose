setwd("~/Documents/data/Phase-I/debug/plink-new/")


d <- read.table("plink.ld", header = T)
head(d)


ggplot(d[which(d$R2<1 & abs(d$BP_B-d$BP_A)<10000),])+
  geom_segment(aes(xend=BP_A, x=BP_B, y=0,yend=abs(BP_B-BP_A), col=R2)) +
  theme_bw() +
  xlab("location of SV on Chr1 (bp)") +
  ylab("Distance between SV and SNP (bp)") 

#+
  geom_point(aes(x=BP_A, y=10000,alpha=0.5))
