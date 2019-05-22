setwd("~/Documents/data/Phase-I/PostMerge/LD//")


data <- read.table("batch5_chr2.ld", header = T)
head(data)

png(file="BATCH5-LD.png", res=300, height = 1000, width = 2000)
ggplot(data[which(data$R2<1 &  data$R2>.4 &  abs(data$BP_B-data$BP_A)>1),])+
  geom_segment(aes(xend=BP_A, x=BP_A, y=0,yend=abs(BP_B-BP_A), col=R2), alpha=0.11) + #
  theme_bw() +
  xlab("location of SV on Chr2 (bp)") +
  ylab("Distance between SV and SNP (bp)") 
dev.off()
#+

length(unique(data[,3]))

  # geom_point(aes(x=BP_A, y=10000,alpha=0.5))
  
  png(file="BATCH1-LD-summary.png", res=300, height = 1000, width = 2000)
  ggplot(data[which(data$R2<1),])+
   geom_point(aes(x=abs(BP_B-BP_A),y=R2),alpha=0.1) +
    theme_bw() +
    scale_x_continuous(limits = c(0,1000000) ) + 
    xlab("Distance between SV and SNP (bp)") +
    ylab("R2 value") 
  dev.off()
