setwd("~/Documents/data/Phase-I/PostMerge/LD//")


data <- read.table("batch5_chr2.ld", header = T)
head(data)

data.no1 <- data[which(data[,7]<1),]
data.rs <- data.no1 %>% 
  filter(str_detect(SNP_B, "rs"))


data.0.8 <- data.rs[which(data.rs[,7]>.49),]

length(unique(data.0.8[,3]))*4/6837

data <- data.0.8

# png(file="BATCH5-LD.png", res=300, height = 1000, width = 2000)
ggplot(data[which(data$R2<1 &  data$R2>.4 &  abs(data$BP_B-data$BP_A)>1),])+
  geom_segment(aes(xend=BP_A, x=BP_A, y=0,yend=abs(BP_B-BP_A), col=R2), alpha=0.11) + 
  scale_y_continuous(limits = c(0,30000000)) +#
  theme_bw() +
  xlab("location of SV on Chr2 (bp)") +
  ylab("Distance between SV and SNP (bp)") 
# dev.off()
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
