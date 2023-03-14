require(ggplot2); require(scales); require(reshape2);library(tidyquant);library(zoo)


org = read.csv('outliers-rooted.txt',header=F)

mt = read.csv('mono-time.txt.xz',h=F,sep=" ",col.names = 1:3)
nrow(mt); head(mt)

mtm = read.csv('mono-mutime2.txt.xz',h=F,sep=" ",col.names = 1:3)
nrow(mtm); head(mtm)
mt = merge(mtm,mt,by="X1",all.x = T)
mt[c(444,42300:42310),]
names(mt) = c("Gene","mu", "size", "time" ,"size2")
mt$outlier = mt$Gene %in% org$V1
mt[c(444,42300:42310),]
mt[is.na(mt$size) & !is.na(mt$mu),]
head(mt)


#summary(mt[!is.na(mt$size) & mt$size>8 & mt$time/mt$mu<1800 & !mt$outlier,"time"])
#summary(mt[!is.na(mt$size) & mt$size>8 & mt$time/mt$mu<1800 & mt$outlier,"time"])
summary(mt[!is.na(mt$size) & mt$size>8 & !mt$outlier,"mu"])
summary(mt[!is.na(mt$size) & mt$size>8 & mt$outlier,"mu"])

ggplot(aes(y=mu,x=outlier),data=mt[!is.na(mt$size) & !is.na(mt$mu) & mt$size>8,])+
    theme_classic()+
    scale_y_continuous(name="Columbea branch length (sub unit)")+
    geom_violin(draw_quantiles = c(1/4,1/2,3/4))+
  stat_summary()+
  stat_summary(aes(label=round(..y..,4)),geom="text",position = position_nudge(y = 0.0011))+
  coord_cartesian(ylim=c(0,.033))
ggsave("Columbea-mrca-subunit.pdf",width=3.4,height =4.2)

#############################  
 
rtt = read.csv('OtidimorphaeColumbea-root-to-tip.txt.xz',sep="\t",header=F)
annt = read.csv('annotations.txt',sep="\t")
rtt = merge(rtt,annt,by.x="V1",by.y = "taxa")
head(rtt)
ggplot(aes(x=reorder(V1,V2),y=V2,color=order),data=rtt)+
  geom_violin(draw_quantiles = c(1/4,1/2,3/4))+
  #geom_boxplot()+
  stat_summary()+
  #stat_summary()
  theme_classic()+
  scale_y_continuous(name="tip to MRCA distance (mutation units)")+
  scale_x_discrete(name="")+
  scale_color_brewer(palette = "Dark2", name="order")+
  coord_cartesian(ylim=c(0,.33))+
  theme(axis.text.x = element_text(angle=35,vjust = 1,hjust=1),legend.position = "bottom")
ggsave("MRCA-to-tip.pdf",width=10,height =6)


