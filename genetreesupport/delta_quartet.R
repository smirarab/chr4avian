require(ggplot2); require(scales); require(reshape2);library(tidyquant);library(zoo)


rc = read.csv('removed-count.tsv',sep="\t")
rc[is.na(rc)] <- 0
rc$Tree = factor(sub("minus-","",gsub("_([A-Z])"," \\1",rc$Tree)))

q = read.csv('all.stat.xz',sep=" ",header=F)
names(q) = c("Taxa","Gene","main","alt","diff")
q$diff = - q$diff
q$removed = factor(sub("minus-","",gsub("_([A-Z])"," \\1",q$Taxa)))
head(q)
ro = levels(reorder(rc$Tree,100*rc$Count+10*rc$Dove+rc$Cuckoos))

w = read.csv('63K_trees.names_header.txt.xz',sep=" ")
outliers = which((w$Chromosome == "chr4") & ( (w$ws >=25030000 & w$ws <=32670000) |  (w$ws >=33510000  & w$ws <=34470000) | (w$ws >=44130000 & w$ws <=56810000)))

ggplot(aes(x=factor(removed,levels=ro,ordered=T),y=diff),
       data=rbind(data.frame(q,t="All loci"),
                  data.frame(q[!q$Gene %in% outliers,],
                             t="Outlier loci removed")))+
  facet_wrap(~t,nrow=1)+
  stat_summary(aes(color="#2040DE"),fun.data = mean_se,alpha=0.9)+
  stat_summary(aes(color="#20C070"),fun.y = median,alpha=0.9)+
  theme_classic()+
  scale_color_identity(name="",labels=c("mean","median"),guide = "legend")+
  scale_y_continuous(name = "Delta quartet score\n Mirandornithes - Columbea")+
  scale_x_discrete(name="")+
  coord_cartesian(ylim=c(-10000,35000))+
  geom_hline(yintercept = 0,color="red",linetype=1)+
  theme(legend.position = c(.4,.84),
        panel.grid.major.x = element_line(colour = "grey80",linewidth = 0.3),
        axis.text.x = element_blank())
ggsave("taxonsampling.pdf",width=8.8,height = 3.5)

ggplot(aes(x=factor(Tree,levels=ro),y=variable,
           fill=value,label=value),
       data=melt(rc,id.vars = c("Tree","Count")))+
  geom_tile()+
  ylab("")+
  scale_fill_steps2(breaks=1:8,right=F,guide="none")+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_x_discrete(name="")
#geom_text(aes(label=value))+
ggsave("removed-legend-hor.pdf",height=1,width = 4.3)
