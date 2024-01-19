require(tidyverse)
require(ggplot2)
require(reshape2)


g = read.csv('WG.tsv',sep="\t",h=T)
head(g$pvalue)
head(p.adjust(g$pvalue,method = "BH"))
head(g$p.adjust)
g$bgCount = as.numeric(sub("/.*","",g$BgRatio))

bgCountD = unique(as.numeric(sub(".*/","",g$BgRatio)))
gCountD = unique(as.numeric(sub(".*/","",g$GeneRatio)))
g$Enrichment = with(g,(Count/gCountD)/(bgCount/bgCountD))
head(g)

head(g)
ggplot(aes(y=Enrichment,x=reorder(term,p.adjust),
           color=cut(p.adjust,c(0,0.01,0.05,0.1,0.5,1)),
           size=Count),data=g[g$Enrichment>2 &g$Count>=4,])+
  geom_point()+  theme_classic()+
  scale_color_viridis_d(name="BH ajsuted p")+
  xlab("")+scale_y_continuous(name = "Term enrichment in outlier region", trans = "log10")+
  #geom_hline(yintercept = 1,color="red")+
  theme(axis.text.x = element_text(size=6,angle=90,hjust = 1),
        axis.text.y = element_text(size=9),
        legend.position = c(.7,.8),
        legend.direction = "horizontal")
ggsave("genes-GO.pdf",width=10,height = 6)

genes = merge(melt(strsplit(g2$geneID,"/")), 
      data.frame(desc=g2$term,color=g2$p.adjust,L1=1:nrow(g2)),
      by = "L1")

head(genes)    
ggplot(data=genes,aes(y=reorder(desc,value),x=reorder(value,desc),
                      fill=color))+geom_tile(fill="black")+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(size=4,angle=90,hjust = 1),
        axis.text.y = element_text(size=7))


