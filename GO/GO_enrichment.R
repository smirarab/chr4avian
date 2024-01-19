# Script for GO enrichment analysis of the outlier genes
# Authors Guangji Chen (email: chenguangji@genomics.cn)

rm(list = ls())
#R 4.3.1 (2023-06-16 ucrt)

#BiocManager::install("clusterProfiler")
library("clusterProfiler") #version 4.10.0
#BiocManager::install("AnnotationDbi")
require("AnnotationDbi") #version 1.64.1
#install.packages("dplyr") 
library("dplyr") #version 1.1.2
#install.packages("writexl") 
require("writexl") #verison 1.4.2
require("ggplot2") #version3.4.4

#### using the geneName based on NCBI annotation GCF_000002315.5, galGal6
outlierN.list<-read.table("./outlier.geneName.txt")
wgN.list<-read.table("./whole-genome.geneName.txt")

#### using the gprofiler gmt as database annotations BioMart classes: releases/2023-07-27
term2gene<-read.gmt("./gprofiler_full_ggallus.name.gmt") #annotations BioMart classes: releases/2023-07-27
### whole genome
df.go.redo.wg<-enricher(gene=outlierN.list$V1, universe=wgN.list$V1,
                        pvalueCutoff = 1, qvalueCutoff = 1,pAdjustMethod = "BH",
                        TERM2GENE = term2gene)
df.go.redo.wg.add <- as.data.frame(df.go.redo.wg)  %>%
  left_join(go2term(df.go.redo.wg$ID) %>% rename(ID=go_id,term=Term)) %>%
  left_join(go2ont(df.go.redo.wg$ID) %>% rename(ID=go_id,Ont=Ontology)) %>% na.omit()
df.go.redo.wg.add.p <- df.go.redo.wg.add %>% mutate(Count=as.numeric(stringr::str_split_fixed(GeneRatio,"/",2)[,1])) %>%
  select(c("ID","term","Ont","pvalue","p.adjust","qvalue","GeneRatio","Count","geneID")) %>% 
  mutate(GeneRatio =sapply(GeneRatio, function(x) eval(parse(text=x)))) %>% 
  arrange(GeneRatio) %>% mutate(Ont=factor(Ont,levels = c("MF","BP","CC"))) %>% arrange(desc(Ont)) #%>% filter(Ont!="NA",pvalue<0.05)

write_xlsx(df.go.redo.wg.add,"./WG_gprofile2_all_GO_results.xlsx")


### plot the data with gene & GOs
dataset0<-t(df.go.redo.wg.add %>% mutate(Term = paste0(term,"\n(p.adjust = ",round(p.adjust,4),")")) %>% select(c("Term","geneID")))
colnames(dataset0)=dataset0[1,]
dataset0=dataset0[-1,]
dataset0.list<-as.list(dataset0)

total.list1=c()
dataset1.list=c()
for (i in 1:4) {
  x=names(dataset0.list)[[i]]
  dataset1.list[[x]]=unlist(strsplit(dataset0.list[[x]],split = "/"))
  total.list1<-c(total.list1,dataset1.list[[x]])
}

cnetplot(dataset1.list, showCategory = 4,  cex.params = list(category_node = 0.7, gene_node = 0.7, category_label = 0.7, gene_label = 0.7), 
         foldChange=table(total.list1), circular = F, colorEdge = TRUE) + 
  scale_color_gradient2(name='Associated GO counts', low='darkgreen', high='firebrick')+
  labs(title = "Top 4 GO terms (p.adjusted < 0.05)")
ggsave("./Figure_GO_genes.pdf", w=7,h=6)
