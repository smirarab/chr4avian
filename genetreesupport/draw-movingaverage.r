require(ggplot2); require(scales); require(reshape2);library(tidyquant);library(zoo)
library("ggpubr"); require(imager)


renameTaxa = list ("Columbea" = "Columbea" ,
      "Columbiformes+Otidimorphae"  = "ColumbiformesOtidimorphae"     ,
      "Columbimorphae" ="Columbimorphae",                            
      "Columbimorphae+Cuculiformes" = "Columbimorphae-Cuckoo",
      "Columbea+Otidiformae" = "TuBuCuckooColumbea",
      "Coraciimorphae-Coliiformes"="CPBTL",
      "Coraciimorphae" = "CPBTL-Coli",
      "Coraciimorphae+Strigiformes" = "CPBTL-Coli-Strigiformes",
      "Gruiformes+Charadriiformes" = "GruiformesCharadriiformes",
      "Gruiformes+CharadriiformesOpisthocomiformes" = "GruiformesCharadriiformesOpisthocomiformes",
      "Neoaves - Mirandornithes" = "MirandornithesBase",
      "Phaethontimorphae+Aequornithes" = "PhaethontimorphaeAequornithes",
      "Phaethontimorphae+Telluraves" = "PhaethontimorphaeTelluraves",
      "Apterygiformes+Casuariiformes+Rheiformes"= "RheiApterygiCasuari",
      "Apterygiformes+Casuariiformes+Tinamiformes" = "Rheiformesout",
      "Rheiformes+Tinamiformes"    = "RheiTinamu" ,
      "Strigiformes+Accipitriformes" = "StrigiformesAccipitriformes",
      "Strisores" = "Strisores" ,
      "Strisores+Aequornithes" = "StrisoresAequornithes",                     
      "Strisores+Aequornithes+Phaethontimorphae" = "StrisoresAequornithesPhaethontimorphae",
      "Strisores+Telluraves" ="StrisoresTelluraves")

renameT = function(x) {as.vector(sapply(x, function(y) {
  names(which(y ==   renameTaxa ))}))}

######################## quartets


w = read.csv('63K_trees.names_header.txt.xz',sep=" ")

g = read.csv('all-clade.stat.xz',sep=" ",header=F)
names(g) = c("Taxa","Ref","Gene","ref","shared", "ratio")
head(g)
nrow(g)
levels(g$Taxa)
g = merge(g,w,by="Gene")
g = merge(g,dcast(data=g,formula=Taxa~"mean",value.var = "ratio",fun.aggregate = mean),by.x = "Taxa",by.y = "Taxa")
nrow(g)
g$Chr = as.numeric(sub("chr","",sub("_.*","",g$Chromosome)))
g = g[order(g$Chr,g$ws),]
g$Chr = factor(ifelse(is.na(g$Chr), sub("chr","",sub("_.*","",g$Chromosome)),g$Chr),levels = c(1:28,"LGE22C19W28" ,"LGE64" ,"M", "Un", "W" ,"Z"))
g = g[g$Taxa != "Columbiformes",]
g = g[g$ref>0,]
g$p = 1:nrow(g)
#g$Chr = sub("chr","",sub("_.*","",g$Chromosome))
nrow(g)
co= cbind(dcast(g[,c("Chr","p")],Chr~.,fun.aggregate = min),dcast(g[,c("Chr","p")],Chr~.,fun.aggregate = max))
names(co)=c("Chr","start",".","end")
cl= co$Chr %in% c(1:10, 12, 16, 30, "Z", "W")
g = g[order(g$Taxa,g$Chr,g$ws),]

g$Taxa = factor(g$Taxa)
levels(g$Taxa) = renameTaxa

### Select group of interest below
groups=renameT(c("RheiApterygiCasuari","Rheiformesout","RheiTinamu")); n="Rhei"
groups=renameT(c("PhaethontimorphaeTelluraves","PhaethontimorphaeAequornithes")); n="Phaethontimorphae"
groups=renameT(c( "Strisores", "StrisoresAequornithes"      , "StrisoresAequornithesPhaethontimorphae" , "StrisoresTelluraves")); n="Capri"
groups=renameT(c("CPBTL", "CPBTL-Coli" ,  "CPBTL-Coli-Strigiformes" , "StrigiformesAccipitriformes")); n = "landbirdbase"
groups=levels(g$Taxa); n ="all";
groups=renameT(c("Columbimorphae","ColumbiformesOtidimorphae","Columbea","Columbiformes")); n="Columbea"



p1=ggplot(aes(x=p,y=ratio,color=Taxa),data= g[ g$Taxa %in% groups  & !g$Taxa=="Columbiformes" & g$ref>1,])+
  #data=g[grepl("apri",g$Taxa) ,])+
  #geom_point(alpha=0.2,size=0.1)+
  theme_classic()+
  #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
  #scale_color_brewer(palette = "Dark2",name="Clade")+
  #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
  #geom_ma(color="blue",n = 1000,linetype=1)+
  coord_cartesian(ylim=c(0.6,0.99))+
  scale_color_brewer(palette = "Dark2",name="")+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
  #geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="red",linetype=2)+
  scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  scale_y_continuous(name="Branch Quartet Support (BQS)",labels=percent,breaks=c(6:10)/10)+
  theme(legend.position = c(.2,.9),legend.direction = "vertical")+
  geom_segment(aes(y=0.602,yend=0.602,x=start,xend=end),color="black",size=0.75,data=co,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=-0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",size=2.3,data=co[co$end-co$start>15000,])+
  geom_vline(aes(xintercept=start),color="gray50",size=0.3,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray50",size=0.3,data=co, linetype=3)+
  geom_ma(n = 200,linetype=1,alpha=0.6)

p1
ggsave(paste("clades-nodots-one",n,".pdf",sep=""),width=7,height = 4)

  chr="4"
p2 = ggplot(aes(x=ws,y=ratio,color=Taxa),data= g[g$Chr == chr& g$Taxa %in% groups  & !g$Taxa=="Columbiformes"& g$ref>1,])+
    #data=g[grepl("apri",g$Taxa) ,])+
    geom_point(alpha=0.2,size=0.1)+
    theme_classic()+
    #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
    #scale_color_brewer(palette = "Dark2",name="Clade")+
    #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
    geom_ma(n = 50,linetype=1,alpha=0.75,size=0.8)+
    #geom_ma(color="blue",n = 1000,linetype=1)+
    coord_cartesian(ylim=c(0.6,1))+
    scale_color_brewer(palette = "Dark2",name="")+
    stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
    geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="gray50",linetype=3)+
    scale_x_continuous(name="Positions along chromosome 4")+
    scale_y_continuous(name="Branch Quartet Support (BQS)",labels=percent,breaks=c(6:10)/10)+
    theme(legend.position = "none",legend.direction = "horizontal")

p2
  ggsave(paste("chr",chr,"-ma",n,".pdf",sep=""),width=7,height = 4)
  

ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))
ggsave(paste("Combined-",n,".pdf",sep=""),width=7,height = 7)
 
  
  ggplot(aes(x=ws,y=ratio-mean,color=reorder(Taxa,ratio)),
         data= g[g$Taxa %in% groups & g$Chr %in% co[co$end-co$start>15000,"Chr"]& g$ref>1,])+
    #data=g[grepl("apri",g$Taxa) ,])+
    #geom_point(alpha=0.25,size=0.2)+
    theme_bw()+
    facet_wrap(~reorder(substr(Chr,0,2),Ref,FUN = function(X) -sum(X)),nrow=4,scales = "free_x",shrink = T)+
    #scale_color_brewer(palette = "Dark2",name="Clade")+
    #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
    geom_ma(n = 50,linetype=1,alpha=0.7)+
    #geom_ma(color="blue",n = 1000,linetype=1)+
    #coord_cartesian(ylim=c(0.5,0.99))+
    scale_color_brewer(palette = "Dark2",name="")+
    stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
    scale_x_continuous(name="regions (sorted by chicken coord)",breaks=c())+
    scale_y_continuous(name="Delta quartet agreement")+
    theme(legend.position = "bottom",panel.spacing = unit(0,"pt"))
ggsave(paste("facet-ma-delta",n,".pdf",sep=""),width=15,height = 9)
  

ggplot(aes(x=p,y=ratio),data=g[ g$Taxa %in% groups &g$ref>1,])+
  #data=g[grepl("apri",g$Taxa) ,])+
  #geom_point(alpha=0.25,size=0.2)+
  theme_classic2()+facet_wrap(~Taxa,nrow=7)+
  #scale_color_brewer(palette = "Dark2",name="Clade")+
  #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
  geom_ma(n = 200,linetype=1,alpha=0.5)+
  #geom_ma(color="blue",n = 1000,linetype=1)+
  coord_cartesian(ylim=c(0.6,0.99))+
  stat_summary(aes(x=1,yintercept=..y..),geom="hline",color="red",linetype=3)+
  geom_segment(aes(y=0.6,yend=0.6,x=start,xend=end,color=Chr),size=1,data=co,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  geom_text(aes(y=0.62,x=(start+end)/2,label=substr(Chr,0,2)),size=2.5,data=co[co$end-co$start>30000,])+
  scale_color_discrete(guide="none")+
  scale_x_continuous(name="regions (sorted by genome coord)",breaks =c())+
  scale_y_continuous(name="Branch Quartet Support (BQS)")
#scale_color_brewer(palette = "Set3",na.value="grey30")#+
#scale_x_continuous(lim=c(60000,65500))#
ggsave(paste("clades-nodots",n,".pdf",sep=""),width=10,height = 12)
ggsave(paste("clades-nodots",n,".png",sep=""),width=10,height = 12)

######################################## Monophyly
  
  w = read.csv('63K_trees.names_header.txt.xz',sep=" ")
  
  m = read.csv('all-monphyletic.txt.xz',sep="\t",header=F)
  m = cbind(m,1:nrow(w))
  names(m) = c("Taxa","mono", "Gene")
  head(m,n=30)
  nrow(m); 
  nrow(g);
  nrow(m)/nrow(w)
  levels(m$Taxa)
  dpres = read.csv('all-present.txt.xz', sep="\t", header=F)
  names(dpres) = c("Taxa","present", "Gene")
  head(dpres); nrow(dpres)
  m = merge(m,dpres,by=c("Taxa","Gene"))
  m = merge(m,w,by="Gene")
  nrow(m)
  
  m$Chr = as.numeric(sub("chr","",sub("_.*","",m$Chromosome)))
  m = m[order(m$Chr,m$ws),]
  m$Chr = factor(ifelse(is.na(m$Chr), sub("chr","",sub("_.*","",m$Chromosome)),m$Chr),levels = c(1:28,"LGE22C19W28" ,"LGE64" ,"M", "Un", "W" ,"Z"))
  m = m[m$Taxa != "Columbiformes",]
  #m = m[m$ref>0,]
  m$p = 1:nrow(m)
  
  nrow(m)
  co= cbind(dcast(m[,c("Chr","p")],Chr~.,fun.aggregate = min),dcast(m[,c("Chr","p")],Chr~.,fun.aggregate = max))
  names(co)=c("Chr","start",".","end")
  cl= co$Chr %in% c(1:10, 12, 16, 30, "Z", "W")
 
  m$Taxa = factor(m$Taxa)
  levels(m$Taxa) = renameTaxa
  
   ### Select group of interest below
  groups=renameT(c("RheiApterygiCasuari","Rheiformesout","RheiTinamu")); n="Rhei"
  groups=renameT(c("PhaethontimorphaeTelluraves","PhaethontimorphaeAequornithes")); n="Phaethontimorphae"
  groups=renameT(c( "Strisores", "StrisoresAequornithes"      , "StrisoresAequornithesPhaethontimorphae" , "StrisoresTelluraves")); n="Capri"
  groups=renameT(c("CPBTL","CPBTL-Coli",  "CPBTL-Coli-Strigiformes" , "StrigiformesAccipitriformes")); n = "landbirdbase"
  groups=unique(m$Taxa); n ="all";
  groups=renameT(c("Columbimorphae","ColumbiformesOtidimorphae","Columbea","Columbiformes")); n="Columbea"
  
p3 =  ggplot(aes(x=p,y=ifelse(is.na(mono),0,1),color=Taxa),data= m[ m$Taxa %in% groups & m$present >1 ,])+
    #data=g[grepl("apri",g$Taxa) ,])+
    #geom_point(alpha=0.2,size=0.1)+
    theme_classic()+
    #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
    #scale_color_brewer(palette = "Dark2",name="Clade")+
    #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
    #geom_ma(color="blue",n = 1000,linetype=1)+
    #coord_cartesian(ylim=c(0.645,0.99))+
    coord_cartesian(ylim=c(0,0.99))+
    scale_color_brewer(palette = "Dark2",name="")+
    stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
    #geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="red",linetype=2)+
    scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
    scale_y_continuous(name="Monophyletic loci",labels=percent)+
    theme(legend.position = c(.2,.9),legend.direction = "vertical")+
    geom_segment(aes(y=-0.02,yend=-0.02,x=start,xend=end),color="black",size=0.75,data=co,
                 arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
    #geom_text(aes(y=-0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",size=2.3,data=co[co$end-co$start>15000,])+
    geom_vline(aes(xintercept=start),color="gray50",size=0.3,data=co, linetype=3)+
    geom_vline(aes(xintercept=end),color="gray50",size=0.3,data=co, linetype=3)+
    geom_ma(n = 200,linetype=1,alpha=0.6)
 
  p3
  ggsave(paste("monophyly-nodots-one",n,".pdf",sep=""),width=7,height = 4)
  
  
  chr="4"
  p4 = ggplot(aes(x=we,y=ifelse(is.na(mono),0,1),color=Taxa),data= m[ m$Chr == chr & m$Taxa %in% groups & m$present >1 ,])+
    #data=g[grepl("apri",g$Taxa) ,])+
    geom_point(alpha=0.2,size=0.1)+
    theme_classic()+
    #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
    #scale_color_brewer(palette = "Dark2",name="Clade")+
    #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
    geom_ma(n = 50,linetype=1,alpha=0.75,size=0.8)+
    #geom_ma(color="blue",n = 1000,linetype=1)+
    coord_cartesian(ylim=c(0,0.99))+
    scale_color_brewer(palette = "Dark2",name="")+
    #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
    geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="gray50",linetype=3)+
    scale_x_continuous(name="Positions along chromosome 4")+
    scale_y_continuous(name="Monophyletic loci",labels=percent)+
    theme(legend.position = "none",legend.direction = "horizontal")
  
  p4
  ggsave(paste("monophyly-chr",chr,"-ma",n,".pdf",sep=""),width=7,height = 4)
  
  
  
  ggplot(aes(x=p,y=ifelse(is.na(mono),0,1)),data=m[ m$Taxa %in% groups ,])+
    #data=g[grepl("apri",g$Taxa) ,])+
    #geom_point(alpha=0.25,size=0.2)+
    theme_classic2()+facet_wrap(~Taxa,nrow=7)+
    #scale_color_brewer(palette = "Dark2",name="Clade")+
    #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
    geom_ma(n = 200,linetype=1,alpha=0.5)+
    #geom_ma(color="blue",n = 1000,linetype=1)+
    coord_cartesian(ylim=c(-0.03,0.99))+
    stat_summary(aes(x=1,yintercept=..y..),geom="hline",color="red",linetype=3)+
    geom_segment(aes(y=-0.03,yend=-0.03,x=start,xend=end,color=Chr),size=0.5,data=co,
                 arrow=arrow(ends = "both",type = "closed",length =unit(1,"pt")))+
    geom_text(aes(y=-0.04,x=(start+end)/2,label=substr(Chr,0,2)),size=2.5,data=co[co$end-co$start>30000,])+
    scale_color_discrete(guide="none")+
    scale_x_continuous(name="regions (sorted by genome coord)",breaks =c())+
    scale_y_continuous(name="Monothyletic loci")
  #scale_color_brewer(palette = "Set3",na.value="grey30")#+
  #scale_x_continuous(lim=c(60000,65500))#
  ggsave(paste("monophyly-nodots",n,".pdf",sep=""),width=10,height = 12)
  ggsave(paste("monophyly-nodots",n,".png",sep=""),width=10,height = 12)
  
  
  
  ############ Combine
  
  ggarrange(p3, p4, ncol = 1, labels = c("A", "B"))
  ggsave(paste("Combined-Monophyly-",n,".pdf",sep=""),width=7,height = 7)
  
  p3=p3+theme(legend.position = "none")
  ggarrange(p1, p2, p3, p4, ncol = 2, nrow=2, labels = c("A", "B" , "C", "D"))
  ggsave(paste("Combined-both-new",n,".eps",sep=""),width=12,height = 9, device=cairo_ps)



