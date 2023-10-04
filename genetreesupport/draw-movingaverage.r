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
names(g) = c("Taxa","Ref","Gene","ref","difference", "ratio")
head(g)
nrow(g)
levels(g$Taxa)
g = merge(g,w,by="Gene")
nrow(g)
g$Chr = as.numeric(sub("chr","",sub("_.*","",g$Chromosome)))
g = g[order(g$Chr,g$ws),]
g$Chr = factor(ifelse(is.na(g$Chr), sub("chr","",sub("_.*","",g$Chromosome)),g$Chr),levels = c(1:28,"LGE22C19W28" ,"LGE64" ,"M", "Un", "W" ,"Z"))
g = g[g$Taxa != "Columbiformes",]
g = g[g$ref>0,]
g$p = 1:nrow(g)
g = merge(g,dcast(data=g,formula=Taxa~"mean",value.var = "ratio",fun.aggregate = mean),by.x = "Taxa",by.y = "Taxa")
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


ggplot(aes(x=p,y=ratio-mean,group=reorder(Taxa,mean),color=reorder(Taxa,mean)),
       data=g[ g$Taxa %in% groups &g$ref>1,])+
  theme_classic2()+
  #facet_wrap(~Taxa)+
  geom_hline(yintercept=0,linetype=3)+
  geom_ma(n = 200,linetype=1,alpha=1/3)+
  scale_color_viridis_d(guide="none",option = "C")+
  scale_x_continuous(name="regions (sorted by genome coord)",breaks =c())+
  scale_y_continuous(name="BQS above mean")
ggsave(paste("relative-BQS",n,".png",sep=""),width=16,height = 15)

ggplot(aes(x=p,y=ratio,color=Taxa),data= g[ g$Taxa %in% groups  & !g$Taxa=="Columbiformes" & g$ref>1,])+
  #data=g[grepl("apri",g$Taxa) ,])+
  #geom_point(alpha=0.2,size=0.1)+
  theme_classic()+
  #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
  #scale_color_brewer(palette = "Dark2",name="Clade")+
  #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
  #geom_ma(color="blue",n = 1000,linetype=1)+
  #coord_cartesian(ylim=c(0.6,0.99))+
  scale_color_brewer(palette = "Dark2",name="")+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
  #geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="red",linetype=2)+
  scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  #scale_y_continuous(name="Branch Quartet Support (BQS)",labels=percent,breaks=c(6:10)/10)+
  theme(legend.position = c(.2,.9),legend.direction = "vertical")+
  geom_segment(aes(y=0.602,yend=0.602,x=start,xend=end),color="black",size=0.75,data=co,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=-0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",size=2.3,data=co[co$end-co$start>15000,])+
  geom_vline(aes(xintercept=start),color="gray50",size=0.3,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray50",size=0.3,data=co, linetype=3)+
  geom_ma(n = 200,linetype=1,alpha=0.6)


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

######################################33
  
r = read.csv('clade-rec.stat.xz',sep=" ",h=F)
tail(r)



r = dcast( data=r,V4+V1+V2~V3)
head(r)
r$T1 = r$`1` - r$`4` 
r$T2 = r$`2` - r$`4` 
r$T3 = r$`3` - r$`4` 


names(r)[1:3] = c("Gene","C1","C2")
r$clade = factor(interaction(r$C1,r$C2,sep="."))
levels(r$clade) = list( 
  "Columbea" =  "Columbimorphae.Phoenicopteriformes"   ,
  "N62" = "Columbiformes.OtherColumbimorphae",
  "N61" = "Otidimorphae.Columbimorphae" , 
  "N60" = "ElementavesTelluraves.Columbaves" , 
  "N45" =  "Telluraves.Elementaves"   , 
  "N5" =  "Cuculiformes.Otidiformes"   ,  
  "N27" ="CPBTL-Coli.StrigiformesAccipitriformes"   , 
  "N23" = "CPBTL.Coliiformes" ,        
  "N26" ="Accipitriformes.Strigiformes" ,                                                    
  "N44" =  "StrisoresAequornithesPhaethontimorphae.GruiformesCharadriiformesOpisthocomiformes",
  "N39" = "PhaethontimorphaeAequornithes.Strisores"  ,
  "N37" = "Aequornithes.Phaethontimorphae"  ,  
  "N43" =  "GruiformesCharadriiformes.Opisthocomiformes"  ,  
  "N42" =  "Charadriiformes.Gruiformes"   ,     
  "N49" = "Casuariiformes.Apterygiformes",
  "N53" = "Tinamiformes.Rheiformes"    
)

r$`C1,C2|S,O` = (r$T2+r$T3-r$T1)/(r$T1+r$T2+r$T3) 
r$`C1,S|C2,O` = (r$T1+r$T3-r$T2)/(r$T1+r$T2+r$T3) 
r$`C2,S|C1,O` = (r$T1+r$T2-r$T3)/(r$T1+r$T2+r$T3)  


tail(r)
r = merge(r,w,by="Gene")
nrow(r)
r$Chr = as.numeric(sub("chr","",sub("_.*","",r$Chromosome)))
r = r[order(r$Chr,r$ws),]
r$Chr = factor(ifelse(is.na(r$Chr), sub("chr","",sub("_.*","",r$Chromosome)),r$Chr),levels = c(1:28,"LGE22C19W28" ,"LGE64" ,"M", "Un", "W" ,"Z"))
r = r[!is.na(r$`C1,C2|S,O`),]
r$p = 1:nrow(r)
#r$Chr = sub("chr","",sub("_.*","",r$Chromosome))
nrow(r)
co= cbind(dcast(r[,c("Chr","p")],Chr~.,fun.aggregate = min),dcast(r[,c("Chr","p")],Chr~.,fun.aggregate = max))
names(co)=c("Chr","start",".","end")
cl= co$Chr %in% c(1:8, 10, 12, 16,24, "Z")
r = r[order(r$Chr,r$ws),]




head(r)
ggplot(aes(x=p,y=value,color=clade, 
           linetype=variable,size=variable,alpha=variable),
       data= melt(r[r$Chr %in% c(1:15,17:24,26:28,"Z"),c(1:3,21,11:13,14:20)],measure.vars = 5:7))+
  theme_classic()+#facet_wrap(~interaction(C1,C2,sep="/"))+
  geom_hline(color="grey30",yintercept = 1/3)+
  scale_color_viridis_d(name="",option = "B")+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline")+
  scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=1)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=1)+
  geom_ma(aes(group=interaction(Chr,variable,C1,C2)),n = 100)+
  scale_size_manual(name="topology",values=c(1,0.5,0.5)*0.75)+
  scale_linetype_manual(name="topology",values=c(1,2,3))+
  scale_alpha_manual(name="topology",values=c(0.6,0.33,0.33))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))
  #facet_wrap(~grepl("Col",paste(C1,C2)))
ggsave("all-clades-QQS.png",width=14,height = 9) 


ggplot(aes(x=p,y=value,#color=interaction(C1,C2,sep="/"), 
           color=variable,size=variable,alpha=variable),
       data= melt(r[r$Chr %in% c(1:15,17:24,26:28,"Z"),c(1:3,21,11:13,14:20)],measure.vars = 5:7))+
  theme_classic()+facet_wrap(~clade,ncol=4)+
  geom_hline(color="grey80",yintercept = 1/3,linetype=3)+
  scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,variable,C1,C2)),n = 100,linetype=1)+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=1,size=1)+
  scale_size_manual(name="topology",values=c(0.85,0.5,0.5)*0.8)+
  scale_color_manual(name="topology",values=c("#107030","#C05050","#4080E0"))+
  scale_alpha_manual(name="topology",values=c(0.6,0.33,0.33))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))
ggsave("all-clades-QQS-facet.png",width=14,height = 15) 


ggplot(aes(x=p,y=value,#color=interaction(C1,C2,sep="/"), 
           color=variable,size=variable,alpha=variable),
       data= melt(r[r$Chr %in% c(1:15,17:24,26:28,"Z") &r$clade =="N53",c(1:3,21,11:13,14:20)],measure.vars = 5:7))+
  theme_classic()+
  geom_hline(color="grey80",yintercept = 1/3,linetype=3)+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=1)+
  scale_x_continuous(name="Loci",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,variable,C1,C2)),n = 100,linetype=1)+
  scale_size_manual(name="topology",values=c(0.85,0.5,0.5)*0.8, labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_color_manual(name="topology",values=c("#107030","#C05050","#4080E0"), labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_alpha_manual(name="topology",values=c(0.6,0.33,0.33), labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))
ggsave("all-clades-QQS-tinamus-rhea.pdf",width=7,height = 5) 

ggplot(aes(x=stats::filter(`C1,C2|S,O`-mean(`C1,C2|S,O`), filter = rep(1/20,20) , method="convolution")/
            sd(stats::filter(`C1,C2|S,O`,filter = rep(1/20,20) , method="convolution"),na.rm = T)),
       data= r[r$Chr %in% c(1:15,17:24,26:28,"Z") & r$clade =="N44",c(1:3,21,11:13,14:20)])+
  theme_classic()+geom_density()+
  #scale_y_log10()+
  stat_function(fun=dnorm,color="red")


ggplot(aes(x=p.adjust(twowayp(stats::filter(`C1,C2|S,O`-mean(`C1,C2|S,O`), filter = rep(1/20,20) , method="convolution")/
    sd(stats::filter(`C1,C2|S,O`,filter = rep(1/20,20) , method="convolution"),na.rm = T)),method="BH")),
       data= r[r$Chr %in% c(1:15,17:24,26:28,"Z") & r$clade =="N44",c(1:3,21,11:13,14:20)])+
  theme_classic()+geom_histogram()+
  scale_x_log10("BH corrected p-values (N23)")+scale_y_log10()+geom_vline(xintercept = 0.05,color="red")

ggplot(aes(x=p.adjust(1-pnorm(filter(`C1,C2|S,O`-mean(`C1,C2|S,O`), filter = rep(1/20,20) , method="convolution")/
            sd(filter(`C1,C2|S,O`,filter = rep(1/20,20) , method="convolution"),na.rm = T)),method="BH")),
       data= r[r$Chr %in% c(1:15,17:24,26:28,"Z") & r$clade =="Columbea",c(1:3,21,11:13,14:20)])+
  theme_classic()+geom_histogram()+
  scale_x_log10("BH corrected p-values (Columbea)")+scale_y_log10()+geom_vline(xintercept = 0.05,color="red")


library(dplyr)
pvalues = function(x, chr,window = 21 ) {
  #mg = which(!is.na(x))
  #x = x[mg]
  #chr = chr[!mg]
  bounds = which(chr[-1] != chr[-length(chr)])
  bounds = 0:(length(bounds)-1)+bounds
  for (i in bounds) x <- c(head(x, i-1), NA, tail(x, -(i-1)))
  b2=which(!is.na(x))
  twowayp = function (x) pmin(1-pnorm(x) , pnorm(x))
  #print(mean(x,na.rm = T))
  ps=p.adjust(twowayp(stats::filter(x-mean(x,na.rm = T), filter = rep(1/window,window) , method="convolution")/
          sd(stats::filter(x,filter = rep(1/window,window) , method="convolution"),na.rm = T)),method="BH")
  #print(c(length(ps), length(ps[b2])))
  ps[b2]
}



r[r$Chr %in% c(1:30),] %>%  group_by(clade) %>%
  #filter(clade=="N44") %>%
  mutate(pvalue = pvalues(`C1,C2|S,O`,Chr,20))  %>% 
  group_by(clade) %>%
  summarise(rejected = sum(pvalue<0.01,na.rm = T), n= sum(pvalue<2,na.rm = T))


r %>%  group_by(clade) %>%
  mutate(pvalue = pvalues(`C1,C2|S,O`))  %>% 
  ggplot(aes(x=pvalue))+
  theme_classic()+geom_histogram()+facet_wrap(~clade)+
  scale_x_log10("BH corrected p-values")+scale_y_log10()+geom_vline(xintercept = 0.05,color="red")
  
  
pr = r[r$Chr %in% c(1:30),] %>%  group_by(clade) %>%
  mutate(pvalue = pvalues(`C1,C2|S,O`,Chr,20))

cl= co$Chr %in% c(1:8, 10, 12, 15,20,28)
ggplot(aes(x=p,y=pvalue,#color=interaction(C1,C2,sep="/"), 
             #color=cut(pvalue,c(0,0.0001,0.01,0.05,1)),
           color=clade,
             alpha=abs(log10(pvalue)),
              size=cut(pvalue,c(0,0.0001,0.01,0.05,0.1,1))),
       data=pr)+
  theme_classic()+
  scale_color_viridis_d(name="",option = "B")+
  geom_hline(color="grey80",yintercept = c(0,0.0001,0.01,0.05),linetype=1)+
   scale_x_continuous(name="Genomic position",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(.19,.15))+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_point()+
  scale_alpha(guide="none")+
  scale_size_manual(name="significance",values=c(0.8,0.7,0.5,0.3,0.1)*0.7, guide="none")+
  #scale_color_manual(name="significance",values=c("#107030","#C05050","#4080E0","grey80"),)+
  #scale_alpha_manual(name="significance",values=c(0.8,0.6,0.3,0.10,0.01),guide="none" )+
  scale_y_continuous(name="p-value",trans = "log10")+
  guides(color=guide_legend(nrow=4,byrow=TRUE))
ggsave("manhattan.png",width=9,height = 6)

ggplot(aes(x=p,y=-log10(pvalue),#color=interaction(C1,C2,sep="/"), 
           #color=cut(pvalue,c(0,0.0001,0.01,0.05,1)),
           color=clade,
           alpha=abs(log10(pvalue)),
           size=cut(pvalue,c(0,0.0001,0.01,0.05,0.1,1))),
       data=pr)+
  theme_classic()+
  scale_color_viridis_d(name="",option = "C")+
  geom_hline(color="grey80",yintercept = -log10(c(0,0.01)),linetype=1)+
  scale_x_continuous(name="Chromosome",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(.19,.83))+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_point()+
  scale_alpha(guide="none")+
  coord_cartesian(xlim=c(30000,876000))+
  scale_size_manual(name="significance",values=c(0.8,0.7,0.5,0.3,0.1)*0.7, guide="none")+
  #scale_color_manual(name="significance",values=c("#107030","#C05050","#4080E0","grey80"),)+
  #scale_alpha_manual(name="significance",values=c(0.8,0.6,0.3,0.10,0.01),guide="none" )+
  scale_y_continuous(name=expression(-log[10](P)))+
  guides(color=guide_legend(nrow=4,byrow=TRUE))
ggsave("manhattan-2.png",width=9,height = 5.2)



ggplot(aes(x=p,y=-log10(pvalue),#color=interaction(C1,C2,sep="/"), 
           #color=cut(pvalue,c(0,0.0001,0.01,0.05,1)),
           color=clade,
           alpha=abs(log10(pvalue)),
           size=cut(pvalue,c(0,0.0001,0.01,0.05,0.1,1))),
       data=pr[])+
  theme_classic()+
  scale_color_manual(name="",guide="none",
                     values=c(rep("#FC0005",2),"#EF712B","#F59F1B",rep("#162EC4",12)))+
  geom_hline(color="grey80",yintercept = -log10(c(0,0.01)),linetype=1)+
  scale_x_continuous(name="",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(.19,.83),axis.title.x = element_blank())+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co[1:28,], linetype=1)+
  geom_point()+
  scale_alpha(guide="none")+
  coord_cartesian(xlim=c(30000,876000))+
  scale_size_manual(name="significance",values=c(0.8,0.7,0.5,0.3,0.1)*0.7, guide="none")+
  #scale_color_manual(name="significance",values=c("#107030","#C05050","#4080E0","grey80"),)+
  #scale_alpha_manual(name="significance",values=c(0.8,0.6,0.3,0.10,0.01),guide="none" )+
  scale_y_continuous(name=expression(-log[10](P)))
ggsave("manhattan-3.png",width=9,height = 5.2)

