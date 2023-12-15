require(ggplot2); require(scales); require(reshape2);library(tidyquant);library(zoo)
library("ggpubr"); require(imager)
library(dplyr)
library(cowplot)
library(magick)
library(forcats)


renameTaxa = list ("Columbea (J2014)" = "Columbea" ,
      "Columbimorphae+Otidimorphae (S2023)"  = "ColumbiformesOtidimorphae"     ,
      "Columbimorphae (J2014 and S2023)" ="Columbimorphae",                            
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

renameTaxa2 = list ("Columbea (J2014)" = "Columbea" ,
                   "N61"  = "ColumbiformesOtidimorphae"     ,
                   "N62" ="Columbimorphae",             
                   "N60" = "MirandornithesBase",
                   "N45" =  "ElementavesTelluraves"   , 
                   "N5" =  "CuculiformesOtidiformes"   ,
                   "N23" = "CPBTL-Coli",
                   "N27" ="CPBTL-Coli-StrigiformesAccipitriformes" ,
                   "N26" = "StrigiformesAccipitriformes", 
                   "N42" = "GruiformesCharadriiformes",
                   "N44" =  "StrisoresAequornithesPhaethontimorphae.GruiformesCharadriiformesOpisthocomiformes",
                   "N39" = "StrisoresAequornithesPhaethontimorphae",
                   "N37" = "PhaethontimorphaeAequornithes",
                   "N43" =  "GruiformesCharadriiformes.Opisthocomiformes"  ,  
                   "N42" =  "Charadriiformes.Gruiformes"   ,     
                   "N49" = "Casuariiformes.Apterygiformes",
                   "N53"    = "RheiTinamu" ,
                   "REM" = "PhaethontimorphaeTelluraves",
                   "REM"= "RheiApterygiCasuari",
                   "REM" = "Rheiformesout",
                   "REM" = "Columbimorphae-Cuckoo",
                   "REM" = "TuBuCuckooColumbea",
                   "REM"="CPBTL",
                   "REM" = "CPBTL-Coli-Strigiformes",
                   "REM" = "GruiformesCharadriiformesOpisthocomiformes",
                   "REM" = "Strisores" ,
                   "REM" = "StrisoresAequornithes",                     
                   "REM" ="StrisoresTelluraves")

renameT = function(x) {as.vector(sapply(x, function(y) {
  names(which(y ==   renameTaxa ))}))}

renameT2 = function(x) {as.vector(sapply(x, function(y) {
  names(which(y ==   renameTaxa2 ))}))}

###################################### BQS

w = read.csv('63K_trees.names_header.txt.xz',sep=" ")
head(w); nrow(w);

g = read.csv('all-clade.stat.xz',sep=" ",header=F)
g = rbind(g, read.csv('new-all-clade.stat.xz',sep=" ",header=F))

names(g) = c("Taxa","ignore","Gene","ref","difference", "ratio")
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
levels(g$Taxa) = renameTaxa2

### Select group of interest below
groups=renameT2(c("RheiApterygiCasuari","Rheiformesout","RheiTinamu")); n="Rhei"
groups=renameT2(c("PhaethontimorphaeTelluraves","PhaethontimorphaeAequornithes")); n="Phaethontimorphae"
groups=renameT2(c( "Strisores", "StrisoresAequornithes"      , "StrisoresAequornithesPhaethontimorphae" , "StrisoresTelluraves")); n="Capri"
groups=renameT2(c("CPBTL", "CPBTL-Coli" ,  "CPBTL-Coli-Strigiformes" , "StrigiformesAccipitriformes")); n = "landbirdbase"


head(g)



groups=renameT2(c("Columbimorphae","ColumbiformesOtidimorphae","Columbea","Columbiformes")); n="Columbea"

p1=ggplot(aes(x=p,y=ratio,color=Taxa),data= g[ g$Taxa %in% groups  & !g$Taxa=="Columbiformes" & g$ref>1,])+
  #data=g[grepl("apri",g$Taxa) ,])+
  #geom_point(alpha=0.2,size=0.1)+
  theme_classic()+
  #facet_wrap(~reorder(Taxa,ratio),nrow=1)+
  #scale_color_brewer(palette = "Dark2",name="Clade")+
  #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
  #geom_ma(color="blue",n = 1000,linetype=1)+
  coord_cartesian(ylim=c(0.65,1))+
  scale_color_brewer(palette = "Dark2",name="")+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
  #geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="red",linetype=2)+
  scale_x_continuous(name="Chromosomes",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  scale_y_continuous(name="Branch Quartet Support (BQS)",labels=percent,breaks=c(6:10)/10)+
  theme(legend.position = c(.24,.92),legend.direction = "vertical",legend.margin = margin(0,0,0,0),legend.box.margin =  margin(0,0,0,0))+
  geom_segment(aes(y=0.634,yend=0.634,x=start,xend=end),color="black",size=0.75,data=co,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.7,"pt")))+
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
    coord_cartesian(ylim=c(0.65,1))+
    scale_color_brewer(palette = "Dark2",name="")+
    stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",linetype=3)+
    geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="gray50",linetype=3)+
    scale_x_continuous(name="Positions along chromosome 4")+
    scale_y_continuous(name="Branch Quartet Support (BQS)",labels=percent,breaks=c(6:10)/10)+
    theme(legend.position = "none",legend.direction = "horizontal")
p2
  ggsave(paste("chr",chr,"-ma",n,".pdf",sep=""),width=7,height = 4)
  

ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
ggsave(paste("Combined-",n,".pdf",sep=""),width=14,height = 5)
ggsave(paste("Combined-BQS",n,".png",sep=""),width=14,height = 5)

  

groups=setdiff(levels(g$Taxa),"REM"); n ="all-new";

ggplot(aes(x=p,y=ratio),data=g[ g$Taxa %in% groups &g$ref>1,])+
  #data=g[grepl("apri",g$Taxa) ,])+
  #geom_point(alpha=0.25,size=0.2)+
  theme_classic2()+facet_wrap(~Taxa,nrow=4)+
  #scale_color_brewer(palette = "Dark2",name="Clade")+
  #geom_ma(color="green",n = 50,linetype=1,alpha=0.45)+
  geom_ma(n = 200,linetype=1,alpha=0.5)+
  #geom_ma(color="blue",n = 1000,linetype=1)+
  coord_cartesian(ylim=c(0.6,0.99))+
  stat_summary(aes(x=1,yintercept=..y..),geom="hline",color="red",linetype=3)+
  #geom_segment(aes(y=0.6,yend=0.6,x=start,xend=end),size=1,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.62,x=(start+end)/2,label=substr(Chr,0,2)),size=2.5,data=co[co$end-co$start>30000,])+
  geom_text(aes(y=0.60,x=(start+end)/2,label=substr(Chr,0,2)),color="black",alpha=1,
            size=2.5,data=co[co$end-co$start>13000,])+
  geom_segment(aes(y=-0.615,yend=-0.015,x=start,xend=end),color="black",size=0.75,data=co,alpha=1,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.7,"pt")))+
  scale_color_discrete(guide="none")+
  scale_x_continuous(name="Chromosomes",breaks =c())+
  scale_y_continuous(name="Branch Quartet Support (BQS)")
#scale_color_brewer(palette = "Set3",na.value="grey30")#+
#scale_x_continuous(lim=c(60000,65500))#
ggsave(paste("clades-nodots",n,".pdf",sep=""),width=10,height = 12)
#ggsave(paste("clades-nodots",n,".png",sep=""),width=10,height = 12)


### Not used:
ggplot(aes(x=p,y=ratio-mean,group=reorder(Taxa,mean),color=reorder(Taxa,mean)),
       data=g[ g$Taxa %in% groups &g$ref>1,])+
  theme_classic2()+
  #facet_wrap(~Taxa)+
  geom_hline(yintercept=0,linetype=3)+
  geom_ma(n = 200,linetype=1,alpha=1/3)+
  scale_color_viridis_d(guide="none",option = "C")+
  scale_x_continuous(name="regions (sorted by genome coord)",breaks =c())+
  scale_y_continuous(name="BQS above mean")
ggsave(paste("relative-BQS",n,".png",sep=""),width=12,height = 10)

ggplot(aes(x=p,y=ratio-mean,group=reorder(Taxa,mean),color=reorder(Taxa,mean)),
       data=g[ g$Taxa %in% groups &g$ref>1,])+
  theme_classic2()+
  facet_wrap(~Taxa)+
  geom_hline(yintercept=0,linetype=3)+
  geom_ma(n = 200,linetype=1,alpha=1/3)+
  scale_color_viridis_d(guide="none",option = "C")+
  scale_x_continuous(name="regions (sorted by genome coord)",breaks =c())+
  scale_y_continuous(name="BQS above mean")
ggsave(paste("relative-BQS-facet",n,".png",sep=""),width=16,height = 16)



###################################### QQS
  

r = read.csv('clade-rec.stat.xz',sep=" ",h=F)
tail(r)

r = dcast( data=r,V4+V1+V2~V3)
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

tail(r)
r$T1 = r$`1` - r$`4` 
r$T2 = r$`2` - r$`4` 
r$T3 = r$`3` - r$`4` 


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



print( r %>% group_by(Chr,clade) %>% summarise(n=n()) %>% group_by(Chr)  %>% summarise(a=min(n)), n=50)

r3 = melt(r[r$Chr %in% c(1:15,17:21,23,24,"Z"),c(
  "Gene"  , "p", "clade",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O","Chromosome", "ws"        
  ,"we"      ,   "size"     ,  "wstart" ,    "Chr" )],
  measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" ))

print(r3 %>% group_by(variable,clade) %>% summarize(qqsm = mean(value,na.rm=T)),n=100)

print(r %>% group_by(clade) %>% summarize(qqsm = mean(T1+T2+T3,na.rm=T)),n=100)


ggplot(aes(x=p,y=value,color=clade, 
           linetype=variable,size=variable,alpha=variable),
       data= r3)+
  theme_classic()+#facet_wrap(~interaction(C1,C2,sep="/"))+
  geom_hline(color="grey30",yintercept = 1/3)+
  scale_color_viridis_d(name="",option = "B")+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline")+
  scale_x_continuous(name="Chromosomes",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(0.2,0.8),legend.margin = margin(0,0,0,0),legend.box.margin =  margin(0,0,0,0))+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=1)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=1)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 200)+
  scale_size_manual(name="topology",values=c(1,0.5,0.5)*0.75)+
  scale_linetype_manual(name="topology",values=c(1,2,3))+
  scale_alpha_manual(name="topology",values=c(0.6,0.33,0.33))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))+
  guides(linetype=guide_legend(nrow=3,byrow=TRUE),
         size=guide_legend(nrow=3,byrow=TRUE),
         alpha=guide_legend(nrow=3,byrow=TRUE),
         color=guide_legend(nrow=3,byrow=TRUE))
  #facet_wrap(~grepl("Col",paste(C1,C2)))
#ggsave("all-clades-QQS.pdf",width=14,height = 9) 
ggsave("all-clades-QQS.png",width=12,height = 4.8) 


ggplot(aes(x=p,y=value,color=variable, 
           size=variable,alpha=variable),
       data= r3)+
  theme_classic()+facet_wrap(~clade,ncol=4)+
  geom_hline(color="grey30",yintercept = 1/3)+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline")+
  scale_x_continuous(name="Chromosomes",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 200,linetype=1)+
  scale_size_manual(name="topology",values=c(1,0.5,0.5)*0.5)+
  scale_color_manual(name="topology",values=c("black","#40D080","#D04080"))+
  scale_alpha_manual(name="topology",values=c(0.8,0.5,0.5))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))+
  guides(size=guide_legend(nrow=1,byrow=TRUE),
         alpha=guide_legend(nrow=1,byrow=TRUE),
         color=guide_legend(nrow=1,byrow=TRUE))
#facet_wrap(~grepl("Col",paste(C1,C2)))
ggsave("all-clades-QQS-facet.pdf",width=14,height = 14) 
ggsave("all-clades-QQS-facet.png",width=14,height = 14) 

r3[r3$clade %in% c("N61", "N62" ,"Columbea"),] %>% 
  mutate(clade = fct_recode(clade, 
                            "S2023:\n(Columbaves)(Elementaves+Telluraves)|\n(Mirandornithes),(Outgroups)" = "N61",
                            "J2014:\n(Columbimorphae)(Mirandornithes)|\n(Passerea),(Outgroups)" = "Columbea",
                            "J2014 nad S2023:\n(Columbiformes)(Mesitornithiformes+Pterocliformes)|\n(Otidimorphae),(Other birds)" = "N62")) %>%
ggplot(aes(x=p,y=value,color=variable, 
           size=variable,alpha=variable))+
  theme_classic()+facet_wrap(~clade,ncol=3)+
  geom_hline(color="#33a02c",yintercept = 1/3,linetype=3)+
  stat_summary(aes(x=min(p),yintercept=..y..),geom="hline")+
  scale_x_continuous(name="Chromosomes",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 100,linetype=1)+
  scale_size_manual(name="topology",values=c(1,1,0.5)*0.6)+
  scale_color_manual(name="topology",values=c("grey20","#a6cee3","#fb9a99"))+
  #scale_color_brewer(name="topology",palette = "Dark2")+
  scale_alpha_manual(name="topology",values=c(0.8,0.8,0.5))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))+
  guides(size=guide_legend(nrow=1,byrow=TRUE),
         alpha=guide_legend(nrow=1,byrow=TRUE),
         color=guide_legend(nrow=1,byrow=TRUE))
ggsave("all-clades-QQS-Columbea-facet.pdf",width=15,height = 5) 
ggsave("all-clades-QQS-Columbea-facet.png",width=15,height = 5) 


MAS=200
p6=melt(r[r$Chr %in% 
              (r %>% filter(clade %in%  c("N61","Columbea") ) %>% 
                 group_by(Chr,clade) %>% 
                 summarise(n=n()) %>% group_by(Chr)  %>% 
                 summarise(a=min(n)) %>% filter(a>MAS))$Chr
            & r$clade %in% c("N61","Columbea") ,c(
    "Gene"  , "p", "clade",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O","Chromosome", "ws"        
    ,"we"      ,   "size"     ,  "wstart" ,    "Chr" )],
    measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" )) %>% 
    filter(variable == "C1,C2|S,O") %>% 
    mutate( clade = fct_recode(clade, "S2023" = "N61", 
                              "J2014 (Columbea)" = "Columbea"
                            #"3rd alt" = "C2,S|C1,O"
                                  )
            ) %>%
 ggplot(aes(x=p,y=value,color=clade))+
  theme_classic()+
  geom_hline(color="grey20",yintercept = 1/3,linetype=3)+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",size=0.8)+
  scale_x_continuous(name="Chromosomes",breaks = c() )+# =co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(.15,.92),legend.direction = "vertical")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  geom_text(aes(y=0.00,x=(start+end)/2,label=substr(Chr,0,2)),color="black",alpha=1,
            size=2.5,data=co[co$end-co$start>13000,])+
  geom_segment(aes(y=-0.015,yend=-0.015,x=start,xend=end),color="black",size=0.75,data=co,alpha=1,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.7,"pt")))+
  geom_vline(aes(xintercept=start),color="gray60",size=0.2,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray60",size=0.2,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = MAS,linetype=1,alpha=0.6)+
  scale_size_manual(name="",values=c(1,1,0.5)*0.8)+
  scale_color_manual(name="",values=c("#1b9e77","#d95f02","#e7298a"))+
  scale_alpha_manual(name="",values=c(0.6,0.6,0.4))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))+
  coord_cartesian(ylim=c(0.025,1))

p6
p7 = melt(r[r$Chr ==4    & r$clade %in% c("N61","Columbea") ,c(
              "Gene"  , "p", "clade",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O","Chromosome", "ws"        
              ,"we"      ,   "size"     ,  "wstart" ,    "Chr" )],
          measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" )) %>% 
  filter(variable == "C1,C2|S,O") %>% 
  mutate( clade = fct_recode(clade, "S2023" = "N61", 
                             "J2014 (Columbea)" = "Columbea"
                             #"3rd alt" = "C2,S|C1,O"
  )
  ) %>%
  ggplot(aes(x=we,y=value,color=clade, 
             size=variable))+
  theme_classic()+
  geom_hline(color="grey20",yintercept = 1/3,linetype=3)+
  scale_x_continuous(name="Positions along chromosome 4")+# =co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "none")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 50,linetype=1, alpha=0.7)+
  scale_size_manual(name="",values=c(1,1,0.5)*0.8)+
  scale_color_manual(name="",values=c("#1b9e77","#d95f02","#e7298a"))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))+
  coord_cartesian(ylim=c(0.025,1))+
  geom_point(alpha=0.2,size=0.1)+
  geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="gray50",linetype=3)
  #geom_vline(xintercept = c(488932, 497238 ,497757,498734,509454,522053),color="gray50",linetype=3)
  
p7


########## To compute monophyly, we simply test to see if the support for a topology == 1
p8=melt(r[r$Chr %in% 
            (r %>% filter(clade %in%  c("N61","Columbea") ) %>% group_by(Chr,clade) %>% summarise(n=n()) %>% group_by(Chr)  %>% summarise(a=min(n)) %>% filter(a>MAS))$Chr
          & r$clade %in% c("N61","Columbea") ,c(
            "Gene"  , "p", "clade",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O","Chromosome", "ws"        
            ,"we"      ,   "size"     ,  "wstart" ,    "Chr" )],
        measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" )) %>% 
  filter(variable == "C1,C2|S,O") %>% 
  mutate( clade = fct_recode(clade, "S2023" = "N61", 
                             "J2014 (Columbea)" = "Columbea"
                             #"3rd alt" = "C2,S|C1,O"
  )) %>%
  ggplot(aes(x=p,y=ifelse(value==1,1,0),color=clade))+
  theme_classic()+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline",size=0.8)+
  scale_x_continuous(name="Chromosomes",breaks = c() )+# =co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = c(.15,.92),legend.direction = "vertical")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  geom_text(aes(y=-0.015,x=(start+end)/2,label=substr(Chr,0,2)),color="black",alpha=1,
            size=2.5,data=co[co$end-co$start>13000,])+
  geom_segment(aes(y=-0.029,yend=-0.029,x=start,xend=end),color="black",size=0.75,data=co,alpha=1,
               arrow=arrow(ends = "both",type = "closed",length =unit(1.7,"pt")))+
  geom_vline(aes(xintercept=start),color="gray60",size=0.2,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray60",size=0.2,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = MAS,linetype=1,alpha=0.6)+
  scale_size_manual(name="",values=c(1,1,0.5)*0.8)+
  scale_color_manual(name="",values=c("#1b9e77","#d95f02","#e7298a"))+
  scale_alpha_manual(name="",values=c(0.6,0.6,0.4))+
  scale_y_continuous(name="Monophyletic portion (n = 200)")+
  coord_cartesian(ylim=c(0.0,1))

p8

p9 = melt(r[r$Chr ==4    & r$clade %in% c("N61","Columbea") ,c(
  "Gene"  , "p", "clade",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O","Chromosome", "ws"        
  ,"we"      ,   "size"     ,  "wstart" ,    "Chr" )],
  measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" )) %>% 
  filter(variable == "C1,C2|S,O") %>% 
  mutate( clade = fct_recode(clade, "S2023" = "N61", 
                             "J2014 (Columbea)" = "Columbea"
                             #"3rd alt" = "C2,S|C1,O"
  ) ) %>%
  ggplot(aes(x=we,y=ifelse(value==1,1,0),color=clade, size=variable))+
  theme_classic()+
  geom_hline(color="grey20",yintercept = 1/3,linetype=3)+
  scale_x_continuous(name="Positions along chromosome 4")+# =co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "none")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 50,linetype=1, alpha=0.7)+
  scale_size_manual(name="",values=c(1,1,0.5)*0.8)+
  scale_color_manual(name="",values=c("#1b9e77","#d95f02","#e7298a"))+
  scale_y_continuous(name="Monophyletic portion (n = 50)")+
  coord_cartesian(ylim=c(0.025,1))+
  geom_point(alpha=0.2,size=0.1)+
  geom_vline(xintercept = c(34470000,25030000,32670000,33510000,56810000,44130000),color="gray50",linetype=3)
#geom_vline(xintercept = c(488932, 497238 ,497757,498734,509454,522053),color="gray50",linetype=3)

p9



p10 = ggplot(aes(x=p,y=ifelse(value==1,1,0)),
       data= r3[r3$variable == "C1,C2|S,O",])+
  theme_classic()+facet_wrap(~clade,ncol=4)+
  #stat_summary(aes(x=min(p),yintercept=..y..),geom="hline")+
  scale_x_continuous(name="Chromosomes",breaks=co[cl,"end"],labels=co[cl,"Chr"])+
  theme(legend.position = "bottom")+
  #theme(legend.position = c(.2,.8),legend.direction = "vertical")+
  #geom_segment(aes(y=0.040,yend=0.040,x=start,xend=end),color="black",size=0.75,data=co,
  #             arrow=arrow(ends = "both",type = "closed",length =unit(1.5,"pt")))+
  #geom_text(aes(y=0.035,x=(start+end)/2,label=substr(Chr,0,2)),color="black",linetype=1,
  #          size=2.3,data=co[co$end-co$start>1500,])+
  geom_vline(aes(xintercept=start),color="gray80",size=0.1,data=co, linetype=3)+
  geom_vline(aes(xintercept=end),color="gray80",size=0.1,data=co, linetype=3)+
  geom_ma(aes(group=interaction(Chr,clade,variable)),n = 200,linetype=1)+
  scale_y_continuous(name="Monophyletic portion (n=200)",labels = percent)
#facet_wrap(~grepl("Col",paste(C1,C2)))
ggsave("all-clades-mono-QQS-facet.pdf",p10,width=12,height = 15) 
ggsave("all-clades-mono-QQS-facet.png",p10,width=12,height = 15) 


ggplot(aes(x=p,y=value,#color=interaction(C1,C2,sep="/"), 
           color=variable,size=variable,alpha=variable),
       data= r3[r3$clade =="N53",])+
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
  geom_ma(aes(group=interaction(Chr,variable)),n = 100,linetype=1)+
  scale_size_manual(name="topology",values=c(0.85,0.5,0.5)*0.8, labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_color_manual(name="topology",values=c("#107030","#C05050","#4080E0"), labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_alpha_manual(name="topology",values=c(0.6,0.33,0.33), labels=c("tinamous,rheas","tinamous,kiwis+emus","rheas,kiwis+emus"))+
  scale_y_continuous(name="Quadripartiton Quartet Support (QQS)",
                     labels=percent,breaks=c(1/3,1/2,1/4,0,3/4,1))
#ggsave("all-clades-QQS-tinamus-rhea.pdf",width=7,height = 5) 


########################### Combine

# Note: panelE is created in coalHMM folder with some manual annotations and added here. 
p5=ggplot() +
  theme_void() +
  coord_cartesian(xlim = c(0.045, 0.96), ylim = c(0.07, .92))+
  draw_image(image_read_pdf('panelE.pdf',density=500), scale = 1, x = 0, y = 0)
ggarrange(
  ggarrange(p6, p7, p8, p9,ncol = 2, nrow=2, labels = c("A", "B" , "C", "D"),label.y = c(1,1,1.05,1.05)),
  p5,
  ncol = 1, nrow=2, heights=c(0.75,0.25),labels=c(NA,"E"),label.y = 1.07
)

ggsave(paste("Figure2-n.pdf",sep=""),width=12*1.07,height = 13*1.07)

print("done!")


################# P-value calculations : rationale

v0 = ggplot(aes(x=p,y=stats::filter(`C1,C2|S,O`, filter = rep(1/20,20) , method="convolution")),data=r)+
  theme_classic()+
  geom_point(aes(color="N37"), data= r[r$Chr %in% c(1:30) & r$clade =="N37",c(1:3,21,11:13,14:20)],alpha=0.07,size=0.6)+
  geom_point(aes(color="Columbea"),data= r[r$Chr %in% c(1:30) & r$clade =="Columbea",c(1:3,21,11:13,14:20)],alpha=0.07,size=0.6)+
  scale_y_continuous("Mean QQS over 20 consecutive loci")+
  scale_x_continuous(name="Genomic positions")+
  scale_color_manual(name="",values=c("red","blue","black"))+ theme(legend.position = c(.8,.9))

v0
v1 = ggplot(aes(x=stats::filter(`C1,C2|S,O`-mean(`C1,C2|S,O`), filter = rep(1/20,20) , method="convolution")/
            sd(stats::filter(`C1,C2|S,O`,filter = rep(1/20,20) , method="convolution"),na.rm = T)),data=r
    )+
  theme_classic()+
  geom_density(aes(color="N37"), data= r[r$Chr %in% c(1:30) & r$clade =="N37",c(1:3,21,11:13,14:20)])+
  geom_density(aes(color="Columbea"),data= r[r$Chr %in% c(1:30) & r$clade =="Columbea",c(1:3,21,11:13,14:20)])+
  scale_x_continuous("Z-scores computed from QQS averaged over 20 consecutive loci")+
  stat_function(fun=dnorm,aes(color="Normal Distribution"),linetype=2,size=1.2)+
  scale_y_continuous(name="density")+
  scale_color_manual(name="",values=c("red","blue","black"))+ theme(legend.position = c(.8,.9))
v1

twowayp = function (x) pmin(1-pnorm(x) , pnorm(x))
v2 = ggplot(aes(x=p.adjust(twowayp(stats::filter(`C1,C2|S,O`-mean(`C1,C2|S,O`), filter = rep(1/20,20) , method="convolution")/
    sd(stats::filter(`C1,C2|S,O`,filter = rep(1/20,20) , method="convolution"),na.rm = T)),method="BH")),
    data=r)+
  theme_classic()+
  geom_histogram(aes(fill="N37"),data=r[r$Chr %in% c(1:15,17:24,26:28,"Z") & r$clade =="N37",c(1:3,21,11:13,14:20)],alpha=0.5)+
  geom_histogram(aes(fill="Columbea"),data=r[r$Chr %in% c(1:15,17:24,26:28,"Z") & r$clade =="Columbea",c(1:3,21,11:13,14:20)],alpha=0.5)+
  scale_x_log10("BH corrected p-values")+scale_y_log10()+geom_vline(xintercept = 0.01,color="black",linetype=2)+
  scale_fill_manual(name="",values=c("red","blue"))+ theme(legend.position = c(.2,.9))
v2

ggarrange(v0,
  ggarrange(v1,v2,ncol = 2, nrow=1, labels = c("B", "C")),
  ncol = 1, nrow=2, labels = c("A",NA))

ggsave(paste("p-value-defs.pdf",sep=""),width=10,height = 11)
#ggsave(paste("p-value-defs.png",sep=""),width=10,height = 11,dpi = 400)


######### The P-value function ####################

pvalues = function(x, chr,window = 21 ) {
  # Make sure the moving average is broken across boundaries of chromosomes. 
  bounds = which(chr[-1] != chr[-length(chr)])
  bounds = 0:(length(bounds)-1)+bounds
  for (i in bounds) x <- c(head(x, i-1), NA, tail(x, -(i-1)))
  b2=which(!is.na(x))
  # We are testing for either extreme
  twowayp = function (x) pmin(1-pnorm(x) , pnorm(x))
  # BH test correction + two-tailed test, centering around mean. 
  # Sorry for the hard to parse expression. The heavy lifting is done by filter
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


## Calculate p-values for the Manhattan plot
pr = r[r$Chr %in% c(1:30),] %>%  group_by(clade) %>%
  mutate(pvalue = pvalues(`C1,C2|S,O`,Chr,20))
tail(pr)

cl= co$Chr %in% c(1:6,8,10, 12,15,20,25)

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
  theme(legend.position = c(.19,.83),axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=13))+
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
ggsave("manhattan.png",width=9,height = 5.2,dpi = 600)




# Calculate Monophly across chromosomes
r %>% 
  filter(clade == "Columbea" ) %>%
  mutate(outlier = 
           Chr==4 & ((25030000<= ws & ws <= 32670000) | (33510000<= ws & ws <=34470000) |(44130000<= ws & ws <=56810000)) )%>%
  mutate(m= ifelse(`C1,C2|S,O`==1,1,0)) %>%
  select(Chr,m,outlier) %>%
  group_by(Chr,outlier) %>%
  summarise(mc = sum(m),t=n()) %>%
  filter(t > 50) %>%
  ggplot(aes(x=reorder(Chr,-t),y=mc,fill=cut(t,c(50,200,1000,2000,5000,15000)),color=!outlier))+
  theme_classic()+geom_bar(stat="identity",size=0.4)+
  #scale_fill_continuous(trans="log10",name="number of loci")+
  scale_fill_viridis_d(name="Total number of loci")+
  scale_y_continuous(name="Loci recovering columbea as monophyletic")+
  scale_x_discrete(name="Chromosome (ordered by the number of loci)")+
  theme(legend.position = c(.7,0.6))+
  scale_color_manual(values=c("red","grey30"),labels=c("Outlier region","Other regions" ),name="")
ggsave("monophyly-disperssion-50.pdf",width=6,height = 4)



#### Number of outlier loci
nrow(r[ r$clade == "Columbea" ,] %>%
       mutate(outlier = 
                Chr==4 & ((25030000<= ws & ws <= 32670000) | (33510000<= ws & ws <=34470000) |(44130000<= ws & ws <=56810000)) ) %>%
       select(Gene,clade,Chromosome, ws, we, wstart, outlier) %>% filter(outlier == T))

### List monophyletic clades
write.csv(x=
            r[r$`C1,C2|S,O` == 1 & r$clade == "Columbea" ,] %>%
            mutate(outlier = 
                     Chr==4 & ((25030000<= ws & ws <= 32670000) | (33510000<= ws & ws <=34470000) |(44130000<= ws & ws <=56810000)) ) %>%
            select(Gene,clade,Chromosome, ws, we, wstart, outlier), file="../root-to-tip/outliers.txt",row.names = FALSE)

########################### Some miscellaneous testing; obsolete; ignore. 


maproot = function(y) {
  rootfail = c(21952, 32668, 42529, 43658, 45524, 54762, 60281, 61138)
  vapply(y,function(x) x+sum(x>rootfail),FUN.VALUE = c(1))
}


summary(mt[!is.na(mt$size2) & !is.na(mt$mu) & mt$size2>0 ,"outlier"])

mt %>% filter(is.na(size2) & !is.na(size))
  
  setdiff(mtm[which(!is.na(mtm$X2)),]$X1,monos)
setdiff((mtm[which(!is.na(mtm$X2)),]$X1),mt$X1)
mtm[12477,]
head(monos)
head(monolen$V1)

g %>% filter(Gene ==444 & Taxa == "Columbea (J2014)")
r %>% filter(Gene ==8503 & clade == "Columbea" )
m %>% filter(Gene ==758 & Taxa == "Columbea (J2014)" )

setdiff(r[r$`C1,C2|S,O` == 1 & r$clade == "Columbea" ,"Gene"],
        m[!is.na(m$mono) & m$Taxa == "Columbea (J2014)","Gene"])

mt[ mt$Gene %in%
setdiff(maproot((mt %>% filter(!is.na(size)))$Gene),m[!is.na(m$mono) & m$Taxa == "Columbea (J2014)","Gene"])
,]

mt[ mt$Gene %in%
      setdiff(maproot((mt %>% filter(!is.na(size)))$Gene),
              r[r$`C1,C2|S,O` == 1 & r$clade == "Columbea" ,"Gene"])
    ,]

mt[ maproot(mt$Gene) %in%
      setdiff(
              r[r$`C1,C2|S,O` == 1 & r$clade == "Columbea" ,"Gene"],
              maproot((mt %>% filter(!is.na(size)))$Gene))
    ,]
  