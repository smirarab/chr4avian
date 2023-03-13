require(ggplot2); require(scales); require(reshape2)

l = read.csv('chr4-summary-ELB.tsv',sep="\t")
head(l)

nrow(l)


l$labelnew = ifelse( with(l,
              (mRNA_begin>25030000 & mRNA_end<32670000) |
              (mRNA_begin>33510000 & mRNA_end<34470000) |
              (mRNA_begin>44130000 & mRNA_end<56810000) ),"outlier","nonoutlier")
 
ggplot(aes(x=mRNA_begin,xend=mRNA_end,y=labelnew,yend=labelnew),data=l)+
  geom_segment()+theme_classic()


ggplot(aes(x=mRNA_begin),data=l)+
  stat_ecdf()+theme_classic()+
  annotate("rect", xmin = 25030000, xmax = 32670000, ymin = 0, ymax =1, alpha = .1)+
  annotate("rect", xmin = 33510000, xmax = 34470000, ymin = 0, ymax =1, alpha = .1)+
  annotate("rect", xmin = 44130000, xmax = 56810000, ymin = 0, ymax =1, alpha = .1)+
  geom_abline(color="blue",linetype=2,slope = 1/90161979)+
  xlab("Position along chromosome 4 of chicken (GalGal4)")+
  ylab("ECDF of gene start positions")
ggsave("genes-distribution.pdf",width=4.5,height = 4)

table(l$label)
table(l$labelnew)
summary(l$mRNA_end)

(32670000-25030000+34470000-33510000+56810000-44130000)/90211544
table(l$labelnew)[2]/(table(l$labelnew)[2]+table(l$labelnew)[1])

l$bestOmega= apply(cbind(l[,3+(1:3)*3],apply(l[,c("Tree1_two_ratio_lnL","Tree2_two_ratio_lnL","Tree3_two_ratio_lnL")],1,which.max)),1,function(x) x[x[4]])
l$pvalue= apply(cbind(l[,18+(1:3)*2],apply(l[,c("Tree1_two_ratio_lnL","Tree2_two_ratio_lnL","Tree3_two_ratio_lnL")],1,which.max)),1,function(x) x[x[4]])
table(l$pvalue<0.05,l$labelnew)
table(l$pvalue<0.01,l$labelnew)
table(l$pvalue<0.001,l$labelnew)

ggplot(aes(x=(mRNA_begin+mRNA_end)/2,
           color=as.factor(apply(l[,c("Tree1_two_ratio_lnL","Tree2_two_ratio_lnL","Tree3_two_ratio_lnL")],1,which.max)),
           y=bestOmega,
           shape=pvalue<0.01
           ),data=l)+
  geom_hline(yintercept = 1,color="grey30",linetype=2)+
  geom_point(size=1)+theme_classic()+
  annotate("rect", xmin = 25030000, xmax = 32670000, ymin = 0, ymax =1000, alpha = .1)+
  annotate("rect", xmin = 33510000, xmax = 34470000, ymin = 0, ymax =1000, alpha = .1)+
  annotate("rect", xmin = 44130000, xmax = 56810000, ymin = 0, ymax =1000, alpha = .1)+
  #stat_summary(aes(x=25030000/2),color="black",data = l[l$mRNA_begin< 25030000,])+
  #stat_summary(aes(x=50300000/2),color="black")+
  xlab("Position along chromosome 4 of chicken (GalGal4)")+
  scale_shape_manual(name=expression(p<0.01),values=c(21,19))+
  scale_y_log10(name=expression(omega~"of the focal branch in the highest likelihood gene tree"))+
  scale_color_manual(values=c("#118EFF","#E6000E","#900890"),
                     name="Best Tree",labels=c("T1","T2","T3"))+
  theme(legend.position = c(.82,.14),legend.direction = "horizontal",
        legend.box.margin = margin(0),
        legend.margin = margin(0))
ggsave("omega.pdf",width=7.4,height = 5.2)
nrow(l)


