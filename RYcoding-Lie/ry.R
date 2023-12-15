require(reshape2)
require(tidyverse)


###############################
r = read.csv('rycoding-qqs.txt',sep=" ",h=F)
r = dcast( data=r,V4+V1+V2~V3)
names(r)[1:3] = c("Gene","Set","RY")
head(r)

# Compute QQS
r$T1 = r$`1` - r$`4`
r$T2 = r$`2` - r$`4`
r$T3 = r$`3` - r$`4`
r$`C1,C2|S,O` = (r$T2+r$T3-r$T1)/(r$T1+r$T2+r$T3)
r$`C1,S|C2,O` = (r$T1+r$T3-r$T2)/(r$T1+r$T2+r$T3)
r$`C2,S|C1,O` = (r$T1+r$T2-r$T3)/(r$T1+r$T2+r$T3)
head(r)

# Brace yoruself! We have wrong indices in the r dataframe for outlier, non-outlier, etc. 
# We need to fix this. 

# read mapping between indices of outliers region genes, 
# non-outlier genes that happen to be on chr4, and control genes outside. 
outlierinds = read.csv("outlier-indices.txt",head=F,sep=" ")
nonoutlierinds = read.csv("non-outlier-indices.txt",head=F,sep=" ")
controlinds = read.csv("control-indices.txt",head=F,sep=" ")
head(controlinds)
# Use these mappings to find "correct" gene names.
r = rbind(merge(r[r$Set == "control",],controlinds,by.x="Gene",by.y="V1"),
      merge(r[r$Set == "outlier_chr4",],outlierinds,by.x="Gene",by.y="V1"),
      merge(r[r$Set == "non-outlier_chr4",],nonoutlierinds,by.x="Gene",by.y="V1"))
names(r)[1] = "Index"
names(r)[14] = "Gene"
names(r)[2] = "OrigSet"
head(r)

# Remove case with 0 relevant quartets
r = r[!is.na(r$`C2,S|C1,O`),]

# Now that we have correct gene names, use outliers.txt to 
#      correctly label rows as either outlier or not
trueoutliers = read.csv('outliers.txt',h=F)$V1
r$Set = ifelse(r$Gene %in% trueoutliers, "outlier", "control")

# Subsample exactly 500 out of 682 outlier genes at random
r[r$Gene %in% sample(unique(r[r$Set == "outlier","Gene"]),500),]$Set = "outlier-selected"

r3 = melt(r[,c(
  "Gene"  , "Set", "RY",   "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" )],
  measure.vars = c( "C1,C2|S,O" , "C1,S|C2,O" , "C2,S|C1,O" ))


r3 %>% group_by(variable, Set, RY) %>% 
  filter(!is.na(value)) %>% 
  summarise(means=mean(value),n=n()) %>% 
  #filter(variable == "C1,C2|S,O")  %>%
  pivot_wider(names_from = RY, values_from = means)

r  %>% group_by(OrigSet, Set, RY) %>% summarize(n=n())


################ Just checking that incorrectly labelled stuff are correctly detected
incorrectoutliers = which(! (outlierinds$V2  %in% trueoutliers)) 
incorrectoutliers == IncorrectlylabelledAsOutlier
incorrectnonoutliers = sort(c(nonoutlierinds[nonoutlierinds$V2 %in% trueoutliers,"V1"] , 
                              which(!(1:635 %in% nonoutlierinds$V1))))
incorrectcontrolinds = sort(c(controlinds[controlinds$V2 %in% trueoutliers,"V1"] , 
                              which(!(1:635 %in% controlinds$V1))))

presumedoutliers = read.csv("outliers.chr4.txt")
table(cut(presumedoutliers$ws,c(0,25030000,32670000,33510000,34470000,44130000,56810000, 156810000)))
names(presumedoutliers) = c("ws")
IncorrectlylabelledAsOutlier = with(presumedoutliers, which(! (ws >= 25030000 & ws <= 32670000 | 
                                                                 ws >= 33510000 & ws <= 34470000 | 
                                                                 ws >= 44130000 & ws <= 56810000)))


