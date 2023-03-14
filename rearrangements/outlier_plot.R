library(tidyverse)
library(ggthemes)
library(glue)

totab_4 <- read_csv('totab_chr4.csv.xz') %>% mutate(tag = 'chr4')
totab_4_n <- totab_4 %>% 
  separate(Description, c('ID', 'Species', 'Chr'), 
           sep = '\\.', 
           extra = 'merge', 
           remove = F)
totab_4_n$Species= factor(totab_4_n$Species)

levels(totab_4_n$Species) = list("Chicken"= "5_GalGal6",
                                 "Flamingo"="1_bPhoRub2",
                                 "Dove\n(Columbimorphae)" = "2_bStrTur1",
                                 "Sandgrouse\n(Columbimorphae)" = "1_bPteGut1",
                                 "Turaco\n(Otidimorphae)" = "1_bTauEry1",  
                                 "Cuckoo\n(Otidimorphae)" = "1_bCucCan1" , 
                                 "Stork\n(other)" = "1_bCicMag1" ,
                                 "Finch\n(other)" = "2_bTaeGut1")

totab_4_n %>% 
  group_by(sp, Block) %>% 
  mutate(
    Gal_start = Start[Species == 'Chicken'],
    Gal_end = End[Species == 'Chicken']
  ) %>% 
  ungroup() %>% 
  #filter(Species != 'Chicken') %>% 
  mutate(
    outlier =
      (Gal_start>25555144 & Gal_end< 33202185) |
      (Gal_start>34230207 & Gal_end< 34999560) |
      (Gal_start>44689897 & Gal_end< 57179311)
  ) %>%
  ggplot() +
  geom_segment(aes(x = Start, xend = End, 
                   y = Chr, yend = Chr,
                   color = outlier),
               size=3) + 
  theme_classic() + ylab("") + xlab("Position") + 
  scale_color_brewer(palette = "Set2",name="Outlier?")+
  theme(legend.position = c(-.08,.95))+
  facet_wrap(~Species, scales = 'free',nrow=2)

ggsave("outlier-regions.png",width=14,height = 7)
ggsave("outlier-regions.pdf",width=14,height = 7)

