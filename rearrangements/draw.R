library(tidyverse)
library(ggthemes)
library(glue)

totab_4 <- read_csv('totab_CM030196.1_CICMAG.csv') %>% mutate(tag = 'chr4')
totab_5 <- read_csv('totab_CM030197.1_CICMAG.csv') %>% mutate(tag = 'chr5')
totab <- bind_rows(totab_4, totab_5)

totab2 <- totab %>% 
  # Separate description column into species and chromosome
  separate(Description, c('ID', 'Species', 'Chr'), sep = '\\.', extra = 'merge', remove = F) %>% 
  # Add sp tag with target species name
  group_by(tag, sp, Block) %>% 
  mutate(sp = Species[Species!='1_bCicMag1']) %>% 
  ungroup() %>% 
  select(Start, End, Species, Block, Chr, Length, Description, sp, Size, tag) %>% 
  pivot_longer(-c(Block, Species, Chr, Length, Description, sp, Size, tag)) %>% 
  # Arrange to make sense for geom_polygon
  arrange(sp, Block, Species, ifelse(Species!='1_bCicMag1', name, rev(name)))  %>% 
  group_by(tag, sp, Block) %>% 
  # Add order and target chr information information
  mutate(colo = unique(Chr[Species!='1_bCicMag1'])) %>% 
  group_by(tag, sp, colo) %>% 
  # Calculate total length per target chr
  mutate(tot = sum(Length)/2) %>% 
  group_by(tag, sp) %>% 
  # Rank target chr by total length
  mutate(rank = dense_rank(desc(tot))) %>% 
  # Get only two largest chr
  filter(rank %in% c(1, 2)) %>% 
  # Order chromosomes, with ref chr in middle
  mutate(
    idx = as.numeric(factor(Description)),
    idx = (idx-min(idx[Species != '1_bCicMag1']))/(max(idx[Species != '1_bCicMag1'])-min(idx[Species != '1_bCicMag1'])),
    idx = ifelse(Species == '1_bCicMag1', 1.5, rank)
    ) %>% 
  group_by(tag, sp, Description) %>% 
  mutate(max_chr = max(value)) %>% 
  ungroup() %>% 
  # Flip some chromosomes to ease interpretation
  mutate(
    value = ifelse(
      (Description %in% c("GCA_901699155.2_bStrTur1.1.LR594554.2")) & (tag == 'chr4'),
      max_chr-value,
      value),
    value = ifelse(
      !(Description %in% c("GCA_009819775.1_bPhoRub2.pri.chr5")) & (tag == 'chr5'),
      max_chr-value,
      value),
    sp = factor(sp, c("5_GalGal6", "1_bPhoRub2", "2_bStrTur1", "1_bPteGut1", 
                      "1_bTauEry1", "1_bCucCan1", "2_bTaeGut1"))
  ) %>%
  {.}

rainbowplot <- totab2 %>% 
  # Two of the chromosomes of 1_bPteGut1 were fused to ease interpretation
  mutate(
    value = ifelse(
      Description == "GCA_009769525.1_bPteGut1.pri.chr15",
      value-2605+55653895,
      value
    ),
    value = ifelse(
      Description == "GCA_009769525.1_bPteGut1.pri.chr20",
      value-0+47064605,
      value
    ),
    rank = ifelse(
      sp == "1_bPteGut1",
      1,
      rank
    ),
    Block = ifelse(
      sp == "1_bPteGut1",
      -Block,
      Block
    )
  ) %>%
  # Filter the largest of the target chromosomes
  filter(rank == 1) %>% 
  select(sp, Block, Species, Length, name, value, tag) %>% 
  # Change name of target name to a generic name
  mutate(Species = ifelse(Species == '1_bCicMag1', '1_bCicMag1', 'target')) %>% 
  pivot_wider(names_from = c(Species, name), values_from = c(value, Length)) %>% 
  group_by(tag, sp) %>% 
  # Flip blocks that are reversed based on the reference sequence
  mutate(
    value_target_End_2 = value_target_End-min(c(value_target_Start, value_target_End)),
    value_target_Start_2 = value_target_Start-min(c(value_target_Start, value_target_End)),
    value_target_End = value_target_End_2,
    value_target_Start = value_target_Start_2
  ) %>%
  group_by(tag, sp, Block) %>% 
  # Calculate midpoint of blocks for dotplot
  mutate(
    # col = (value_1_bCicMag1_End+value_1_bCicMag1_Start)/2,
    # x = (value_target_End+value_target_Start)/2,
    col = (value_target_End+value_target_Start)/2,
    x = (value_1_bCicMag1_End+value_1_bCicMag1_Start)/2,
    size = abs(value_target_End-value_target_Start)
    ) %>% 
  group_by(tag, sp) %>% 
  arrange(x) %>% 
  # mutate(sp = str_sub(sp, 8, 13)) %>% 
  ungroup() %>% 
  select(x, col, sp, size, tag) %>% 
  # bind_rows(
  #   tibble(
  #     x = seq(0, 85409045, as.numeric(res)),
  #     col = x,
  #     Order = factor('Phoenicopteriformes', order_order),
  #     sp = '1_bCicMag1',
  #     size = as.numeric(res)
  #   )
  # ) %>% 
  {.}

levels(rainbowplot$sp) = list("Chicken"= "5_GalGal6",
                              "Flamingo\n(Mirandornithes)" = "1_bPhoRub2",
     "Dove\n(Columbimorphae)" = "2_bStrTur1",
     "Sandgrouse\n(Columbimorphae)" = "1_bPteGut1",
      "Cuckoo\n(Otidimorphae)" = "1_bCucCan1" , 
     "Turaco\n(Otidimorphae)" = "1_bTauEry1",  
     #"Stork\n(other)" = "1_bCicMag1" ,
     "Finch\n(other)" = "2_bTaeGut1")
rainbowplot = rainbowplot[!grepl("Finch",rainbowplot$sp) & 
                            !grepl("Chicke",rainbowplot$sp) &
                            rainbowplot$tag!="chr5",]
cuts <- rainbowplot %>% 
  # Arrange table by target coordinates per species
  arrange(tag, sp, col) %>% 
  # Add the target index
  mutate(idx_target = 1:n()) %>% 
  # Arrange table by reference coordinates per species
  arrange(tag, sp, x) %>% 
  # For each species
  group_by(tag, sp) %>% 
  mutate(
    # Lag the target index by one
    lag_idx_target = lag(idx_target),
    # If the target index is different from the lagged index ± 1
    cond = (lag_idx_target != idx_target+1) | (lag_idx_target != idx_target-1)
  ) %>% 
  filter((lag_idx_target != idx_target+1) & (lag_idx_target != idx_target-1)) %>% 
  ungroup()

rainbowplot %>% #filter(grepl("Cuck",sp))  %>% 
  ggplot() +
  #geom_vline(aes(xintercept = x/1000000),
  #           color = 'grey', data = cuts, size = 0.1) +
  geom_point(aes(x/1000000, col/1000000, color = sp), size = 0.1) +
  #geom_point(aes(x/1000000, col/1000000), size = 0.1, color = 'black',
  #           data = filter(rainbowplot, between(x, 5.5e7, 9e7))) +
  facet_wrap(.~sp,nrow=2) +
  xlab('Position in Stork coordinates (Mb)') +
  ylab('Position in target coordinates (Mb)') +
  theme_few() +
  theme(legend.position = 'none') +
  scale_color_manual(values=c("#E6000E","#E6000E","#E6000E","#118EFF","#118EFF","#118EFF","#118EFF"))+
  #Dove
  #geom_vline(xintercept = c(51587158,53698656,55785465,60614553,60956969.5,64301766.5,74856316,82157457.5,82414653,82500793,83325162,83839255)/1000000)+
  #Flamingo
  #geom_vline(xintercept = c(57081129,57115396,60625380,64242837,82256821,82387464)/1000000)+
  #Turaco
  #geom_vline(xintercept = c(377990.5,739885.5,55909486,56046799,60732528,60756246,60949927.5,64493073,82527764.5,82705752,83628778,83663684,83680957.5,83734358)/1000000)+
  #Sandgrouse
  #geom_vline(xintercept = c(51457258,53825629.5,55757140.5,56061738.5,57006115,57082300.5,60531580,60832023.5,64390875.5,64578860.5,64870822.5,82177487.5,82476572,83765263.5)/1000000)+
  #Cuckoo
  #geom_vline(xintercept = c(64689420,65536386.5)/1000000)+
  coord_fixed()+
  annotate("rect", xmin = 60750296/1000000, xmax = 84036155/1000000, ymin = 0/1000000, ymax =86610582/1000000, alpha = .1)

ggsave('synteny_point_stork_new.png', width=12.5, height=3.5)
ggsave('synteny_point_stork_new.pdf', width=6, height=5)

write.table(x=rainbowplot %>% arrange(sp,tag,x) %>%
              mutate(sp=sub('\n',' ',sp)) ,
            file="figure1b.tsv",sep ="\t")


