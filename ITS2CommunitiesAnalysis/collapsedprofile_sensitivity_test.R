# SENSITIVITY TEST FOR ITS2 PROFILE DATA #

# clear working environment
rm(list=ls())

library(tidyverse)
library(vegan)


setwd("/Users/lumosmaximma/Desktop/coral/its2wd/profiles")
profs1 = read.csv("profs.csv")
profs = profs1[,-15]

# sensitivity test for effect of response in all recovered vs resistant: not sensitive
its2profNormRecov = subset(profs, response == "recovered")
its2profNormResist = subset(profs, response == "resistant")
its2profNormResp = rbind(its2profNormRecov, its2profNormResist)

nreps = 100
nbal = min(nrow(its2profNormRecov), nrow(its2profNormResist))

adon.F.list = list()
adon.R2.list = list()
adon.P.list = list()

for(n in 1:nreps){
  subbal = its2profNormResp %>%
    group_by(response) %>%
    sample_n(size = nbal)
  
  set.seed(123)
  bray1 = vegdist(subbal[,8:ncol(subbal)], method = "bray")
  test = adonis(bray1~response, data = subbal)
  set.seed(NULL)
  
  adon.F.list[[n]] = test$aov.tab$F.Model[1]
  adon.R2.list[[n]] = test$aov.tab$R2[1]
  adon.P.list[[n]] = test$aov.tab$`Pr(>F)`[1]
  print(n)
}

adon.F = do.call(rbind,adon.F.list)
adon.R2 = do.call(rbind,adon.R2.list)
adon.P = do.call(rbind,adon.P.list)
adon.stats = do.call(cbind, list(adon.F, adon.R2, adon.P))
adon.stats

# sensitivity test for reef zone of all colonies: not sensitive
nreps = 100
subzone = profs %>% count(reef_zone) %>% select(c(n))
nbal = min(subzone) 

adon.F.list = list()
adon.R2.list = list()
adon.P.list = list()

for(n in 1:nreps){
  subbal = profs %>%
    group_by(reef_zone) %>%
    sample_n(size = nbal)
  
  set.seed(123)
  bray1 = vegdist(subbal[,8:ncol(subbal)], method = "bray")
  test = adonis(bray1~reef_zone, data = subbal)
  set.seed(NULL)
  
  adon.F.list[[n]] = test$aov.tab$F.Model[1]
  adon.R2.list[[n]] = test$aov.tab$R2[1]
  adon.P.list[[n]] = test$aov.tab$`Pr(>F)`[1]
  print(n)
}

adon.F = do.call(rbind,adon.F.list)
adon.R2 = do.call(rbind,adon.R2.list)
adon.P = do.call(rbind,adon.P.list)
adon.stats = do.call(cbind, list(adon.F, adon.R2, adon.P))
adon.stats

# did not do sensitivity test for differences in resistant colonies between May and Oct bc betadisper was not significant

# did not do sensitivity test for differences in recovered colonies between Aug and Oct bc betadisper was not significant

# sensitivity test for differences based on reef zone in October samples: not sensitive
after_blALL = subset(profs, month == "October")

nreps = 100
balreszoneall_a = after_blALL %>% count(reef_zone) %>% select(c(n))
nbal = min(balreszoneall_a)

adon.F.list = list()
adon.R2.list = list()
adon.P.list = list()

for(n in 1:nreps){
  subbal = after_blALL %>%
    group_by(reef_zone) %>%
    sample_n(size = nbal)
  
  set.seed(123)
  bray1 = vegdist(subbal[,8:ncol(subbal)], method = "bray")
  test = adonis(bray1~reef_zone, data = subbal)
  set.seed(NULL)
  
  adon.F.list[[n]] = test$aov.tab$F.Model[1]
  adon.R2.list[[n]] = test$aov.tab$R2[1]
  adon.P.list[[n]] = test$aov.tab$`Pr(>F)`[1]
  print(n)
}

adon.F = do.call(rbind,adon.F.list)
adon.R2 = do.call(rbind,adon.R2.list)
adon.P = do.call(rbind,adon.P.list)
adon.stats = do.call(cbind, list(adon.F, adon.R2, adon.P))
adon.stats

# sensitivity test for differences based on bleaching status in May crest, all colonies: not sensitive
may_crest_prof = subset(profs, month == "May" & reef_zone == "crest")

nreps = 100
balmaycrestall = may_crest_prof %>% count(bleaching_status) %>% select(c(n))
nbal = min(balmaycrestall)

adon.F.list = list()
adon.R2.list = list()
adon.P.list = list()

for(n in 1:nreps){
  subbal = may_crest_prof %>%
    group_by(bleaching_status) %>%
    sample_n(size = nbal)
  
  set.seed(123)
  bray1 = vegdist(subbal[,8:ncol(subbal)], method = "bray")
  test = adonis(bray1~bleaching_status, data = subbal)
  set.seed(NULL)
  
  adon.F.list[[n]] = test$aov.tab$F.Model[1]
  adon.R2.list[[n]] = test$aov.tab$R2[1]
  adon.P.list[[n]] = test$aov.tab$`Pr(>F)`[1]
  print(n)
}

adon.F = do.call(rbind,adon.F.list)
adon.R2 = do.call(rbind,adon.R2.list)
adon.P = do.call(rbind,adon.P.list)
adon.stats = do.call(cbind, list(adon.F, adon.R2, adon.P))
adon.stats

# did not do sensitivity test for differences in crest colonies over time bc betadisper was not significant

# did not do sensitivity test for differences in forereef colonies over time bc betadisper was not significant  

# did not do sensitivity test for differences in backreef colonies over time bc betadisper was not significant 

# did not do sensitivity test for differences in community based on bleaching status in backreef bc betadisper was not significant 

# did not do sensitivity test for differences in bleaching status in all may colonies bc betadisper was not significant

# sensitivity test for differences in reef zone in may samples: not sensitive
may_prof = subset(profs, month == "May")
nreps = 100
balmaym = may_prof %>% count(reef_zone) %>% select(c(n))
nbal = min(balmaym)

adon.F.list = list()
adon.R2.list = list()
adon.P.list = list()

for(n in 1:nreps){
  subbal = may_prof %>%
    group_by(reef_zone) %>%
    sample_n(size = nbal)
  
  set.seed(123)
  bray1 = vegdist(subbal[,8:ncol(subbal)], method = "bray")
  test = adonis(bray1~reef_zone, data = subbal)
  set.seed(NULL)
  
  adon.F.list[[n]] = test$aov.tab$F.Model[1]
  adon.R2.list[[n]] = test$aov.tab$R2[1]
  adon.P.list[[n]] = test$aov.tab$`Pr(>F)`[1]
  print(n)
}

adon.F = do.call(rbind,adon.F.list)
adon.R2 = do.call(rbind,adon.R2.list)
adon.P = do.call(rbind,adon.P.list)
adon.stats = do.call(cbind, list(adon.F, adon.R2, adon.P))
adon.stats

