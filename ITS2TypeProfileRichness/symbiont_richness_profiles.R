rm(list=ls())

library(tidyverse)
library("reshape2")
library(vegan)
library(pairwiseAdonis)
library(plyr)
library(lme4) 
library(nlme)
library(MuMIn)

setwd("/Users/lumosmaximma/Desktop/coral/its2wd/profiles")

profs = read.csv("profs.csv")
profs = profs[,1:14]
profs$response = as.factor(profs$response)
profs$depth_m = as.factor(profs$depth_m)
profs$month = as.factor(profs$month)
profs$reef_zone = as.factor(profs$reef_zone)
profs$bleaching_status = as.factor(profs$bleaching_status)

# calculating profile richness for each sample (number of symbiont type profiles present)
richness_prof = ddply(profs, ~Sample, function(x) {
  data.frame(Richness = sum(x[8:ncol(profs)]>0))  
}) %>% select(c(Richness))
richprofs = cbind(profs, richness_prof)
richprofs$Individual = substr(richprofs$Sample,1,3)

#richness differences between all recovered and resistant: not significant 0.1962
its2profNormRecov = subset(richprofs, response == "recovered")
its2profNormResist = subset(richprofs, response == "resistant")
its2profNormResp = rbind(its2profNormRecov, its2profNormResist)
results1 = glmer(Richness~response+depth_m+(1|Individual), data = its2profNormResp, family = poisson)
summary(results1)
exp(confint(results1))
#richness differences between reef zones for all colonies: crest and forereef not significantly different frome each other, other comparisons are though
results2 = glmer(Richness~reef_zone+(1|Individual), data = richprofs, family = poisson)
summary(results2)
exp(confint(results2))
#NOTE: had to remove depth_m from the above model because of rank deficiency (also reef_zone and depth are singularities... does not change the fact that reef_zone p values were all highly significant
results3 = glmer(Richness~relevel(reef_zone,ref="forereef")+(1|Individual), data = richprofs, family = poisson)
summary(results3) 
exp(confint(results3))
#richness differences between paired May and October in resistant colonies: not significant 0.853
richprofsresist_pair = richprofs %>%
  group_by(Individual) %>%
  filter(response == "resistant" & n()>1) #selecting resistant paired colonies
results6 = glm(Richness~month, data = richprofsresist_pair, family = poisson)
summary(results6)
exp(confint(results6))
#richness differences between paired August and October in recovered colonies: not significant 0.292
richprofsrecov_pair = richprofs %>%
  group_by(Individual) %>%
  filter(response == "recovered" & n()>1) #selecting recovered paired colonies
results7 = glm(Richness~month, data = richprofsrecov_pair, family = poisson)
summary(results7)
exp(confint(results7))
#richness differences between reef zones for all October samples: forereef and crest not different but other comparisons are
after_blALL = subset(richprofs, month == "October")
results12.1 = glm(Richness~reef_zone, data = after_blALL, family = poisson)
summary(results12.1)
results12.2 = glm(Richness~relevel(reef_zone,ref = "forereef"), data = after_blALL, family = poisson)
summary(results12.2)
#NOTE: did not include depth_m even though collinear with reef zone because it did not change p value or estimates of effect
#richness differences between bleaching statuses in May crest samples: significant 0.000845
may_crest = subset(richprofs, month == "May" & reef_zone == "crest")
results14 = glm(Richness~bleaching_status, data = may_crest, family = poisson)
summary(results14)
exp(confint(results14))
#richness differences for crest (actually shallow forereef) based on month (May vs Oct): not significant 0.188
crest = subset(richprofs, reef_zone == "crest")
results15 = glm(Richness~month, data = crest, family = poisson)
summary(results15)
exp(confint(results15))
#richness differences for forereef (actually deep forereef) based on month (Aug vs Oct): not significant 0.3313
forereef = subset(richprofs, reef_zone == "forereef")
results16 = glm(Richness~month, data = forereef, family = poisson)
summary(results16)
exp(confint(results16))
#richness differences for backreef based on month (May vs Oct): not significant 0.439
backreef = subset(richprofs, reef_zone == "backreef")
results17 = glm(Richness~month, data = backreef, family = poisson)
summary(results17)
exp(confint(results17))
#richness differences between bleaching statuses in May backreef samples (unknown strategies): not significant 0.248
may_br = subset(richprofs, month == "May" & reef_zone == "backreef")
results19 = glm(Richness~bleaching_status, data = may_br, family = poisson)
summary(results19)
exp(confint(results19))
#richness differences between bleaching statuses in May all samples: health significant, reef zone not, but marginally
may_health = subset(richprofs, month == "May")
results20 = glm(Richness~bleaching_status+reef_zone, data = may_health, family = poisson)
summary(results20)
exp(confint(results20))
#chose to include bleaching status and reef zone because there are significant differences in reef zone so wanted to take that into account
#richness differences between reef zones for all May samples: forereef and crest not different but other comparisons are
mayall = subset(richprofs, month == "May")
results24 = glm(Richness~reef_zone, data = mayall, family = poisson)
summary(results14)

#size analysis: is there an impact of colony size on richness?
sizenas = read.csv("richprofssize.csv")
size = na.omit(sizenas)
size$response = as.factor(size$response)
size$depth_m = as.factor(size$depth_m)
size$month = as.factor(size$month)
size$reef_zone = as.factor(size$reef_zone)
size$bleaching_status = as.factor(size$bleaching_status)
size$Richness = as.numeric(size$Richness)
size$colony_area = as.numeric(size$colony_area)
results21.1 = lme(Richness~colony_area+reef_zone, data = size, random = (~1|Individual))
summary(results21.1) #colony area does not impact richness (0.9110), but reef zone does
results21.2 = lme(Richness~colony_area+relevel(reef_zone, ref="crest"), data = size, random = (~1|Individual))
summary(results21.2)
results22.1 = lme(colony_area~reef_zone, data = size, random = (~1|Individual))
results22.2 = lme(colony_area~relevel(reef_zone, ref="crest"), data=size, random = (~1|Individual))
summary(results22.1) #significant difference in size between reef zones
summary(results22.2)
results23 = lm(Richness~colony_area, data = size)
summary(results23) 

#colony size and dominant symbiont type
size$MaxSym = colnames(size)[max.col(size[,10:16],ties.method = "first")+9]
size$MaxSym = as.factor(size$MaxSym)
sizesym = ggplot(data = size)+
  geom_boxplot(mapping = aes(x = MaxSym, y = colony_area))
sizesym

#RICHNESS PLOTS
richprofs$month = ordered(richprofs$month, levels = c("May", "August", "October"))
richprofs$Richness = as.factor(richprofs$Richness)
Richcountszone = data.frame(richness = factor(), zone = factor(), Month = factor(), counts = double())
for(Richness3 in levels(richprofs$Richness)) {
  for(reef_zone3 in levels(richprofs$reef_zone)) {
    for(month3 in levels(richprofs$month)) {
        countall1 = nrow(richprofs %>% filter(Richness == Richness3, reef_zone == reef_zone3, month == month3))
        if(countall1 != 0) {
          Richcountszone = Richcountszone %>% add_row(richness = Richness3, zone = reef_zone3, Month = month3, counts = countall1)
        } 
      }
    }
  }
normalize_month_counts = TRUE #TRUE or FALSE to either normalize by samples collected per month or not
scale_factor = 1
Richcountszone$Month = as.factor(Richcountszone$Month)
Richcountszone$zone = as.factor(Richcountszone$zone)
if(normalize_month_counts){
  month_n = richprofs %>% dplyr::group_by(month) %>% dplyr::count(reef_zone)
  scale_factor = 23
  for(x in levels(Richcountszone$Month)){
    for(y in levels(Richcountszone$zone)){
      Richcountszone$counts[Richcountszone$Month == x & Richcountszone$zone == y] = Richcountszone$counts[Richcountszone$Month == x & Richcountszone$zone == y]/month_n$n[month_n$month == x & month_n$reef_zone == y]
    }
  }
  
}
zone.labs <- c("backreef", "shallow forereef", "deep forereef")
names(zone.labs) <- c("backreef", "crest", "forereef")

#circles normalized to number of samples collected each month
profzone = ggplot(data = Richcountszone)+
  geom_point(mapping = aes(x = Month, y = richness, color = Month, fill = Month), shape = 21, color = "black", size = scale_factor*(Richcountszone$counts))+
  facet_wrap(~ zone, scale = "free_x", labeller = labeller(zone = zone.labs))+
  labs(x = "Month", y = "Richness")+
  theme_classic(base_size=12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  scale_fill_manual(values = c("#ADD8E6","#71649D","#E59F2E"))+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("May" = "May", "August" = "Aug", "October" = "Oct"))
profzone
ggsave("profzone.pdf", width=5.25, height=4.1125, units="in")

# resistant and recovered colonies, separated by reef zone and then month and colored by response (shows the collinearity)
its2profNormResp$Richness = as.factor(its2profNormResp$Richness)
Richcountsresp = data.frame(richness = factor(), zone = factor(), Month = factor(), Response = factor(), counts = double())
for(Richness2 in levels(its2profNormResp$Richness)) {
  for(reef_zone2 in levels(its2profNormResp$reef_zone)[-1]) {
    for(month2 in levels(its2profNormResp$month)) {
      for(response2 in levels(its2profNormResp$response)[2:3]) {
        countall = nrow(its2profNormResp %>% filter(Richness == Richness2, reef_zone == reef_zone2, month == month2, response == response2))
        if(countall != 0) {
          Richcountsresp = Richcountsresp %>% add_row(richness = Richness2, zone = reef_zone2, Month = month2, Response = response2, counts = countall)
        } 
      }
    }
  }
}
normalize_month_counts2 = TRUE #TRUE or FALSE to either normalize by samples collected per month or not
scale_factor = 1
Richcountsresp$Month = as.factor(Richcountsresp$Month)
Richcountsresp$zone = as.factor(Richcountsresp$zone)
if(normalize_month_counts2){
  month_n2 = its2profNormResp %>% dplyr::group_by(month) %>% dplyr::count(reef_zone)
  scale_factor = 23
  for(x in levels(Richcountsresp$Month)){
    for(y in levels(Richcountsresp$zone)){
      Richcountsresp$counts[Richcountsresp$Month == x & Richcountsresp$zone == y] = Richcountsresp$counts[Richcountsresp$Month == x & Richcountsresp$zone == y]/month_n2$n[month_n2$month == x & month_n2$reef_zone == y]
    }
  }
  
}
zone.labs <- c("backreef", "shallow forereef", "deep forereef")
names(zone.labs) <- c("backreef", "crest", "forereef")
profresp = ggplot(data = Richcountsresp)+
  #geom_boxplot(data = richprofs, mapping = aes(x = Month, y = richness), width = 0.6)+
  facet_wrap(~ zone, scale = "free_x", labeller = labeller(zone = zone.labs))+
  labs(x = "Month", y = "Richness")+
  theme_classic(base_size=12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  geom_point(mapping = aes(x = Month, y = richness, color = Response, fill = Response), shape = 21, color = "black", size = scale_factor*Richcountsresp$counts)+
  scale_fill_manual(values = c("#960200","#FFD046"))+
  scale_x_discrete(labels=c("May" = "May", "August" = "Aug", "October" = "Oct"))
profresp
ggsave("profresp.pdf", width=6, height=4.7, units="in")

# healthy vs bleached in may, separated by reef zone, with facet wrap
may_health$Richness = as.factor(may_health$Richness)
Richcountshealth = data.frame(richness = factor(), zone = factor(), health = factor(), Response = factor(), counts = double())
for(Richness4 in levels(may_health$Richness)) {
  for(reef_zone4 in levels(may_health$reef_zone)) {
    for(health4 in levels(may_health$bleaching_status)) {
        countall2 = nrow(may_health %>% filter(Richness == Richness4, reef_zone == reef_zone4, bleaching_status == health4))
        if(countall2 != 0) {
          Richcountshealth = Richcountshealth %>% add_row(richness = Richness4, zone = reef_zone4, health = health4, counts = countall2)
        } 
      }
    }
}
normalize_month_counts3 = TRUE #TRUE or FALSE to either normalize by samples collected per month or not
scale_factor = 1
Richcountshealth$health = as.factor(Richcountshealth$health)
Richcountshealth$zone = as.factor(Richcountshealth$zone)
if(normalize_month_counts3){
  month_n3 = may_health %>% dplyr::group_by(bleaching_status) %>% dplyr::count(reef_zone)
  scale_factor = 23
  for(x in levels(Richcountshealth$health)){
    for(y in levels(Richcountshealth$zone)){
      Richcountshealth$counts[Richcountshealth$health == x & Richcountshealth$zone == y] = Richcountshealth$counts[Richcountshealth$health == x & Richcountshealth$zone == y]/month_n3$n[month_n3$bleaching_status == x & month_n3$reef_zone == y]
    }
  }
  
}
profhealth = ggplot(data = Richcountshealth)+
  facet_wrap(~zone, scale = "free_x")+
  labs(x = "Bleaching status", y = "Richness")+
  theme_classic(base_size=12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  geom_point(mapping = aes(x = health, y = richness, color = health, fill = health), shape = 21, color = "black", size = scale_factor*Richcountshealth$counts)+
  scale_fill_manual(values = c("#BFACC8","#2A4470"), name = "Bleaching status")+
  theme(legend.position = "none")
profhealth
ggsave("profhealth.pdf", width=5, height=4.7, units="in")

# richness over colony area
richsizeprof = ggplot(data = size)+
  geom_point(mapping = aes(x = colony_area, y = Richness, color = reef_zone), size = 4.5, alpha = 0.7)+
  labs(x = "Colony area", y = "Richness")+
  theme_classic(base_size=12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  scale_fill_manual(values = c("#BFACC8","#2A4470"), name = "Reef zone")+
  scale_color_manual(values = c("#71649D","#ADD8E6","#E59F2E"))
richsizeprof
ggsave("richsizeprof.pdf", width=5.5, height=4.5, units = "in")
