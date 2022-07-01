# ITS2 PROFILE DATA AND COMMUNITY STRUCTURE ANALYSES #

#clear working environment
rm(list=ls())

library(tidyverse)
library(reshape2)
library(MCMC.OTU)
library(vegan)
library(edgeR)
library(pairwiseAdonis)
library(lme4) 

setwd("/Users/lumosmaximma/Desktop/coral/its2wd/profiles")

#ITS2 sequences analysis
#upload output data from SymPortal: absolute profiles and metadata; cleaned data is the SymPortal output with everything except the clades removed; also removed samples that failed SymPortal 
its2profs = read.delim("20210730T152455.profiles.absolute.clean.txt", header = TRUE, check.names = FALSE) #deleted samples 213Mc, 219Mc, 276Mc, 321Af bc no profiles (too low reads)
head(its2profs)

#import additional sample metadata and merge with the SymPortal output above
its2profMetaData = read.csv("ahya_metadata_its2_2.csv")
its2prof = cbind(its2profs[1], its2profMetaData[,2:length(its2profMetaData)], its2profs[,2:length(its2profs)])
head(its2prof)

#transpose df so columns become rows and rows become columns, while removing metadata
its2profTransposed = t(its2prof[,7:length(its2prof[1,])])
its2profList = DGEList(counts = its2profTransposed)
head(its2profList$samples)

#weighted trimmed mean of M-values normalization in edgeR
its2profNorm = calcNormFactors(its2profList, method = "TMM")
head(its2profNorm$samples)
its2TMM2 = t(cpm(its2profNorm, normalized.lib.sizes = TRUE))
its2profNorm = cbind(its2prof[,c(1:6)], its2TMM2)
head(its2profNorm)

# #sanity check: looks good!
# apply(its2profPerc[, c(7:(ncol(its2profPerc)))], 1, function(x) {sum(x, na.rm = T)})

# collapse profiles
Aprof = its2profNorm[,7:10] %>%
  rowSums()%>%
  as.data.frame()
names(Aprof)[1] = 'A1'
Aprof$A1ee = its2profNorm[,11]

CprofC3C11 = its2profNorm[,18:21] %>%
  rowSums()%>%
  as.data.frame()
names(CprofC3C11)[1] = 'C3/C115/C116'
CprofC3ae = its2profNorm[,c(13, 15:17)] %>%
  rowSums()%>%
  as.data.frame()
names(CprofC3ae)[1] = 'C3ae'
CprofC1 = its2profNorm[,c(12, 14)] %>%
  rowSums()%>%
  as.data.frame()
names(CprofC1)[1] = 'C1'

Dprof = its2profNorm[,c(22:25, 27)] %>%
  rowSums()%>%
  as.data.frame()
names(Dprof)[1] = 'D1'
Dprof$D6 = its2profNorm[,26] 

its2profNorm = cbind(its2profNorm[,1:6], Aprof, CprofC1, CprofC3ae, CprofC3C11, Dprof)
its2profNorm$Individual = substr(its2profNorm$Sample,1,3)
#write.csv(its2profNorm, file = "profs.csv")

# construct ITS2 profiles barplot
its2profNorm$depth_m = as.factor(its2profNorm$depth_m)
meltedprof = melt(its2profNorm)
resistantprof = subset(meltedprof, response == "resistant")
recoveredprof = subset(meltedprof, response == "recovered")

backreefprof = subset(meltedprof, reef_zone == "backreef")
crestprof = subset(meltedprof, reef_zone == "crest")
forereefprof = subset(meltedprof, reef_zone == "forereef")

meltedmay = subset(meltedprof, month == "May")
mhealthyprof = subset(meltedmay, bleaching_status == "healthy")
mbleachedprof = subset(meltedmay, bleaching_status == "bleached")

mcrest = subset(meltedmay, reef_zone == "crest")
mbr = subset(meltedmay, reef_zone == "backreef")

resistant_pair = its2profNorm %>%
  group_by(Individual) %>%
  filter(response == "resistant" & n()>1) #selecting resistant paired colonies
meltresist = melt(resistant_pair)

recovered_pair = its2profNorm %>%
  group_by(Individual) %>%
  filter(response == "recovered" & n()>1) #selecting recovered paired colonies
meltrecov = melt(recovered_pair)

#picking color palette for stacked bar plots
#colorFuncpurp = colorRampPalette(c("#d9c4ff","#36175E"))
#plot(rep(1,6),col=colorFuncpurp(6),pch=19,cex=3)
#colorFuncyel = colorRampPalette(c("#ffd28d","#ff8000"))
#plot(rep(1,10),col=colorFuncyel(10),pch=19,cex=3)
#colorFuncgreen = colorRampPalette(c("#99E385", "#2c4b24"))  
#plot(rep(1,5),col=colorFuncgreen(5),pch=19,cex=3)

gradientgop = c("#99E385","#5FA44C",
                "#FFDF98","#FF8A10","#D95D00",
                "#D9C4FF","#775C9E")                

# stacked percent bar plot of resistant colony profiles by month
resistantbarpC = ggplot(resistantprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ month, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 type profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
resistantbarpC
ggsave("resistantbarpC.pdf", width = 4, height = 2, unit = "in")

# stacked percent bar plot of recovered colony profiles by month
recoveredprof$month = ordered(recoveredprof$month, levels = c("May", "August", "October"))
recoveredbarpC = ggplot(recoveredprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ month, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 type profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
recoveredbarpC
ggsave("recoveredbarpC.pdf", width = 4, height = 2, unit = "in")

#stacked bar plot of profiles from colonies in the backreef by month
backreefprof$month = ordered(backreefprof$month, levels = c("May", "August", "October"))
backreefbarpC = ggplot(backreefprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ month, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
backreefbarpC
ggsave("backreefbarpC.pdf", width = 4.05, height = 2, unit = "in")

#stacked bar plot of profiles from colonies in the crest (shallow forereef) by month
crestprof$month = ordered(crestprof$month, levels = c("May", "August", "October"))
crestbarpC = ggplot(crestprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ month, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
crestbarpC
ggsave("crestbarpC.pdf", width = 4.05, height = 2, unit = "in")

#stacked bar plot of profiles from colonies on the forereef (deep forereef) by month
forereefprof$month = ordered(forereefprof$month, levels = c("May", "August", "October"))
forereefbarpC = ggplot(forereefprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ month, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
forereefbarpC
ggsave("forereefbarpC.pdf", width = 4.05, height = 2, unit = "in")

#stacked bar plot of profiles from healthy colonies in may
mhealthybarpC = ggplot(mhealthyprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ reef_zone, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
mhealthybarpC
ggsave("mhealthybarpC.pdf", width = 4, height = 2, unit = "in")

#stacked bar plot of profiles from bleached colonies in may
mbleachedbarpC = ggplot(mbleachedprof, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ reef_zone, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
mbleachedbarpC
ggsave("mbleachedbarpC.pdf", width = 4, height = 2, unit = "in")

#stacked bar plot of profiles from crest may colonies
mcrestbarpC = ggplot(mcrest, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ bleaching_status, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
mcrestbarpC
ggsave("mcrestbarpC.pdf", width = 5.5, height = 2, unit = "in")

#stacked bar plot of profiles from bleached colonies in may
mbrbarpC = ggplot(mbr, aes(x = Sample, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ bleaching_status, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
mbrbarpC
ggsave("mbrbarpC.pdf", width = 3.5, height = 2, unit = "in")

#stacked bar plot of resistant colonies over time, grouped by individual
resistantpairbarpC = ggplot(meltresist, aes(x = month, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ Individual, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 type profile")+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
resistantpairbarpC
ggsave("resistantpairbarpC.pdf", width = 6, height = 2.5, unit = "in")

#stacked bar plot of paired recovered colonies grouped by individual
recoveredpairbarpC = ggplot(meltrecov, aes(x = month, y = value, fill = variable))+
  geom_bar(position = "fill",stat = "identity", color = "black", size = 0.25)+
  facet_grid(~ Individual, scales = "free", space = "free")+
  ylab("Proportion")+
  labs(fill = "ITS2 type profile")+
  scale_fill_manual(name = "ITS2 Type Profile", values = gradientgop)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none")
recoveredpairbarpC
ggsave("recoveredpairbarpC.pdf", width = 5.5, height = 2.5, unit = "in")

#PROFILE STATISTICS: PERMANOVAs
#NOTE: did not include random effect of individuals because it did not change any of the model output values
its2profNorm = its2profNorm %>% relocate(Individual, .before = A1)
resistant_pair = resistant_pair %>% relocate(Individual, .before = A1)
recovered_pair= recovered_pair %>% relocate(Individual, .before = A1)
its2profNormTAXA = its2profNorm[,-c(1:7)]
its2profNormENV = its2profNorm[,-c(8:length(its2profNorm))]

#differences in community based on response (resistant vs. recovered): significant 0.001
its2profNormRecov = subset(its2profNorm, response == "recovered")
its2profNormResist = subset(its2profNorm, response == "resistant")
its2profNormResp = rbind(its2profNormRecov, its2profNormResist)
distprof2 = vegdist(its2profNormResp[,8:ncol(its2profNormResp)], method = "bray")
dispprof2 = betadisper(distprof2,its2profNormResp$response)
anova(dispprof2) #significant 0.004794
set.seed(123)
adonis(its2profNormResp[,8:ncol(its2profNormResp)]~response+depth_m, data=its2profNormResp[,1:6]) 
#differences in community based on reef zone for ALL colonies: significant 0.001
distprof3 = vegdist(its2profNormTAXA, method = "bray")
dispprof3 = betadisper(distprof3,its2profNormENV$reef_zone)
anova(dispprof3) #significant 2.2e-16
set.seed(123)
adonis(its2profNormTAXA~reef_zone, data=its2profNormENV) 
pairwise.adonis(its2profNormTAXA, factors = its2profNormENV$reef_zone, sim.method = "bray", p.adjust.m = "BH", perm = 9999)
#May vs. October in resistant colonies: not significant 0.481
distprof5 = vegdist(resistant_pair[,8:ncol(resistant_pair)], method = "bray")
dispprof5 = betadisper(distprof5,resistant_pair$month)
anova(dispprof5) #not significant 0.214
set.seed(123)
adonis(resistant_pair[, c(8:length(resistant_pair))]~month, data=resistant_pair, permutations=999,method="bray")
#August vs. October in recovered colonies: not significant 0.706
recovered_pair = read.csv("norm_profs_AO.csv")
distprof6 = vegdist(recovered_pair[,8:ncol(recovered_pair)], method = "bray")
dispprof6 = betadisper(distprof6,recovered_pair$month)
anova(dispprof6) #not significant 0.7042
set.seed(123)
adonis(recovered_pair[, c(8:length(recovered_pair))]~month, data=recovered_pair, permutations=999,method="bray")
#reef zone in October samples: significant 0.001
after_blALL = subset(its2profNorm, month == "October")
distprof10 = vegdist(after_blALL[,8:length(after_blALL)], method = "bray")
dispprof10 = betadisper(distprof10, after_blALL$reef_zone)
anova(dispprof10) #significant 2.219e-14
set.seed(123)
adonis(after_blALL[, c(8:length(after_blALL))]~reef_zone, data=after_blALL, permutations=999, method="bray")
pairwise.adonis(after_blALL[,8:length(after_blALL)], factors = after_blALL$reef_zone, sim.method = "bray", p.adjust.m = "BH", perm = 9999)
#bleaching status in May crest samples: significant 0.001
may_crest_prof = subset(its2profNorm, month == "May" & reef_zone == "crest")
distprof12 = vegdist(may_crest_prof[8:length(may_crest_prof)], method = "bray")
dispprof12 = betadisper(distprof12, may_crest_prof$bleaching_status)
anova(dispprof12) #significant 1.636-5
set.seed(123)
adonis(may_crest_prof[, c(8:length(may_crest_prof))]~bleaching_status, data = may_crest_prof, permutations = 999, method="bray")  
#differences in community over time (May vs Oct) on the crest (actually shallow forereef): not significant 0.523
crest = subset(its2profNorm, reef_zone == "crest")
distprof14 = vegdist(crest[8:length(crest)], method = "bray")
dispprof14 = betadisper(distprof14, crest$month)
anova(dispprof14) #not significant 0.9603
set.seed(123)
adonis(crest[, c(8:length(crest))]~month, data = crest, permutations = 999, method="bray")  
#differences in community over time (Aug vs Oct) on the forereef (actually deeper forereef): not significant 0.228
forereef = subset(its2profNorm, reef_zone == "forereef")
distprof15 = vegdist(forereef[8:length(forereef)], method = "bray")
dispprof15 = betadisper(distprof15, forereef$month)
anova(dispprof15) #not significant 0.2077
set.seed(123)
adonis(forereef[, c(8:length(forereef))]~month, data = forereef, permutations = 999, method="bray")
#differences in community over time (May vs Oct) on the backreef: not significant 0.846
backreef = subset(its2profNorm, reef_zone == "backreef")
distprof16 = vegdist(backreef[8:length(backreef)], method = "bray")
dispprof16 = betadisper(distprof16, backreef$month)
anova(dispprof16) #not significant 0.1204
set.seed(123)
adonis(backreef[, c(8:length(backreef))]~month, data = backreef, permutations = 999, method="bray")   
#differences in community based on bleaching status in backreef: not significant 0.807
may_br_prof = subset(its2profNorm, month == "May" & reef_zone == "backreef")
distprof18 = vegdist(may_br_prof[8:length(may_br_prof)], method = "bray")
dispprof18 = betadisper(distprof18, may_br_prof$bleaching_status)
anova(dispprof18) #not significant 0.4801
set.seed(123)
adonis(may_br_prof[, c(8:length(may_br_prof))]~bleaching_status, data = may_br_prof, permutations = 999, method="bray")  
#differences in community based on bleaching status in may: health not significant 0.182, reef zone significant 0.001
may_prof = subset(its2profNorm, month == "May")
distprof19 = vegdist(may_prof[8:length(may_prof)], method = "bray")
dispprof19 = betadisper(distprof19, may_prof$bleaching_status)
anova(dispSprof19) #not significant 0.1443
set.seed(123)
adonis(may_prof[, c(8:length(may_prof))]~bleaching_status+reef_zone+bleaching_status:reef_zone, data = may_prof, permutations = 999, method="bray")  
#differences in reef zone in may samples: significant 0.001
dispprof20 = betadisper(distprof19, may_prof$reef_zone)
anova(dispprof20) #significant 1.293e-15
set.seed(123)
adonis(may_prof[, c(8:length(may_prof))]~reef_zone, data = may_prof, permutations = 999, method="bray")  

#NON METRIC MULTIDIMENSIONAL SCALING
#all colonies
set.seed(123)
nmds_allprof = metaMDS(its2profNorm[8:length(its2profNorm)], maxit=9999, trymax = 2000, k = 2, autotransform = F) #solution reached, stress: 0.1162671
nmds_alldf = as.data.frame(scores(nmds_allprof))
nmds_alldf$Sample = its2profNorm$Sample
nmds_alldf$response = its2profNorm$response
nmds_alldf$reef_zone = its2profNorm$reef_zone
nmds_alldf$month = its2profNorm$month

#all recovered and resistant colonies
set.seed(123)
nmds_respprof = metaMDS(its2profNormResp[8:length(its2profNormResp)], maxit=9999, trymax = 10000, k = 2, distance = "bray", autotransform = T) #solution reached, stress: 0.03857115
nmds_respdf = as.data.frame(scores(nmds_respprof))
nmds_respdf$Sample = its2profNormResp$Sample
nmds_respdf$response = its2profNormResp$response
nmds_respdf$reef_zone = its2profNormResp$reef_zone
nmds_respdf$month = its2profNormResp$month

#may crest samples (for bleached vs healthy)
set.seed(123)
nmds_maycrest = metaMDS(may_crest_prof[8:length(may_crest_prof)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = T) #solution reached, stress: 0.07188213
nmds_mcdf = as.data.frame(scores(nmds_maycrest))
nmds_mcdf$Sample = may_crest_prof$Sample
nmds_mcdf$response = may_crest_prof$response
nmds_mcdf$reef_zone = may_crest_prof$reef_zone
nmds_mcdf$month = may_crest_prof$month
nmds_mcdf$bleaching_status = may_crest_prof$bleaching_status

#May samples
mayss = subset(its2profNorm, month == "May")
set.seed(123)
nmds_mayall = metaMDS(mayss[8:length(mayss)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = F) #solution reached, stress: 0.05598264
nmds_mayalldf = as.data.frame(scores(nmds_mayall))
nmds_mayalldf$Sample = mayss$Sample
nmds_mayalldf$response = mayss$response
nmds_mayalldf$reef_zone = mayss$reef_zone
nmds_mayalldf$month = mayss$month
nmds_mayalldf$bleaching_status = mayss$bleaching_status

#October samples
set.seed(123)
octss = subset(its2profNorm, month == "October")
nmds_octall = metaMDS(octss[8:length(octss)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = F) #solution reached, stress: 0.06634123
nmds_octalldf = as.data.frame(scores(nmds_octall))
nmds_octalldf$Sample = octss$Sample
nmds_octalldf$response = octss$response
nmds_octalldf$reef_zone = octss$reef_zone
nmds_octalldf$month = octss$month

#after bleaching samples, all colonies
set.seed(123)
nmds_afterall = metaMDS(after_blALL[8:length(after_blALL)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = F) #solution reached, stress: 0.06634123
nmds_aalldf = as.data.frame(scores(nmds_afterall))
nmds_aalldf$Sample = after_blALL$Sample
nmds_aalldf$response = after_blALL$response
nmds_aalldf$reef_zone = after_blALL$reef_zone
nmds_aalldf$month = after_blALL$month

#backreef samples
set.seed(123)
nmds_backreef = metaMDS(backreef[8:length(backreef)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = F) #solution reached, stress: 0.0.09957208
nmds_backreefdf = as.data.frame(scores(nmds_backreef))
nmds_backreefdf$Sample = backreef$Sample
nmds_backreefdf$response = backreef$response
nmds_backreefdf$reef_zone = backreef$reef_zone
nmds_backreefdf$month = backreef$month

#may backreef samples
set.seed(123)
nmds_br_prof = metaMDS(may_br_prof[8:length(may_br_prof)], distance = "bray", maxit=1000, trymax = 500, k = 2, autotransform = F) #solution reached, stress: 0.137684
nmds_brprofdf = as.data.frame(scores(nmds_br_prof))
nmds_brprofdf$Sample = may_br_prof$Sample
nmds_brprofdf$response = may_br_prof$response
nmds_brprofdf$reef_zone = may_br_prof$reef_zone
nmds_brprofdf$month = may_br_prof$month
nmds_brprofdf$bleaching_status = may_br_prof$bleaching_status

#STRUCTURE PLOTS
#resistant vs recovered
nmds1 = ordiplot(nmds_respdf[1:2], display = "sites", xlab = "NMDS1", ylab = "NMDS2" )
points(nmds1, "sites", pch = 19, col = "#960200", select = nmds_respdf$response == "recovered")
points(nmds1, "sites", pch = 19, col = "#FFD046", select = nmds_respdf$response == "resistant")
ordispider(nmds1, nmds_respdf$response, col = c("#960200","#FFD046"))
ordihull(nmds1, nmds_respdf$response, draw = c("polygon"), col = c("#960200","#FFD046"), alpha=0.2, lwd = 0.0000000001, lty = 0)

#reef zone for all colonies
nmds2 = ordiplot(nmds_alldf[1:2], display = "sites", xlab = "NMDS1", ylab = "NMDS2" )
points(nmds2, "sites", pch = 19, col = "#009A9A", select = nmds_alldf$reef_zone == "backreef")
points(nmds2, "sites", pch = 19, col = "#81C14B", select = nmds_alldf$reef_zone == "crest")
points(nmds2, "sites", pch = 19, col = "#026440", select = nmds_alldf$reef_zone == "forereef")
ordispider(nmds2, nmds_alldf$reef_zone, col = c("#009A9A","#81C14B", "#026440"))
ordihull(nmds2, nmds_alldf$reef_zone, draw = c("polygon"), col = c("#009A9A","#81C14B", "#026440"), alpha=0.2, lwd = 0.0000000001, lty = 0)

#healthy vs bleached, may crest
nmds4 = ordiplot(nmds_mcdf[1:2], display = "sites", xlab = "NMDS1", ylab = "NMDS2" )
points(nmds4, "sites", pch = 19, col = "#BFACC8", select = nmds_mcdf$bleaching_status == "bleached")
points(nmds4, "sites", pch = 19, col = "#2A4470", select = nmds_mcdf$bleaching_status == "healthy")
ordispider(nmds4, nmds_mcdf$bleaching_status, col = c("#BFACC8","#2A4470"))
ordihull(nmds4, nmds_mcdf$bleaching_status, draw = c("polygon"), col = c("#BFACC8","#2A4470"), alpha=0.2, lwd = 0.0000000001, lty = 0 )

#bleaching status for may backreef colonies
nmds19 = ordiplot(nmds_brprofdf[1:2], display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2" )
points(nmds19, "sites", pch = 19, col = "#BFACC8", select = nmds_brprofdf$bleaching_status == "bleached")
points(nmds19, "sites", pch = 19, col = "#2A4470", select = nmds_brprofdf$bleaching_status == "healthy")
ordispider(nmds19, nmds_brprofdf$bleaching_status, col = c("#BFACC8","#2A4470"))
ordihull(nmds19, nmds_brprofdf$bleaching_status, draw = c("polygon"), col = c("#BFACC8","#2A4470"), alpha=0.2, lwd = 0.0000000001, lty = 0 )

#bleaching status for all may colonies, with color as bleaching status and shape as reef zone
nmds_mayalldf$style_zone <- paste(nmds_mayalldf$bleaching_status, nmds_mayalldf$reef_zone, sep="_")
nmds_mayalldf$style_zone = as.factor(nmds_mayalldf$style_zone)
nmds24 = ordiplot(nmds_mayalldf[1:2], display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2" )
points(nmds24, "sites", pch = 1, col = "#BFACC8", select = nmds_mayalldf$style_zone == "bleached_backreef")
points(nmds24, "sites", pch = 2, col = "#BFACC8" , select = nmds_mayalldf$style_zone == "bleached_crest")
points(nmds24, "sites", pch = 1, col = "#2A4470", select = nmds_mayalldf$style_zone == "healthy_backreef")
points(nmds24, "sites", pch = 2, col = "#2A4470", select = nmds_mayalldf$style_zone == "healthy_crest")
ordispider(nmds24, nmds_mayalldf$bleaching_status, col = c("#BFACC8","#2A4470"))
ordihull(nmds24, nmds_mayalldf$bleaching_status, draw = c("polygon"), col = c("#BFACC8","#2A4470"), alpha=0.2, lwd = 0.0000000001, lty = 0)

#reef zone in all after (October) samples
nmds25 = ordiplot(nmds_octalldf[1:2], display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2" )
points(nmds25, "sites", pch = 19, col = "#009A9A", select = nmds_octalldf$reef_zone == "backreef")
points(nmds25, "sites", pch = 19, col = "#81C14B", select = nmds_octalldf$reef_zone == "crest")
points(nmds25, "sites", pch = 19, col = "#026440", select = nmds_octalldf$reef_zone == "forereef")
ordispider(nmds25, nmds_octalldf$reef_zone, col = c("#009A9A","#81C14B","#026440"))
ordihull(nmds25, nmds_octalldf$reef_zone, draw = c("polygon"), col = c("#009A9A","#81C14B","#026440"), alpha=0.2, lwd = 0.0000000001, lty = 0)

#reef zone in all may samples, not very useful bc no deeper forereef samples
nmds30 = ordiplot(nmds_mayalldf[1:2], display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2" )
points(nmds30, "sites", pch = 19, col = "#009A9A", select = nmds_mayalldf$reef_zone == "backreef")  
points(nmds30, "sites", pch = 19, col = "#81C14B", select = nmds_mayalldf$reef_zone == "crest")
points(nmds30, "sites", pch = 19, col = "#026440", select = nmds_mayalldf$reef_zone == "forereef")
ordispider(nmds30, nmds_mayalldf$reef_zone, col = c("#009A9A","#81C14B","#026440"))
ordihull(nmds30, nmds_mayalldf$reef_zone, draw = c("polygon"), col = c("#009A9A","#81C14B","#026440"), alpha=0.2, lwd = 0.0000000001, lty = 0)  

