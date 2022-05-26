# BLEACHING OVER TIME GRAPHS #

#clear working environment
rm(list=ls())

setwd("/Users/lumosmaximma/Desktop/coral/its2wd")
library(tidyverse)

scores = read.csv("BleachingScores.csv")
scores$month = ordered(scores$month, levels = c("May", "August", "October"))
scores$color_bin=as.factor(scores$color_bin)
scores$month=as.factor(scores$month)

scores$color_bin=as.numeric(scores$color_bin)
scoresbin = scores %>% group_by(month) %>% count(color_bin) %>% as.data.frame()
scoresbin = scoresbin %>% add_row(month = "October", color_bin = 1, n = 0)
scoresbin = scoresbin %>% add_row(month = "October", color_bin = 2, n = 0)
scoresbin = scoresbin %>% group_by(month) %>% mutate(freq = round(n/sum(n),3))

#bar chart alternative to the above
zone.labs <- c("backreef", "shallow forereef", "deeper forereef")
names(zone.labs) <- c("backreef", "crest", "forereef")
barmonthzone <- ggplot(data = scores, aes(x = color_bin, fill = month, stat = "stack")) + 
  geom_bar(aes(y = ..prop..), position = "stack", width = 0.8)+ #alpha=0.4
  facet_grid(. ~ reef_zone, scale = "free_y", labeller = labeller(reef_zone = zone.labs))+
  theme_classic(base_size=12)+
  labs(x = "Bleaching Score", y = "Proportion of Colonies")+
  scale_fill_manual(name = "Month", values = c("#71649D","#ADD8E6","#E59F2E"))+
  scale_color_manual(name = "Month", values = c("#71649D","#ADD8E6","#E59F2E"))+ # alternate color palette "#6CBB57", "#419B9A", "#E59F2E"
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
barmonthzone
ggsave("barmonthzone.pdf", width=10, height=3.2, units="in")
