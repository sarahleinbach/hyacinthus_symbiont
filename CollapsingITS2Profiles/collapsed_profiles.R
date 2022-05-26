# COLLAPSING ITS2 PROFILES VISUALIZATION #

# These graphs were generated using Bray-Curtis dissimilarities provided as a direct output of 
# SymPortal and subsequently used to inform the decision of how to collapse ITS2 profiles. #C 

#clear working environment
rm(list=ls())

library(tidyverse)

setwd("/Users/lumosmaximma/Desktop/coral/its2wd/profiles")

#importing PCoA coordinates for each clade, csvs output from SymPortal
Acoords = read.csv("Acoords.csv")
Acoords = Acoords[-6,-2]
Ccoords = read.csv("Ccoords.csv")
Ccoords = Ccoords[-11,-2]
Dcoords = read.csv("Dcoords.csv")
Dcoords = Dcoords[-7,-2]

#Symbiodinium profile collapse
Apcoa = ggplot(data = Acoords)+
  geom_point(mapping = aes(x = PC1, y = PC2), size = 3)+
  geom_text(mapping = aes(x = PC1, y = PC2), label = Acoords$sample)+
  labs(x = "PC1 (92.9%)", y = "PC2 (4.4%)")+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
Apcoa
ggsave("Apcoa.pdf", width=4.5, height=4.5, units = "in")

#Cladocopium profile collapse
Cpcoa = ggplot(data = Ccoords)+
  geom_point(mapping = aes(x = PC1, y = PC2), size = 3)+
  geom_text(mapping = aes(x = PC1, y = PC2), label = Ccoords$sample)+
  labs(x = "PC1 (33.1%)", y = "PC2 (19.5%)")+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
Cpcoa
ggsave("Cpcoa.pdf", width=4.5, height=4.5, units = "in")

#Durusdinium profile collapse
Dpcoa = ggplot(data = Dcoords)+
  geom_point(mapping = aes(x = PC1, y = PC2), size = 3)+
  geom_text(mapping = aes(x = PC1, y = PC2), label = Dcoords$sample)+
  labs(x = "PC1 (70.4%)", y = "PC2 (18.6%)")+
  theme_classic(base_size = 12)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
Dpcoa
ggsave("Dpcoa.pdf", width=4.5, height=4.5, units = "in")

