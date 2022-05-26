# ITS2 SEQUENCE DATA VENN DIAGRAMS #

#clear working environment
rm(list=ls())

library(tidyverse)
library(VennDiagram)
library(eulerr)
library(gridExtra)

setwd("/Users/lumosmaximma/Desktop/coral/its2wd")
seqs = read.csv("seqs.csv")

# IMPORTANT NOTE: the first set of Venn diagrams were generated using ASV data, which was used in preliminary analyses, but 
# not included in the final manuscript. To view code for figures used in the manuscript, scroll down to line 130.

#create lists of ASVs present in each grouping of samples
#by response
resistant = subset(seqs, response == "resistant")
  resistant = resistant[,8:length(resistant)]
  resistanttaxa = as.data.frame(colSums(resistant))
  resistanttaxa = subset(resistanttaxa, colSums(resistant) != 0)
  resistanttaxa = rownames(resistanttaxa)
recovered = subset(seqs, response == "recovered")
  recovered = recovered[,8:length(recovered)]
  recoveredtaxa = as.data.frame(colSums(recovered))
  recoveredtaxa = subset(recoveredtaxa, colSums(recovered) != 0)
  recoveredtaxa = rownames(recoveredtaxa)  

#by reef zone  
crestsub = subset(seqs, reef_zone == "crest")
  crestsub = crestsub[,8:length(crestsub)]
  cresttaxa = as.data.frame(colSums(crestsub))
  cresttaxa = subset(cresttaxa, colSums(crestsub) != 0)
  cresttaxa = rownames(cresttaxa)
forereefsub = subset(seqs, reef_zone == "forereef")
  forereefsub = forereefsub[,8:length(forereefsub)]
  forereeftaxa = as.data.frame(colSums(forereefsub))
  forereeftaxa = subset(forereeftaxa, colSums(forereefsub) != 0)
  forereeftaxa = rownames(forereeftaxa)
backreefsub = subset(seqs, reef_zone == "backreef")
  backreefsub = backreefsub[,8:length(backreefsub)]
  backreeftaxa = as.data.frame(colSums(backreefsub))
  backreeftaxa = subset(backreeftaxa, colSums(backreefsub) != 0)
  backreeftaxa = rownames(backreeftaxa)
  
#by bleaching status in May, crest only
may_crest_subh = subset(seqs, month == "May" & reef_zone == "crest" & bleaching_status == "healthy")
  may_crest_subh = may_crest_subh[,8:length(may_crest_subh)]
  may_crest_subhtaxa = as.data.frame(colSums(may_crest_subh))
  may_crest_subhtaxa = subset(may_crest_subhtaxa, colSums(may_crest_subh) != 0)
  may_crest_subhtaxa = rownames(may_crest_subhtaxa)
may_crest_subb = subset(seqs, month == "May" & reef_zone == "crest" & bleaching_status == "bleached")
  may_crest_subb = may_crest_subb[,8:length(may_crest_subb)]
  may_crest_subbtaxa = as.data.frame(colSums(may_crest_subb))
  may_crest_subbtaxa = subset(may_crest_subbtaxa, colSums(may_crest_subb) != 0)
  may_crest_subbtaxa = rownames(may_crest_subbtaxa)

#by bleaching status in May, backreef only  
may_br_subh = subset(seqs, month == "May" & reef_zone == "backreef" & bleaching_status == "healthy")
  may_br_subh = may_br_subh[,8:length(may_br_subh)]
  may_br_subhtaxa = as.data.frame(colSums(may_br_subh))
  may_br_subhtaxa = subset(may_br_subhtaxa, colSums(may_br_subh) != 0)
  may_br_subhtaxa = rownames(may_br_subhtaxa)
may_br_subb = subset(seqs, month == "May" & reef_zone == "backreef" & bleaching_status == "bleached")
  may_br_subb = may_br_subb[,8:length(may_br_subb)]
  may_br_subbtaxa = as.data.frame(colSums(may_br_subb))
  may_br_subbtaxa = subset(may_br_subbtaxa, colSums(may_br_subb) != 0)
  may_br_subbtaxa = rownames(may_br_subbtaxa)  
  
#by bleaching status in May, backreef only  
may_subh = subset(seqs, month == "May" & bleaching_status == "healthy")
  may_subh = may_subh[,8:length(may_subh)]
  may_subhtaxa = as.data.frame(colSums(may_subh))
  may_subhtaxa = subset(may_subhtaxa, colSums(may_subh) != 0)
  may_subhtaxa = rownames(may_subhtaxa)
may_subb = subset(seqs, month == "May" & bleaching_status == "bleached")
  may_subb = may_subb[,8:length(may_subb)]
  may_subbtaxa = as.data.frame(colSums(may_subb))
  may_subbtaxa = subset(may_subbtaxa, colSums(may_subb) != 0)
  may_subbtaxa = rownames(may_subbtaxa) 
  
#recovered/resistant Venn diagram
recov.resis = length(intersect(recoveredtaxa, resistanttaxa))
recov.only = length(recoveredtaxa) - recov.resis  
resis.only = length(resistanttaxa) - recov.resis
responsevenn = euler(c("recovered" = (recov.only), "resistant" = (resis.only), "recovered&resistant" = (recov.resis)))
plot(responsevenn, quantities = TRUE, col = c("#960200","#FFD046"), fill = c("#960200","#FFD046"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(resistanttaxa,(intersect(recoveredtaxa,resistanttaxa)))

#reef zones Venn diagram
crest.back = intersect(cresttaxa, backreeftaxa)
crest.back.fore = intersect(crest.back, forereeftaxa)
crest.back.only = length(crest.back) - length(crest.back.fore)
crest.back.fore.length = length(crest.back.fore)
crest.fore = intersect(cresttaxa, forereeftaxa)
crest.fore.only = length(crest.fore) - length(crest.back.fore)
back.fore = intersect(backreeftaxa, forereeftaxa)
back.fore.only = length(back.fore) - length(crest.back.fore)
crest = length(cresttaxa) - crest.back.only - crest.fore.only - crest.back.fore.length
backreef = length(backreeftaxa) - back.fore.only - crest.back.only - crest.back.fore.length
forereef = length(forereeftaxa) - back.fore.only - crest.fore.only - crest.back.fore.length
zonevenn = euler(c("crest" = (crest), "backreef" = (backreef), "forereef" = (forereef), "crest&backreef" = (crest.back.only), "crest&forereef" = (crest.fore.only), "backreef&forereef" = (back.fore.only), "crest&backreef&forereef" = (crest.back.fore.length)))
plot(zonevenn, quantities = TRUE, col = c("#81C14B","#009A9A", "#026440"), fill = c("#81C14B","#009A9A", "#026440"), font=1, cex=1, alpha=0.5, lwd=c(2,2,2))

#bleached vs healthy May crest
health.bleach = length(intersect(may_crest_subbtaxa, may_crest_subhtaxa))
health.only = length(may_crest_subhtaxa) - health.bleach
bleach.only = length(may_crest_subbtaxa) - health.bleach
maycrestvenn = euler(c("healthy" = (health.only), "bleached" = (bleach.only), "healthy&bleached" = (health.bleach)))
plot(maycrestvenn, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(may_crest_subhtaxa,(intersect(may_crest_subbtaxa, may_crest_subhtaxa)))

#bleached vs healthy May backreef
health.bleachbr = length(intersect(may_br_subbtaxa, may_br_subhtaxa))
health.onlybr = length(may_br_subhtaxa) - health.bleachbr
bleach.onlybr = length(may_br_subbtaxa) - health.bleachbr
maybrvenn = euler(c("healthy" = (health.onlybr), "bleached" = (bleach.onlybr), "healthy&bleached" = (health.bleachbr)))
plot(maybrvenn, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(may_br_subhtaxa,(intersect(may_br_subbtaxa, may_br_subhtaxa)))

#bleached vs healthy May 
health.bleachall = length(intersect(may_subbtaxa, may_subhtaxa))
health.onlyall = length(may_subhtaxa) - health.bleachall
bleach.onlyall = length(may_subbtaxa) - health.bleachall
mayallvenn = euler(c("healthy" = (health.onlyall), "bleached" = (bleach.onlyall), "healthy&bleached" = (health.bleachall)))
plot(mayallvenn, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))


# VENN DIAGRAMS WITH COLLAPSED PROFILES #

# NOTE: these are the visualizations included in the manuscript.

#clear working environment
rm(list=ls())

library(tidyverse)
library(VennDiagram)
library(eulerr)
library(gridExtra)

setwd("/Users/lumosmaximma/Desktop/coral/its2wd/profiles")
profs = read.csv("profs.csv")
profs = profs[,-15]

#create lists of profiles present in each grouping of samples
#by response
resistant2 = subset(profs, response == "resistant")
resistant2 = resistant2[,8:length(resistant2)]
resistanttaxa2 = as.data.frame(colSums(resistant2))
resistanttaxa2 = subset(resistanttaxa2, colSums(resistant2) != 0)
resistanttaxa2 = rownames(resistanttaxa2)
recovered2 = subset(profs, response == "recovered")
recovered2 = recovered2[,8:length(recovered2)]
recoveredtaxa2 = as.data.frame(colSums(recovered2))
recoveredtaxa2 = subset(recoveredtaxa2, colSums(recovered2) != 0)
recoveredtaxa2 = rownames(recoveredtaxa2)  

#by reef zone  
crestsub2 = subset(profs, reef_zone == "crest")
crestsub2 = crestsub2[,8:length(crestsub2)]
cresttaxa2 = as.data.frame(colSums(crestsub2))
cresttaxa2 = subset(cresttaxa2, colSums(crestsub2) != 0)
cresttaxa2 = rownames(cresttaxa2)
forereefsub2 = subset(profs, reef_zone == "forereef")
forereefsub2 = forereefsub2[,8:length(forereefsub2)]
forereeftaxa2 = as.data.frame(colSums(forereefsub2))
forereeftaxa2 = subset(forereeftaxa2, colSums(forereefsub2) != 0)
forereeftaxa2 = rownames(forereeftaxa2)
backreefsub2 = subset(profs, reef_zone == "backreef")
backreefsub2 = backreefsub2[,8:length(backreefsub2)]
backreeftaxa2 = as.data.frame(colSums(backreefsub2))
backreeftaxa2 = subset(backreeftaxa2, colSums(backreefsub2) != 0)
backreeftaxa2 = rownames(backreeftaxa2)

#by bleaching status in May, crest only
may_crest_subh2 = subset(profs, month == "May" & reef_zone == "crest" & bleaching_status == "healthy")
may_crest_subh2 = may_crest_subh2[,8:length(may_crest_subh2)]
may_crest_subhtaxa2 = as.data.frame(colSums(may_crest_subh2))
may_crest_subhtaxa2 = subset(may_crest_subhtaxa2, colSums(may_crest_subh2) != 0)
may_crest_subhtaxa2 = rownames(may_crest_subhtaxa2)
may_crest_subb2 = subset(profs, month == "May" & reef_zone == "crest" & bleaching_status == "bleached")
may_crest_subb2 = may_crest_subb2[,8:length(may_crest_subb2)]
may_crest_subbtaxa2 = as.data.frame(colSums(may_crest_subb2))
may_crest_subbtaxa2 = subset(may_crest_subbtaxa2, colSums(may_crest_subb2) != 0)
may_crest_subbtaxa2 = rownames(may_crest_subbtaxa2)

#by bleaching status in May, backreef only  
may_br_subh2 = subset(profs, month == "May" & reef_zone == "backreef" & bleaching_status == "healthy")
may_br_subh2 = may_br_subh2[,8:length(may_br_subh2)]
may_br_subhtaxa2 = as.data.frame(colSums(may_br_subh2))
may_br_subhtaxa2 = subset(may_br_subhtaxa2, colSums(may_br_subh2) != 0)
may_br_subhtaxa2 = rownames(may_br_subhtaxa2)
may_br_subb2 = subset(profs, month == "May" & reef_zone == "backreef" & bleaching_status == "bleached")
may_br_subb2 = may_br_subb2[,8:length(may_br_subb2)]
may_br_subbtaxa2 = as.data.frame(colSums(may_br_subb2))
may_br_subbtaxa2 = subset(may_br_subbtaxa2, colSums(may_br_subb2) != 0)
may_br_subbtaxa2 = rownames(may_br_subbtaxa2)  

#by bleaching status in May, backreef only  
may_subh2 = subset(profs, month == "May" & bleaching_status == "healthy")
may_subh2 = may_subh2[,8:length(may_subh2)]
may_subhtaxa2 = as.data.frame(colSums(may_subh2))
may_subhtaxa2 = subset(may_subhtaxa2, colSums(may_subh2) != 0)
may_subhtaxa2 = rownames(may_subhtaxa2)
may_subb2 = subset(profs, month == "May" & bleaching_status == "bleached")
may_subb2 = may_subb2[,8:length(may_subb2)]
may_subbtaxa2 = as.data.frame(colSums(may_subb2))
may_subbtaxa2 = subset(may_subbtaxa2, colSums(may_subb2) != 0)
may_subbtaxa2 = rownames(may_subbtaxa2) 

#recovered/resistant Venn diagram
recov.resis2 = length(intersect(recoveredtaxa2, resistanttaxa2))
recov.only2 = length(recoveredtaxa2) - recov.resis2
resis.only2 = length(resistanttaxa2) - recov.resis2
responsevenn2 = euler(c("recovered" = (recov.only2), "resistant" = (resis.only2), "recovered&resistant" = (recov.resis2)))
plot(responsevenn2, quantities = TRUE, col = c("#960200","#FFD046"), fill = c("#960200","#FFD046"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(resistanttaxa2,(intersect(recoveredtaxa2,resistanttaxa2)))

#reef zones Venn diagram
crest.back2 = intersect(cresttaxa2, backreeftaxa2)
crest.back.fore2 = intersect(crest.back2, forereeftaxa2)
crest.back.only2 = length(crest.back2) - length(crest.back.fore2)
crest.back.fore.length2 = length(crest.back.fore2)
crest.fore2 = intersect(cresttaxa2, forereeftaxa2)
crest.fore.only2 = length(crest.fore2) - length(crest.back.fore2)
back.fore2 = intersect(backreeftaxa2, forereeftaxa2)
back.fore.only2 = length(back.fore2) - length(crest.back.fore2)
crest2 = length(cresttaxa2) - crest.back.only2 - crest.fore.only2 - crest.back.fore.length2
backreef2 = length(backreeftaxa2) - back.fore.only2 - crest.back.only2 - crest.back.fore.length2
forereef2 = length(forereeftaxa2) - back.fore.only2 - crest.fore.only2 - crest.back.fore.length2
zonevenn2 = euler(c("crest" = (crest2), "backreef" = (backreef2), "forereef" = (forereef2), "crest&backreef" = (crest.back.only2), "crest&forereef" = (crest.fore.only2), "backreef&forereef" = (back.fore.only2), "crest&backreef&forereef" = (crest.back.fore.length2)))
plot(zonevenn2, quantities = TRUE, col = c("#81C14B","#009A9A", "#026440"), fill = c("#81C14B","#009A9A", "#026440"), font=1, cex=1, alpha=0.5, lwd=c(2,2,2))

#bleached vs healthy May crest
health.bleach2 = length(intersect(may_crest_subbtaxa2, may_crest_subhtaxa2))
health.only2 = length(may_crest_subhtaxa2) - health.bleach2
bleach.only2 = length(may_crest_subbtaxa2) - health.bleach2
maycrestvenn2 = euler(c("healthy" = (health.only2), "bleached" = (bleach.only2), "healthy&bleached" = (health.bleach2)))
plot(maycrestvenn2, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(may_crest_subhtaxa2,(intersect(may_crest_subbtaxa2, may_crest_subhtaxa2)))

#bleached vs healthy May backreef
health.bleachbr2 = length(intersect(may_br_subbtaxa2, may_br_subhtaxa2))
health.onlybr2 = length(may_br_subhtaxa2) - health.bleachbr2
bleach.onlybr2 = length(may_br_subbtaxa2) - health.bleachbr2
maybrvenn2 = euler(c("healthy" = (health.onlybr2), "bleached" = (bleach.onlybr2), "healthy&bleached" = (health.bleachbr2)))
plot(maybrvenn2, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))
setdiff(may_br_subhtaxa2,(intersect(may_br_subbtaxa2, may_br_subhtaxa2)))

#bleached vs healthy May 
health.bleachall2 = length(intersect(may_subbtaxa2, may_subhtaxa2))
health.onlyall2 = length(may_subhtaxa2) - health.bleachall2
bleach.onlyall2 = length(may_subbtaxa2) - health.bleachall2
mayallvenn2 = euler(c("healthy" = (health.onlyall2), "bleached" = (bleach.onlyall2), "healthy&bleached" = (health.bleachall2)))
plot(mayallvenn2, quantities = TRUE, col = c("#2A4470","#BFACC8"), fill = c("#2A4470","#BFACC8"), font = 1, cex = 1, alpha = 0.5, lwd = c(2,2,2))


