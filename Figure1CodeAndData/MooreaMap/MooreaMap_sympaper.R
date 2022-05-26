library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggspatial)
library(marmap)
library(raster)
library(ggstar)

setwd("/Users/lumosmaximma/Desktop/coral/MooreaMap")
map = st_read("spatialdata/PYF_adm0.shp")
bathy = read.table("spatialdata/bathy_LIDAR_Riegl_820_05m.xyz")
reg = griddify(bathy, nlon =1000, nlat = 1000)
bath.p = rasterToPoints(reg)
bath.df = data.frame(bath.p)
colnames(bath.df) = c("x", "y", "alt")

sitepm = data.frame(longitude = -149.8177, latitude = -17.4731)
sites = data.frame(longitude = c(-149.8177,-149.8170), latitude = c(-17.4731,-17.4751))

#moorea island map with the surrounding reef
mooreareef_sym = ggplot(map)+
  geom_point(data=subset(bath.df, bath.df$alt<3 & bath.df$alt>0), aes(x,y),color = grey(0.92), size = 0.01, show.legend=FALSE) +
  geom_sf(data = map, color = "black", fill = "white") +
  #geom_point(data=sitepm, aes(x = longitude, y = latitude), size = 5)+
  coord_sf(xlim = c(-149.94,-149.74), ylim = c(-17.61,-17.46))+
  annotation_scale(location="br",width_hint=0.3)+
  annotation_north_arrow(location="br",which_north="true",
                         pad_x = unit(0.7,"in"),
                         pad_y = unit(0.25,"in"),
                         style = north_arrow_fancy_orienteering(),
                         height = unit(1.2, "cm"),
                         width = unit(1.2, "cm"))+
  theme_classic(base_size=16)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill=NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
mooreareef_sym
#ggsave("mooreareef_sym.pdf", width=8, height=7.6, units="in")

pacificblack = ggplot()+
  borders("world2", fill = "white", col= "black")+
  theme_classic(base_size=11)+
  coord_cartesian(xlim = c(110, 290), ylim = c(-30,30))+
  geom_star(aes(x = 210, y = -17), size = 8, color = "red", fill = "red")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(colour = "black", size=1))
pacificblack
ggsave("pacificblack.pdf", width=4, height=2, units="in")


