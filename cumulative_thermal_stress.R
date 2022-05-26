#clear working environment
rm(list=ls())

library(plyr)
library(dplyr)
library(lubridate)
library(scales)
library(runner)

setwd("/Users/lumosmaximma/Desktop/coral/MooreaTemp")

#loading data, large files, can be downloaded from the Moorea Coral Reef LTER website, 
#http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.1035
lter0old = read.csv("MCR_LTER00_BottomMountThermistors_2020.csv", header=TRUE)
lter0new = read.csv("LTER00_BottomMountThermistors_NEW.csv", header=TRUE)
lter0 = rbind(lter0old,lter0new)
lter1 = read.csv("MCR_LTER01_BottomMountThermistors_2020.csv", header=TRUE)
lter2old = read.csv("MCR_LTER02_BottomMountThermistors_2020.csv", header=TRUE)
lter2new = read.csv("LTER02_BottomMountThermistors_NEW.csv", header=TRUE)
lter2 = rbind(lter2old, lter2new)
lter3 = read.csv("MCR_LTER03_BottomMountThermistors_2020.csv", header=TRUE)
lter4 = read.csv("MCR_LTER04_BottomMountThermistors_2020.csv", header=TRUE)
lter5 = read.csv("MCR_LTER05_BottomMountThermistors_2020.csv", header=TRUE)
lter6 = read.csv("MCR_LTER06_BottomMountThermistors_2020.csv", header=TRUE)

#combining all data and subsetting for forereef at 10 meters depth
thermal = rbind(lter0, lter1, lter2, lter3, lter4, lter5, lter6)
thermal = subset(thermal, thermal$reef_type_code=="FOR")
thermal$sensor_type = NULL
thermal$sensor_depth_m = as.factor(thermal$sensor_depth_m)
thermal = subset(thermal, thermal$sensor_depth_m=="10")

#formatting dates
thermal$date_use = ymd_hms(thermal$time_local)
thermal$day_mo_yr = format(thermal$date_use, '%Y-%m-%d')
thermal$day = format(thermal$date_use, '%d')
thermal$month = format(thermal$date_use, '%m')
thermal$year = format(thermal$date_use, '%Y')

thermal$day = as.factor(thermal$day)
thermal$month = as.factor(thermal$month)
thermal$year = as.factor(thermal$year)
thermal$date_use = as.Date(thermal$date_use, "%m-%d-%Y")

#selecting date range
thermal_2019 = subset(thermal, date_use>="2018-09-30" & date_use<="2019-10-31")
thermal_2019_mean = ddply(thermal_2019, .(day, month, year), summarize, mean_daily_temp = mean(temperature_c), sd=NA)
thermal_2019_mean$date = as.Date(with(thermal_2019_mean, paste(month, day, year,sep="-")), format="%m-%d-%Y")

#calculating accumulated heat stress
mma_ref = 29 #value based on Pratchett et al., 2013
thermal_2019_mean$hotspot = thermal_2019_mean$mean_daily_temp - mma_ref
thermal_2019_mean$hotspot[thermal_2019_mean$hotspot < 0] = 0
thermal_2019_mean$cumstress = NA
thermal_2019_mean = thermal_2019_mean %>% arrange(ymd(thermal_2019_mean$date))
thermal_2019_mean = thermal_2019_mean %>%
  mutate(
    cumstress = sum_run(
      x = thermal_2019_mean$hotspot, 
      k = 84, 
      idx = as.Date(thermal_2019_mean$date, format = "%Y/%m/%d"))
  )

#fine-tuning date range and formatting
thermal_2019_mean = subset(thermal_2019_mean, date>="2018-11-01" & date<="2019-10-31")
thermal_2019_mean$date = as.POSIXct(thermal_2019_mean$date, tz = "Pacific/Tahiti")

#PLOTTING THERMAL STRESS
#assigning color palette for heat colors
heat_palette = colorRampPalette(c("#ffffff", "#ffffc1", "#ffff84", "#ffff46", "#ffff09", "#ffed00",
                                  "#ffd200", "#ffb700", "#ffa700", "#ff9700", "#ff8600", "#ff7600",
                                  "#ff6600", "#ff5500", "#ff4500", "#ff3300", "#ff2400", "#ff1100",
                                  "#ff0400", "#ea0000", "#d30000", "#bc0000", "#a40000", "#8c0000",
                                  "#750000", "#5e0000", "#460000", "#2f0000", "#170000", "#000000"))

dhw_floor = floor(thermal_2019_mean$cumstress) + 1 
heat_colorspal = heat_palette(40)
dhw_color = heat_colorspal[dhw_floor]

from = thermal_2019_mean$date - 1.75*86400
to = thermal_2019_mean$date + 1.75*86400

polyCurve = function(x, y, from, to, n = 50, miny, col = "red", border = col) {
  drawPoly = function(fun, from, to, n = 50, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n)
    polygon(x = c(Sq[1], Sq, Sq[n]),
            y = c(miny, fun(Sq), miny),
            col = col, border = border)
  }
  lf = length(from)
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col = rep(col, length.out = lf)
  if(length(border) != lf)
    border = rep(border, length.out = lf)
  if(missing(miny))
    miny = min(y)
  interp = approxfun(x = x, y = y)
  mapply(drawPoly, from = from, to = to, col = col, border = border,
         MoreArgs = list(fun = interp, n = n, miny = miny))
  invisible()
}

#set dates for sampling timepoints
may_sample = as.POSIXct("2019-05-21 00:00:00", tz = "Pacific/Tahiti", format="%Y-%m-%d %H:%M:%S")
aug_sample = as.POSIXct("2019-08-02 00:00:00", tz = "Pacific/Tahiti", format="%Y-%m-%d %H:%M:%S")
oct_sample = as.POSIXct("2019-10-06 00:00:00", tz = "Pacific/Tahiti", format="%Y-%m-%d %H:%M:%S")

#plotting thermal stress curve with colors
plot(thermal_2019_mean$date, thermal_2019_mean$cumstress, type="l", 
                   xlab="", ylab="", ylim=c(0,40), cex.axis=1, cex.lab=1.2,
                   yaxs="i", xaxs="i", lwd=0.5, xaxt='n', yaxt='n', col="gray40",
                   panel.last = c(polyCurve(thermal_2019_mean$date, thermal_2019_mean$cumstress,
                                            from = from, to = to, miny = 0, col = dhw_color),
                                  abline(v = may_sample, col = "red", lwd = 2.1, lty = 2),
                                  abline(v = aug_sample, col = "red", lwd = 2.1, lty = 2),
                                  abline(v = oct_sample, col = "red", lwd = 2.1, lty = 2)))

# Add in time axes (multiple axes added to allow for customization)
axis.POSIXct(side=1,at=seq(thermal_2019_mean$date[1],thermal_2019_mean$date[length(thermal_2019_mean$date)],
                           by="month"),tck=-0.03,
             cex.axis=0.93,labels=c("Nov","Dec","Jan","Feb","Mar","Apr","May",
                                    "Jun","Jul","Aug","Sep","Oct"),
             lwd.ticks=1.5,padj=0)
axis.POSIXct(side=1,thermal_2019_mean$date,cex.axis=0.93,tck=0,padj=-1.5, label = FALSE)

Y = c(0, 5, 10, 15, 20, 25, 30, 35, 40)

axis(side = 2, at = Y, cex.axis = 0.93, tck = -0.02, lwd.ticks = 1.5, las = 2, hadj = 0)

#plotting temperature data to thermal stress curve
par(new=T)
plot(thermal_2019_mean$mean_daily_temp, type = 'l', col = "darkgray",
     ylim = c(25, 31.5), lwd = 2.5,
     xlab = "", ylab = "", xaxs = "i", xaxt = "n", yaxt = "n") 
abline(29, 0, col = "black", lwd = 2.5, lty = 1)
title(ylab="Cumulative heat stress (°C-day)", line = 1.2, cex.lab = 1.1)

Z = c(25, 26, 27, 28, 29, 30, 31) 
axis(side = 4, at = Z, cex.axis = 0.93, tck = -0.02, lwd.ticks = 1.5, las = 2, hadj = 0)
mtext("Temperature (°C)",side=4, cex=0.75,line=1.25)

#6.91 x 3.56 in

