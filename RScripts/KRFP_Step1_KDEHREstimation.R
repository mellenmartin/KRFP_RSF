#### script for estimating fisher homeranges using a KDE with reference bandwith plugin approach ####
#### based on script written by W. David Walter and available at the Walter lab website
#### https://ecosystems.psu.edu/research/labs/walter-lab/manual/home-range-estimation/kernel-density-estimation-kde-with-reference-bandwidth-selection-href

library(ks)
library(rgdal)
library(maptools)
library(gpclib)
library(PBSmapping)
library(sf)
library(dplyr)
library(sp)
library(gridExtra)
library(tidyverse)
library(lubridate)


rm(list=ls());gc() #clear the memory
# set working drive
setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/")

##########################################################################################################################################
##########################################################################################################################################
################################################ Format location data and combine individual data ########################################
##### Load and Format data -- do any subsetting and clean you need to here ####

# load csv data #
all_locs <- read.csv("KRFP_AllLocations_210525.csv",header=T, na.strings = c("","NA"))

# create a column separating all_locs ID by year to allow for easier estimates of annual home ranges
all_locs$ID_year <- paste0(all_locs$FisherID,"_",all_locs$FisherYear)

# create datetime for home range analyses --
# important for assessing and incorporating spatiotemporal autocorrelation
# in aKDE estimates or for subsetting GPS data
all_locs$datetime <- paste(all_locs$Date, all_locs$Time24)
all_locs$datetime_format <- as.POSIXct(strptime(all_locs$datetime, "%m/%d/%Y %H:%M:%S"),tz="PST8PDT")

# because we're going to do a simple KDE approach here, want to subset GPS data so we have similar
# autocorrelation across data sets and individuals
#convert to a character and use strptime again
all_locs$TimeDiff<-c(as.numeric(difftime(all_locs$datetime_format[2:length(all_locs$datetime_format)],all_locs$datetime_format[1:(length(all_locs$datetime_format)-1)],units="days")),NA)
all_locs$julian <- yday(all_locs$datetime_format)


# subset to remove locations w/ incorrect dates,  unknown individuals, and juvenile locations
all_locs_sub2 <- filter(all_locs, DoubleCheck != "Yes") 
all_locs_sub2 <- filter(all_locs_sub2,  SexFull != "Unknown")
all_locs_sub2 <- filter(all_locs_sub2, AgeTransition != "Juvenile" )
all_locs_sub2 <- filter(all_locs_sub2, Type != "Den_same" )
all_locs_sub2 <- filter(all_locs_sub2, Type != "GPS" )


# subset to separate dens
dens <- filter(all_locs, Type == "Den_same")

# subset to separate dens
gps <- filter(all_locs, Type == "GPS")

# subsample GPS data 
gps_sub <- gps %>% 
    group_by(ID_year, julian) %>% 
    do(sample_n(.,1))
  
# subsample den locations
dens_sub <- dens %>% 
  group_by(ID_year, ID) %>% 
  filter(row_number() %% 5 == 1)

# re join the non-den locs and re-used den locs
all_locs_sub2 <- bind_rows(all_locs_sub2, gps_sub)
all_locs_sub2 <- bind_rows(all_locs_sub2, dens_sub)

# subset to only individuals with >= 30 locations. Might want to re-run with 
all_locs_sub2 <- all_locs_sub2 %>%
  group_by(ID_year) %>%
  filter(n() >= 30)

all_locs_sub_male <- filter(all_locs_sub2, SexFull == "Male")

all_locs_sub_female2 <- filter(all_locs_sub2, SexFull == "Female")
all_locs_sub_female2 <- filter(all_locs_sub_female2, Type != "Den_same")

all_locs_sub_female2 <- all_locs_sub_female2 %>% 
  group_by(ID_year, julian) %>% 
  do(sample_n(.,1))

all_locs_sub_female2 <- all_locs_sub_female2 %>%
  group_by(ID_year) %>%
  filter(n() >= 35)

all_locs_sub <- bind_rows(all_locs_sub_male, all_locs_sub_female2)

write.csv(all_locs_sub, "KRFP_FilteredHRLocations_210628.csv",row.names = FALSE)
##### KDE estimation w/ plug-in bandwidth selection ####

# create a simplified dataframe for KDE estimation, selecting just the ID_year and x,y coords
indata <- all_locs_sub[, c("ID_year", "easting_NAD83", "northing_NAD83")]

# unique ID
innames <- unique(indata$ID_year)
outnames <- innames

output <- as.data.frame(matrix(0, nrow = length(innames), ncol = 10))
colnames(output) <- c("ID_year", "noFixes", "h11", "h12", "h21", "h22", "iso30areaKm", "iso50areaKm", "iso90areaKm", "iso95areaKm")

# define the "levels" (i.e., contours) of the HR that you want to derive #
levels <- c(30,50,90,95)

# set output directory for KDE loop to write files to #
dirout <- "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/FisherHR_KDEPolygons_pluginBandwidth_210628"

for(i in 1:length(innames)){
  data <- indata[which(indata$ID_year == innames[i]),]
  if(dim(data)[1] != 0){
    data.xy = data[c("easting_NAD83", "northing_NAD83")]
    coordinates (data.xy) <- ~easting_NAD83+northing_NAD83
    sppt <- SpatialPointsDataFrame(coordinates(data.xy), data)
    sppt <- st_as_sf(sppt, coords = c("easting_NAD83", "northing_NAD83"), crs = "+proj=utm +datum=NAD83 +zone=11 +units=m")
    ##proj4string(sppt) <- CRS("+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
    #st_write(sppt, dsn = paste(dirout,outnames[i], sep = "/"), factor2char = TRUE)
    
    output$ID_year[i] <- as.character(outnames[i])
    locs <- cbind(data$easting_NAD83, data$northing_NAD83)
    try(HpiOut <- Hpi(locs, pilot = "samse", binned = TRUE))
    if(is.null(HpiOut) == "FALSE"){
      output$h11[i] <- HpiOut[1,1]
      output$h12[i] <- HpiOut[1,2]
      output$h21[i] <- HpiOut[2,1]
      output$h22[i] <- HpiOut[2,2]
      fhatOut <- kde(x = locs, H = HpiOut)
      
    }
    if(is.null(fhatOut) == "FALSE"){
      for(j in 1:length(levels)){
        fhat.contlev <- contourLevels(fhatOut, cont=c(levels[j]))
        fhat.contlines <- contourLines(x = fhatOut$eval.points[[1]], y = fhatOut$eval.points[[2]], z = fhatOut$estimate, level = fhat.contlev)
        sldf <- ContourLines2SLDF(fhat.contlines)
        proj4string(sldf) <- CRS("+proj=utm +datum=NAD83 +zone=11 +units=m")
        ps <- SpatialLines2PolySet(sldf)
        #attr(ps, "projection") <- "UTM-Zone15"
        sp <- PolySet2SpatialPolygons(ps)
        dataframe <- as.data.frame(matrix(as.character(1, nrow = 1, ncol = 1)))
        spdf <- SpatialPolygonsDataFrame(sp, dataframe, match.ID = TRUE)
        
        if(j == 1){
          pls <- slot(spdf, "polygons")[[1]]
          gpclibPermit()
          xx <- checkPolygonsHoles(pls)
          a <- sapply(slot(xx, "Polygons"), slot, "area")
          h <- sapply(slot(xx, "Polygons"), slot, "hole")
          output$iso30areaKm[i] <- sum(ifelse(h, -a, a))/1000000
          writeOGR(spdf, dirout, paste(outnames[i], "KUD30", sep=""), "ESRI Shapefile")}
        if(j == 2){
          pls <- slot(spdf, "polygons")[[1]]
          gpclibPermit()
          xx <- checkPolygonsHoles(pls)
          a <- sapply(slot(xx, "Polygons"), slot, "area")
          h <- sapply(slot(xx, "Polygons"), slot, "hole")
          output$iso50areaKm[i] <- sum(ifelse(h, -a, a))/1000000
          writeOGR(spdf, dirout, paste(outnames[i], "KUD50", sep=""), "ESRI Shapefile")}
        if(j == 3){
          pls <- slot(spdf, "polygons")[[1]]
          gpclibPermit()
          xx <- checkPolygonsHoles(pls)
          a <- sapply(slot(xx, "Polygons"), slot, "area")
          h <- sapply(slot(xx, "Polygons"), slot, "hole")
          output$iso90areaKm[i] <- sum(ifelse(h, -a, a))/1000000
          writeOGR(spdf, dirout, paste(outnames[i], "KUD90", sep=""), "ESRI Shapefile")}
        if(j == 4){
          pls <- slot(spdf, "polygons")[[1]]
          gpclibPermit()
          xx <- checkPolygonsHoles(pls)
          a <- sapply(slot(xx, "Polygons"), slot, "area")
          h <- sapply(slot(xx, "Polygons"), slot, "hole")
          output$iso95areaKm[i] <- sum(ifelse(h, -a, a))/1000000
          writeOGR(spdf, dirout, paste(outnames[i], "KUD95", sep=""), "ESRI Shapefile")}
      }
    }
  }        
  rm(data, data.xy, sppt, locs, HpiOut, fhatOut, fhat.contlev, fhat.contlines, sldf, ps, sp, dataframe, spdf)  
}
write.csv(output, "PluginOutput_All_210628.csv",row.names = FALSE) ### Write output to folder

