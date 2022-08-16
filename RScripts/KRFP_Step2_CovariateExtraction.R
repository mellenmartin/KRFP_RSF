#### loading and extracting raster and landcover data within fisher home ranges ####

rm(list=ls());gc() #clear the memory

#load required packages
library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(stars)
library(rasterVis)
library(nngeo)
library(landscapemetrics)
library(landscapetools)
#requiredPackages <- c( 'stars', 'gridExtra', 'mapview', 
#                      'rasterVis', 'maptools', 'sp', 
#                      'landscapetools') 
#
##### load in and bind BBMM shapefiles #####

library(sf)
library(dplyr)
wd <- "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/FisherHR_KDEPolygons_pluginBandwidth_210923"
setwd(wd)

# load in HR contours
filenames <- list.files(wd, pattern = "*.shp$")

# Create list to store the contours within the loop
KDE.list <- list()

## Import all KDE shapefiles, add individual identifiers, calculate area, and calculate stream density
for(i in 1:length(filenames)){
  # set CRS for the shapefiles
  crs <-  "+proj=utm +zone=11 +ellps=WGS84 +datum=NAD83 units=m"
  
  library(stringr)
  
  temp.sf <- st_read(filenames[i])
  temp.sf <- st_transform(temp.sf, crs = crs)
  names(temp.sf)[names(temp.sf)=="V1"] <- "ID"
  temp.sf$ID <- paste0(filenames[i])
  
  # add identifiers 
  temp.sf$HR_ID <-  substr(temp.sf$ID,1,nchar(temp.sf$ID)-9)
  temp.sf$contour <- str_sub(temp.sf$ID, - 6, - 5)
  temp.sf$Year <-  sub('.*_',"", temp.sf$HR_ID)
  
  # calculate area
  temp.sf$areakm <- as.numeric(st_area(temp.sf)/1000000)
  
  #### calculate stream length and stream density ####
  
  # read in stream shapefile and home range shapefile
  streams <- st_read("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/NHD_Peren_Intermit_WGS84UTM11.shp", quiet=TRUE)
  streams <- st_transform(streams, crs = crs)
  
  # crop streams by each home range
  streamscrop <- st_intersection(streams, temp.sf)
  
  # calculate stream length 
  temp.sf$streamlength <- as.numeric(sum(st_length(streamscrop$geometry)))
  
  # calculate stream density column
  temp.sf$streamdensity <- temp.sf$streamlength/temp.sf$areakm
  
  library(maptools)
  library(raster)
  
  #### extract land cover and summarize ####
  if (temp.sf$Year < 2015){
    bare_2014 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2014/Bare2014.tif")
    forest_2014 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2014/forest2014.tif")
    shrub_2014 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2014/Shrub2014.tif")
    
    bare_2014 <- projectRaster(bare_2014, crs=crs(temp.sf))
    shrub_2014 <- projectRaster(shrub_2014, crs=crs(temp.sf))
    forest_2014 <- projectRaster(forest_2014, crs=crs(temp.sf))
    # bare_summary
    bare <- raster::extract(bare_2014, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    forest <- raster::extract(forest_2014, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    shrub <- raster::extract(shrub_2014, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)

    temp.sf$bare <- bare$Bare2014
    temp.sf$forest <- forest$forest2014
    temp.sf$shrub <- shrub$Shrub2014
    temp.sf$mort <- 0
    temp.sf$total <- temp.sf$bare+temp.sf$shrub+temp.sf$forest
    temp.sf$barepercent <- temp.sf$bare/temp.sf$total*100
    temp.sf$shrubpercent <- temp.sf$shrub/temp.sf$total*100
    temp.sf$forestpercent <- temp.sf$forest/temp.sf$total*100
    temp.sf$mortpercent <- 0
    
    #### calculate forest connectivity in the HR 
    forest_bin <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2014/DenseForestBinary_2014.tif")
    forest_masked <- raster::mask(forest_bin, temp.sf)
    forest_trimmed <- trim(forest_masked, values = NA)
    
    #turn into landscape
    cohesion <- lsm_c_cohesion(forest_trimmed, directions = 8)
    cohesion <- cohesion[2, , drop = TRUE]
    temp.sf$cohesion <- cohesion$value

    
  } else if (temp.sf$Year == 2015){
    bare_2015 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2015/Bare2015.tif")
    forest_2015 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2015/forest2015.tif")
    shrub_2015 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2015/Shrub2015.tif")
    mort_2015 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2015/TreeMortality2015.tif")
    
    bare_2015 <- projectRaster(bare_2015, crs=crs(temp.sf))
    shrub_2015 <- projectRaster(shrub_2015, crs=crs(temp.sf))
    forest_2015 <- projectRaster(forest_2015, crs=crs(temp.sf))
    mort_2015 <- projectRaster(mort_2015, crs=crs(temp.sf))
    
    bare <- raster::extract(bare_2015, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    forest <- raster::extract(forest_2015, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    shrub <- raster::extract(shrub_2015, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    mort <- extract(mort_2015, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    temp.sf$bare <- bare$Bare2015
    temp.sf$forest <- forest$forest2015
    temp.sf$shrub <- shrub$Shrub2015
    temp.sf$mort <- mort$TreeMortality2015
    temp.sf$total <- temp.sf$bare+temp.sf$shrub+temp.sf$forest+temp.sf$mort 
    temp.sf$barepercent <- temp.sf$bare/temp.sf$total*100
    temp.sf$shrubpercent <- temp.sf$shrub/temp.sf$total*100
    temp.sf$forestpercent <- temp.sf$forest/temp.sf$total*100
    temp.sf$mortpercent <- temp.sf$mort/temp.sf$total*100
    
    #### calculate forest connectivity in the HR 
    forest_bin <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2015/DenseForestBinary_2015.tif")
    forest_masked <- raster::mask(forest_bin, temp.sf)
    forest_trimmed <- trim(forest_masked, values = NA)
    
    #turn into landscape
    cohesion <- lsm_c_cohesion(forest_trimmed, directions = 8)
    cohesion <- cohesion[2, , drop = TRUE]
    temp.sf$cohesion <- cohesion$value
    
  } else if (temp.sf$Year == 2016){
    bare_2016 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2016/Bare2016.tif")
    forest_2016 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2016/forest2016.tif")
    shrub_2016 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2016/Shrub2016.tif")
    mort_2016 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2016/TreeMortality2016.tif")
    
    bare_2016 <- projectRaster(bare_2016, crs=crs(temp.sf))
    shrub_2016 <- projectRaster(shrub_2016, crs=crs(temp.sf))
    forest_2016 <- projectRaster(forest_2016, crs=crs(temp.sf))
    mort_2016 <- projectRaster(mort_2016, crs=crs(temp.sf))
    
    bare <- raster::extract(bare_2016, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    forest <- raster::extract(forest_2016, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    shrub <- raster::extract(shrub_2016, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    mort <- extract(mort_2016, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    temp.sf$bare <- bare$Bare2016
    temp.sf$forest <- forest$forest2016
    temp.sf$shrub <- shrub$Shrub2016
    temp.sf$mort <- mort$TreeMortality2016
    temp.sf$total <- temp.sf$bare+temp.sf$shrub+temp.sf$forest+temp.sf$mort 
    temp.sf$barepercent <- temp.sf$bare/temp.sf$total*100
    temp.sf$shrubpercent <- temp.sf$shrub/temp.sf$total*100
    temp.sf$forestpercent <- temp.sf$forest/temp.sf$total*100
    temp.sf$mortpercent <- temp.sf$mort/temp.sf$total*100
    
    #### calculate forest connectivity in the HR 
    forest_bin <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2016/DenseForestBinary_2016.tif")
    forest_masked <- raster::mask(forest_bin, temp.sf)
    forest_trimmed <- trim(forest_masked, values = NA)
    
    #turn into landscape
    cohesion <- lsm_c_cohesion(forest_trimmed, directions = 8)
    cohesion <- cohesion[2, , drop = TRUE]
    temp.sf$cohesion <- cohesion$value
    
  } else if (temp.sf$Year == 2017){
    bare_2017 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2017/Bare2017.tif")
    forest_2017 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2017/forest2017.tif")
    shrub_2017 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2017/Shrub2017.tif")
    mort_2017 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2017/TreeMortality2017.tif")
    
    bare_2017 <- projectRaster(bare_2017, crs=crs(temp.sf))
    shrub_2017 <- projectRaster(shrub_2017, crs=crs(temp.sf))
    forest_2017 <- projectRaster(forest_2017, crs=crs(temp.sf))
    mort_2017 <- projectRaster(mort_2017, crs=crs(temp.sf))
    
    bare <- raster::extract(bare_2017, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    forest <- raster::extract(forest_2017, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    shrub <- raster::extract(shrub_2017, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    mort <- extract(mort_2017, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    temp.sf$bare <- bare$Bare2017
    temp.sf$forest <- forest$forest2017
    temp.sf$shrub <- shrub$Shrub2017
    temp.sf$mort <- mort$TreeMortality2017
    temp.sf$total <- temp.sf$bare+temp.sf$shrub+temp.sf$forest+temp.sf$mort 
    temp.sf$barepercent <- temp.sf$bare/temp.sf$total*100
    temp.sf$shrubpercent <- temp.sf$shrub/temp.sf$total*100
    temp.sf$forestpercent <- temp.sf$forest/temp.sf$total*100
    temp.sf$mortpercent <- temp.sf$mort/temp.sf$total*100
    
    #### calculate forest connectivity in the HR 
    forest_bin <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2017/DenseForestBinary_2017.tif")
    forest_masked <- raster::mask(forest_bin, temp.sf)
    forest_trimmed <- trim(forest_masked, values = NA)
    
    #turn into landscape
    cohesion <- lsm_c_cohesion(forest_trimmed, directions = 8)
    cohesion <- cohesion[2, , drop = TRUE]
    temp.sf$cohesion <- cohesion$value
    
  } else {
    bare_2018 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2018/Bare2018.tif")
    forest_2018 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2018/forest2018.tif")
    shrub_2018 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2018/Shrub2018.tif")
    mort_2018 <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2018/TreeMortality2018.tif")
    
    bare_2018 <- projectRaster(bare_2018, crs=crs(temp.sf))
    shrub_2018 <- projectRaster(shrub_2018, crs=crs(temp.sf))
    forest_2018 <- projectRaster(forest_2018, crs=crs(temp.sf))
    mort_2018 <- projectRaster(mort_2018, crs=crs(temp.sf))
    
    bare <- raster::extract(bare_2018, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    forest <- raster::extract(forest_2018, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    shrub <- raster::extract(shrub_2018, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    mort <- extract(mort_2018, temp.sf, fun=sum, na.rm=TRUE, df=TRUE)
    temp.sf$bare <- bare$Bare2018
    temp.sf$forest <- forest$forest2018
    temp.sf$shrub <- shrub$Shrub2018
    temp.sf$mort <- mort$TreeMortality2018
    temp.sf$total <- temp.sf$bare+temp.sf$shrub+temp.sf$forest+temp.sf$mort 
    temp.sf$barepercent <- temp.sf$bare/temp.sf$total*100
    temp.sf$shrubpercent <- temp.sf$shrub/temp.sf$total*100
    temp.sf$forestpercent <- temp.sf$forest/temp.sf$total*100
    temp.sf$mortpercent <- temp.sf$mort/temp.sf$total*100
    
    #### calculate forest connectivity in the HR 
    forest_bin <- raster("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/2018/DenseForestBinary_2018.tif")
    forest_masked <- raster::mask(forest_bin, temp.sf)
    forest_trimmed <- trim(forest_masked, values = NA)
    
    #turn into landscape
    cohesion <- lsm_c_cohesion(forest_trimmed, directions = 8)
    cohesion <- cohesion[2, , drop = TRUE]
    temp.sf$cohesion <- cohesion$value
    
    
  }
  
  
  #### add to list ####
  temp <- list(temp.sf)
  # add to list
  KDE.list[[i]] <- temp.sf
}

combineShp <- do.call(what = sf:::rbind.sf, args = KDE.list)
combineShp <- st_transform(combineShp, crs = crs)

#combineShp$Year[combineShp$Year=="_201"]<-2019

# write the combined HR polygons w/ stream info for future use
st_write(combineShp, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/HRShapefiles/CombinedFishers_KUDAll_HRCovs_210925_2.shp")

##### combine polygons for overall estimates ####
library(sf)
library(dplyr)
crs <-  "+proj=utm +zone=11 +ellps=WGS84 +datum=NAD83 units=m"

combineShp <- st_read("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/HRShapefiles/CombinedFishers_KUDAll_HRCovs_210925.shp")
combineShp <- st_transform(combineShp, crs = crs)

# combine polygons buffer polygons by 15km, ~1/2 male HR
combineShp_dissolve <- st_combine(combineShp)
combineShp_buffer <- st_buffer(combineShp_dissolve, dist = 2000)
combineShp_RSF_buffer <- st_buffer(combineShp, dist = 2000)
combineShp_RSF_buffer <- combineShp_RSF_buffer[combineShp_RSF_buffer$contour == 95,]
st_write(combineShp_buffer, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/HRShapefiles/KRFP_FisherHRs_DissolvedBuffered_210924.shp")
st_write(combineShp_RSF_buffer, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/HRShapefiles/KRFP_FisherHRs_RSF_Buffered_210923.shp")


##### extract covariates at used locations for RSF analysis #####
##### load packages 
library(rgdal)
library(sp)
library(raster)    
library(sdm)
library(dplyr)
library(sfheaders)
setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF") 

# read in the filtered locations we used for the HR analyses
useddata <- read.csv("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/LocationData/KRFP_FilteredHRLocations_210923.csv")
polyIDs <- unique(combineShp_RSF_buffer$HR_ID)
useddataIDs <- unique(useddata$ID_year)

#useddata<-useddata[!(useddata$ID_year=="F08_2009"),]

# select relevant columns
useddata <- useddata %>% dplyr::select(LocID, 
                                       Type, 
                                       FisherID, 
                                       Sex, 
                                       AgeTransition, 
                                       Date,
                                       Period,
                                       Time24,
                                       FisherYear,
                                       x_,
                                       y_,
                                       ID_year,
                                       datetime_format,
                                       julian)
useddata$Use <- rep(1, length(useddata$LocID))
useddata <- tibble::rowid_to_column(useddata, "ID")

usedsum <- useddata %>% count(ID_year)

## sample available points by home range ##
library(rgdal)
library(sp)

useddata$ID_year <- factor(useddata$ID_year)


polylist <- list()
loclist <- list()
for (i in 1:length(polyIDs)){ 
  
  # id for polygon "i" in your list of polygos
  id <- polyIDs[i]
  
  # define the coordinate system 
  crs <-  "+proj=utm +zone=11 +ellps=WGS84 +datum=NAD83 units=m"
  
  # subset to polygon that matches od
  used <- useddata[useddata$ID_year == id,]
  #used <- na.omit(used)
  usedIDs <- used$ID
  
  for (j in 1:length(usedIDs)){
    
    ptID <- as.numeric(usedIDs[j])
    
    pt <- used[used$ID == ptID,]
    
    # create spatialpoints layer to generate 10 random locs per loc
    locs <- SpatialPoints(pt[, c( "x_", "y_" )])
    
    proj4string(locs) <- CRS("+proj=utm +zone=11 +datum=NAD83")
    
    poly <- combineShp_RSF_buffer[combineShp_RSF_buffer$HR_ID == id,]
    
    poly <- st_transform(poly, crs = crs)
    
    #b <- buffer(loc, width = 402)
    
    #locs <- st_sample(poly, poly$locs*mult)
    locs <-st_sample(poly,10,type="random")
    
    coords <- st_coordinates(locs)
    avail_locs <- as.data.frame(coords, row.names=FALSE)
    
    names(avail_locs)[1] <- "x_"
    names(avail_locs)[2] <- "y_"
    avail_locs$LocID <- rep(pt$ID)
    avail_locs$Type <- rep("RAND")
    avail_locs$FisherID <- rep(pt$FisherID)  
    avail_locs$Sex <- rep(pt$Sex)
    avail_locs$AgeTransition <- rep(pt$AgeTransition)
    avail_locs$Date <- rep(pt$Date)
    avail_locs$Period <- rep(pt$Period)
    avail_locs$Time24 <- rep(pt$Time24)
    avail_locs$FisherYear <- rep(pt$FisherYear)
    avail_locs$ID_year <- rep(pt$ID_year)
    avail_locs$datetime_format <- rep(pt$datetime_format)
    avail_locs$julian <- rep(pt$julian)
    avail_locs$Use <- rep(0)
    
    
    # add the subsetted dataframe that now has snow depth to the dataframe list you created
    loclist[[j]] <- avail_locs  
  }
  
  polyavail_locs <- do.call(rbind, loclist)
  polylist[[i]] <- polyavail_locs
}

allavail_locs <- do.call(rbind, polylist)
allavail_locs$ID <- allavail_locs$LocID

allavail_locs_sub <- allavail_locs %>% group_by(LocID) %>% sample_n(10)

write.csv(allavail_locs_sub, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_AvailableLocations_210925_2.csv", row.names = FALSE)

## now bind those available locations to the used locations
alldata <- bind_rows(useddata, allavail_locs_sub)

write.csv(alldata, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_UsedandAvailableLocs_210925_2.csv", row.names = FALSE)

##### read in all locations and rasters to extract data ####
alldata <- read.csv("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_UsedandAvailableLocs_210924.csv")

library(rgdal)
library(sp)
library(raster)    

# create identifier for extraction loop
years <- unique(alldata$FisherYear)

extent <- extent(min(alldata$x_), max(alldata$x_), min(alldata$y_), max(alldata$y_))

start <- Sys.time()

# create empty list to store data
dflist <- list()

for (i in 1:length(years)){ 
  
  # define the coordinate system 
  crs <-  "+proj=utm +zone=11 +ellps=WGS84 +datum=NAD83 units=m"
  
  # id for raster "i" in your list of rasters
  id <- as.numeric(years[i])

  subdata <- alldata[alldata$FisherYear == id,]

  # can also add "pattern = "*.tif" if not all of the files in the folder are rasters
  rastlist <-list.files(paste0("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/RSF_",id))
  
  # can also add "pattern = "*.tif" if not all of the files in the folder are rasters
  # rastlist<-list.files(wd, pattern = "*.tif")

  covlist <- list()
  
  for (j in 1:length(rastlist)){ 
    
    # id for raster "i" in your list of rasters
    rastid <- rastlist[j]
    
    # now, bring in snow depth raster that matches "j" in rastlist
    raster <-(paste0("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/RSF_",id,"/",rastid))
    
    raster <- raster(raster)
    
    # project the raster with that coordinate system
    raster <- projectRaster(raster, crs = crs(combineShp))
    
    library(stringr)
    rastername <- sub("^(\\D+).*", "\\1", rastid)
    
    
    # create spatialpoints layer to extract covariate at locations in subset of observations (i.e., "GPS")
    pts <- SpatialPoints( subdata[, c( "x_", "y_" ) ], 
                          proj4string = CRS(proj4string(raster)))
    
    # extract values at all locations (locs is actually both used and random points)
    raster_cov <- raster::extract(raster, pts)
    
    raster_data <- as.data.frame(raster_cov, row.names = NULL)
  
    # rename column
    names(raster_data)[names(raster_data) == "raster_cov"] <- paste0(rastername)
    
    # add the subsetted dataframe to the dataframe list you created
    covlist[[j]] <- raster_data  
    
    gc()
  }
  
  # bind all dataframes in the list
  allcov_data <- do.call(cbind, covlist)
  allcov_data2 <- cbind(subdata,allcov_data)
  
  # add the subsetted dataframe that now has snow depth to the dataframe list you created
  dflist[[i]] <- allcov_data2  
}

end <- Sys.time()

# bind all dataframes in the list
library(data.table)
alldata_covs <- rbindlist(dflist, fill = TRUE)

# save for later
write.csv(alldata_covs, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_UsedAvail_Covariates_210925_2.csv", row.names = FALSE)

##### read in used/avail and compute summary stats and figures #####
library(dplyr)
library(ggplot2)
library(GGally)
library(corrplot)
library(tidyr)

alldata <- read.csv("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_UsedAvail_Covariates_210813.csv")

alldata <- subset(alldata, select = -c(21,33) )
alldata$Sex <- as.character(alldata$Sex)
alldata$Sex[alldata$Sex == "f"] <- "F"


randlocs <- alldata %>% dplyr::filter(Type == "RAND")
usedlocs <- alldata %>% dplyr::filter(Type != "RAND")



# corr plot
df.corr <- dplyr::select(usedlocs, 17:34)
df.corr.out <- ggcorr(df.corr, method = c("spearman"), label = TRUE, hjust = 0.75) # option 1
df.corr.out

# cor table
corr_table <- cor(df.corr, method = c("spearman"), use = "complete.obs")
corr_table <- round(corr_table,2)

write.csv(corr_table, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_RSFCorrTable_210708.csv", row.names = TRUE)

##### home range summaries #####
#### load in summaries for HR percentage summary stats ####
library(sf)
library(dplyr)
crs <-  "+proj=utm +zone=11 +ellps=WGS84 +datum=NAD83 units=m"

combineShp <- st_read("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/Shapefiles/HRShapefiles/CombinedFishers_KUDAll_HRCovs_210925.shp")
combineShp <- st_transform(combineShp, crs = crs)

combineShp$Sex<-  substr(combineShp$HR_ID,1,1)

st_geometry(combineShp) <- NULL
class(combineShp)
## [1] "data.frame"

# fix the "period" column
combineShp$Period <- ifelse(combineShp$Year < 2012, "Pre-Drought", (ifelse(combineShp$Year > 2012 & combineShp$Year < 2015, "Drought", "Tree Mortality")))

write.csv(combineShp, "AnnualHomeRange_IndividualHRSizeLandcover_210926.csv", row.names = FALSE)

indsums <- read.csv("AnnualHomeRange_IndividualHRSizeLandcover_210926.csv")


options(digits=2)
annualsummaries <- indsums %>%
  group_by(Sex, contour, Period) %>%
  summarize(Individuals = n_distinct(HR_ID),
            forest_median = median(forestpercent),
            forest_min = min(forestpercent),
            forest_max = max(forestpercent),
            forest_sd = sd(forestpercent),
            mort_median = median(mortpercent),
            mort_min = min(mortpercent),
            mort_max = max(mortpercent),
            mort_sd = sd(mortpercent),
            shrub_median = median(shrubpercent),
            shrub_min = min(shrubpercent),
            shrub_max = max(shrubpercent),
            shrub_sd = sd(shrubpercent),
            bare_median = median(barepercent),
            bare_min = min(barepercent),
            bare_max = max(barepercent),
            bare_sd = sd(barepercent),
            area_mean = mean(areakm),
            area_median = median(areakm),
            area_min = min(areakm),
            area_max = max(areakm),
            area_sd = sd(areakm))

write.csv(annualsummaries, "AnnualHomeRange_HRSizeLandcoverSummaries_210926.csv", row.names = FALSE)

# read back in
annsums <- read.csv("AnnualHomeRange_HRSizeLandcoverSummaries_210926.csv")
annsums$SexFull <- ifelse(annsums$Sex == "F", "Female", "Male")

#re-order for plotting
##### bar plots and summary plots ####
library(ggplot2)

annsums2 <- annsums[,-c(5:7)]
annsums2 <- annsums2[,-c(1)]

# annual summaries
options(digits=2)
hrsizesum <- annsums2 %>%
  group_by(SexFull, Period, contour) %>%
  summarize(Mean = mean(AreaMean),
            Median = median(AreaMedian),
            Minimum = min(AreaMin),
            Maximum = max(AreaMax),
            SD = max(AreaSD)
  )

write.csv(hrsizesum, "AnnualHomeRange_SizeSummaries_210926.csv", row.names = FALSE)


# re order the drought period for graphing
annsums$Period <- factor(annsums$Period, levels = c("Pre-Drought", "Drought", "Tree Mortality"), ordered = TRUE)


png(filename = "KRFP_HRComposition_210926.png", width = 17, height = 12, units = "in", res = 300)
ggplot(annsums, aes(x = Period, y = Median, fill = CoverType)) +
  geom_col(position = "fill") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  labs(x="Period", y = "Proportion of Home Range Contour")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22),
        strip.text.x=element_text(size=18),
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  facet_wrap(~ SexFull + contour)
dev.off()

ggsave("SexSpecificHRComp_2.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       dpi = 320)



#### home range summary bar plot #####
hrsums <- read.csv("AnnualHomeRange_AllSummaries_210527.csv")
hrsums$Period <- factor(hrsums$Period, levels = c("Pre-Mortality", "Drought", "Mortality"), ordered = TRUE)

ggplot(data=hrsums, aes(x=Year, y=Median)) +
  geom_bar(stat="identity", position = position_dodge())+
  labs(x="Fisher Year", y = bquote("Home Range Contour Area"~km^2))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22),
        strip.text.x=element_text(size=18),
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  facet_wrap(~SexFull + Contour)

ggsave("SexSpecificHRSize_1.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       dpi = 320)


#### summary for tables #####
library(dplyr)
options(digits=2)
usedsummaries <- usedlocs %>%
  summarize(bare = median(Bare, na.rm = TRUE),
            ddi = median(ddi, na.rm = TRUE),
            denseforestbin = median(DenseForestBinary, na.rm = TRUE),
            distdenseforest = median(DistDenseForest, na.rm = TRUE),
            fortypba = median(fortypba, na.rm = TRUE),
            liveforest = median(LiveForest, na.rm = TRUE),
            LST = median(LST_SummerMean, na.rm = TRUE),
            qmddom = median(qmd_dom, na.rm = TRUE),
            qukeba = median(quke_ba, na.rm = TRUE),
            ogsi = median(ogsi, na.rm = TRUE),
            shrub = median(Shrub, na.rm = TRUE),
            streamdist = median(StreamDist, na.rm = TRUE),
            tpi = median(TPI, na.rm = TRUE),
            mortforest = median(TreeMortality, na.rm = TRUE))#,
#bare_sd = sd(Bare, na.rm = TRUE),
#ddi_sd = sd(ddi, na.rm = TRUE),
#denseforestbin_sd = sd(DenseForestBinary, na.rm = TRUE),
#distdenseforest_sd = sd(DistDenseForest, na.rm = TRUE),
#fortypba_sd = sd(fortypba, na.rm = TRUE),
#liveforest_sd = sd(LiveForest, na.rm = TRUE),
#LST_sd = sd(LST_SummerMean, na.rm = TRUE),
#qmddom_sd = sd(qmd_dom, na.rm = TRUE),
#qukeba_sd = sd(quke_ba, na.rm = TRUE),
#sdi_sd = sd(sdi_reineke, na.rm = TRUE),
#shrub_sd = sd(Shrub, na.rm = TRUE),
#streamdist_sd = sd(StreamDist, na.rm = TRUE),
#tpi_sd = sd(TPI, na.rm = TRUE),
#mortforest_sd = sd(TreeMortality, na.rm = TRUE))
round(usedsummaries[c(2:14)],2)            
View(usedsummaries)
usedsummaries2 <- usedsummaries %>% pivot_longer(cols= bare:mortforest,
                                                 names_to="Covariate", 
                                                 values_to = "Used")


usedsdsummaries <- usedlocs %>%
  summarize(bare = sd(Bare, na.rm = TRUE),
            ddi = sd(ddi, na.rm = TRUE),
            denseforestbin = sd(DenseForestBinary, na.rm = TRUE),
            distdenseforest = sd(DistDenseForest, na.rm = TRUE),
            fortypba = sd(fortypba, na.rm = TRUE),
            liveforest = sd(LiveForest, na.rm = TRUE),
            LST = sd(LST_SummerMean, na.rm = TRUE),
            qmddom = sd(qmd_dom, na.rm = TRUE),
            qukeba = sd(quke_ba, na.rm = TRUE),
            sdi = sd(sdi_reineke, na.rm = TRUE),
            shrub = sd(Shrub, na.rm = TRUE),
            streamdist = sd(StreamDist, na.rm = TRUE),
            tpi = sd(TPI, na.rm = TRUE),
            mortforest = sd(TreeMortality, na.rm = TRUE))#,

usedsummariessd <- usedsdsummaries %>% pivot_longer(cols= bare:mortforest,
                                                    names_to="Covariate", 
                                                    values_to = "Used_sd")

usedsummariessd <- inner_join(usedsummaries2, usedsummariessd, by = "Covariate")


randsummaries <- randlocs %>%
  summarize(bare = median(Bare, na.rm = TRUE),
            ddi = median(ddi, na.rm = TRUE),
            denseforestbin = median(DenseForestBinary, na.rm = TRUE),
            distdenseforest = median(DistDenseForest, na.rm = TRUE),
            fortypba = median(fortypba, na.rm = TRUE),
            liveforest = median(LiveForest, na.rm = TRUE),
            LST = median(LST_SummerMean, na.rm = TRUE),
            qmddom = median(qmd_dom, na.rm = TRUE),
            qukeba = median(quke_ba, na.rm = TRUE),
            sdi = median(sdi_reineke, na.rm = TRUE),
            shrub = median(Shrub, na.rm = TRUE),
            streamdist = median(StreamDist, na.rm = TRUE),
            tpi = median(TPI, na.rm = TRUE),
            mortforest = median(TreeMortality, na.rm = TRUE))#,
#bare_sd = sd(Bare, na.rm = TRUE),
#ddi_sd = sd(ddi, na.rm = TRUE),
#denseforestbin_sd = sd(DenseForestBinary, na.rm = TRUE),
#distdenseforest_sd = sd(DistDenseForest, na.rm = TRUE),
#fortypba_sd = sd(fortypba, na.rm = TRUE),
#liveforest_sd = sd(LiveForest, na.rm = TRUE),
#LST_sd = sd(LST_SummerMean, na.rm = TRUE),
#qmddom_sd = sd(qmd_dom, na.rm = TRUE),
#qukeba_sd = sd(quke_ba, na.rm = TRUE),
#sdi_sd = sd(sdi_reineke, na.rm = TRUE),
#shrub_sd = sd(Shrub, na.rm = TRUE),
#streamdist_sd = sd(StreamDist, na.rm = TRUE),
#tpi_sd = sd(TPI, na.rm = TRUE),
#mortforest_sd = sd(TreeMortality, na.rm = TRUE))
round(randsummaries[c(2:14)],2)            
View(randsummaries)
randsummaries2 <- randsummaries %>% pivot_longer(cols= bare:mortforest,
                                                 names_to="Covariate", 
                                                 values_to = "Rand")


randsdsummaries <- randlocs %>%
  summarize(bare = sd(Bare, na.rm = TRUE),
            ddi = sd(ddi, na.rm = TRUE),
            denseforestbin = sd(DenseForestBinary, na.rm = TRUE),
            distdenseforest = sd(DistDenseForest, na.rm = TRUE),
            fortypba = sd(fortypba, na.rm = TRUE),
            liveforest = sd(LiveForest, na.rm = TRUE),
            LST = sd(LST_SummerMean, na.rm = TRUE),
            qmddom = sd(qmd_dom, na.rm = TRUE),
            qukeba = sd(quke_ba, na.rm = TRUE),
            sdi = sd(sdi_reineke, na.rm = TRUE),
            shrub = sd(Shrub, na.rm = TRUE),
            streamdist = sd(StreamDist, na.rm = TRUE),
            tpi = sd(TPI, na.rm = TRUE),
            mortforest = sd(TreeMortality, na.rm = TRUE))#,

randsummariessd <- randsdsummaries %>% pivot_longer(cols= bare:mortforest,
                                                 names_to="Covariate", 
                                                 values_to = "Rand_sd")

randsummariessd <- inner_join(randsummaries2, randsummariessd, by = "Covariate")
usedrandsummaries <- inner_join(usedsummariessd, randsummariessd, by = "Covariate") 

write.csv(usedrandsummaries, "C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/RandUsed_RSFCovSummaries_210712.csv", row.names = FALSE)


#re-order for plotting
##### bar plots and summary plots ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(viridis)
usedlong <- usedlocs[c(17:32)] %>%
  pivot_longer(everything())
names(usedlong) <- c( "Covariate","Value")
usedlong$Type <- rep( "Used", length(usedlong$Covariate))

randlong <- randlocs[c(17:32)] %>%
  pivot_longer(everything())
names(randlong) <- c( "Covariate","Value")
randlong$Type <- rep( "Random", length(randlong$Covariate))


usedavail <- rbind(usedlong,randlong)

bare <- filter(usedavail, Covariate == "Bare")
ddi <-  filter(usedavail, Covariate == "ddi")                 
denseforestbin <- filter(usedavail, Covariate == "DenseForestBinary")
distdenseforest <-  filter(usedavail, Covariate == "DistDenseForest")                 
fortypba <- filter(usedavail, Covariate == "fortypba")
landcoverclass <-  filter(usedavail, Covariate == "LandCoverClassified")                 
landcover5m <- filter(usedavail, Covariate == "LandCovFivem")
liveforest <-  filter(usedavail, Covariate == "LiveForest")                 
lst <- filter(usedavail, Covariate == "LST_SummerMean")
qmddom <-  filter(usedavail, Covariate == "qmd_dom")                 
qukeba <- filter(usedavail, Covariate == "quke_ba")
sdi <-  filter(usedavail, Covariate == "sdi_reineke")                 
shrub <- filter(usedavail, Covariate == "Shrub")
streamdist <-  filter(usedavail, Covariate == "StreamDist")                 
tpi <- filter(usedavail, Covariate == "TPI")
treemort <-  filter(usedavail, Covariate == "TreeMortality")                 


bare_plot <- ggplot(bare, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-10, 105)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Bare", y = "")

ddi_plot <- ggplot(ddi, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-5,15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="DDI", y = "")

denseforestbin_plot <- ggplot(denseforestbin, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-5,5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="DForestBin", y = "")

distdenseforest_plot <- ggplot(distdenseforest, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-100,900) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="DistDForest", y = "")


fortypba_plot <- ggplot(fortypba, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-5,1000) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="40 PBA", y = "")

liveforest_plot <- ggplot(liveforest, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-10,110) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Forest", y = "")

lst_plot <- ggplot(lst, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(280,330) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="LST", y = "")

qmddom_plot <- ggplot(qmddom, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-15,120) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="QMD Dom.", y = "")

qukeba_plot <- ggplot(qukeba, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-5,25) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Quke BA", y = "")

sdi_plot <- ggplot(sdi, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-100,800) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="SDI", y = "")

shrub_plot <- ggplot(shrub, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-15,100) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Shrub", y = "")

streamdist_plot <- ggplot(streamdist, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(-60,600) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Stream Dist.", y = "")

treemort_plot <- ggplot(treemort, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram() +
  xlim(5,110) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="TreeMort", y = "")



library(ggpubr)
cov_plot <- ggarrange(bare_plot, ddi_plot, denseforestbin_plot, distdenseforest_plot, fortypba_plot, 
                      liveforest_plot, lst_plot, qmddom_plot, qukeba_plot, sdi_plot, shrub_plot, streamdist_plot,
                      tpi_plot, treemort_plot)
ggsave("KRFP_RSFCovsSum_210709.jpg", plot = last_plot(), dpi=300)



# re order the drought period for graphing
annsums$Period <- factor(annsums$Period, levels = c("Pre-Mortality", "Drought", "Mortality"), ordered = TRUE)

barplot <- ggplot(annsums, aes(x = Period, y = median, fill = CoverType)) +
  geom_col(position = "fill") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  labs(x="Period", y = "Proportion of Home Range Contour")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22),
        strip.text.x=element_text(size=18),
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  facet_wrap(~ SexFull + Contour)

ggsave("SexSpecificHRComp_2.png",
       plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       dpi = 320)

