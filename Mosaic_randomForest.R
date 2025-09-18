
# creation of a composite climatology (mosaic) for north america, from PRISM and daymet
# use ML to create a synthetic "prism-like" climatology outside the PRISM area. 
# Colin Mahony 
# Jan 1, 2024

library(terra)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "Pr")

# Define the study area (adjust if needed)
studyarea <- ext(c(-150, -110, 54, 72))

# composite DEM for study area
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_akbcus/dem/", sep="")
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))
dem.akbcus <- extend(dem.bc, studyarea) # start with the BC DEM and extend it to the full study area range
dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif") #250m dem downloaded from http://www.cec.org/north-american-environmental-atlas/elevation-2023/
dem.noram <- project(dem.noram, dem.akbcus, method="near") #project 250m source dem to the 800m prism grid. method="near" to preserve elevation variance 
dem.akbcus <- ifel(is.na(dem.akbcus), dem.noram, dem.akbcus)  # replace NA values in the study area raster with the noram DEM value
dem.akbcus <- crop(dem.akbcus, studyarea)
rm(dem.noram)

# create the ocean proximity layer
dem.coarse <- aggregate(dem.akbcus, fact=16)
# plot(dem.coarse)
data(coastsCoarse)
coastsCoarse <- vect(coastsCoarse) #convert to spatvector
# plot(coastsCoarse,add=TRUE,col='blue')
coastal <- distance(dem.coarse, coastsCoarse) 
coastal <- disagg(coastal, fact=16, method="bilinear") 
# plot(coastal)

# use a coastal buffer to create a curved boundary for the daymet component of the predictand. 
coast.pacific <- crop(coastsCoarse, ext(dem.akbcus)-c(0,0,0,10)) #remove the arctic coastline so that we are buffering only the pacific coast
# lines(coast.pacific)
buffer <- buffer(coast.pacific, 1200000)
# lines(buffer)

## ------------------------------------
## loop the variables
e <- 2
m <- 3


for(e in 2:length(elements)){
  for(m in 1:length(monthcodes)){
    
    # load the PRISM mosaic data for the variable
    dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_akbcus/", elements[e], "/", sep="")
    file <- list.files(dir, pattern=paste(".*._", elements[e],"_.*._",monthcodes[m], ".asc", sep=""))
    prism <- rast(paste(dir, file, sep=""))
    prism <- crop(prism, studyarea)
    prism <- project(prism, dem.akbcus) # THIS IS A PATCH. FIX THIS FURTHER UP. DEM SHOULD MATCH PRISM. 
    if(elements[e]=="Pr") values(prism) <- log2(values(prism))
    values(prism)[!is.finite(values(prism))] <- NA
    
    # load the daymet data for the variable
    dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/daymet/"
    daymet <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tmin", "tmax", "pr")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
    daymet <- crop(daymet, prism)
    if(elements[e]=="Pr") values(daymet) <- log2(values(daymet))
    values(daymet)[!is.finite(values(daymet))] <- NA
    
    # load the worldclim data for the variable
    dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/WorldClim/"
    worldclim <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tmin", "tmax", "prec")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
    worldclim <- crop(worldclim, prism)
    if(elements[e]=="Pr") values(worldclim) <- log2(values(worldclim))
    values(worldclim)[!is.finite(values(worldclim))] <- NA
    
    # load the chelsa data for the variable
    dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/CHELSA/"
    chelsa <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tasmin", "tasmax", "pr")[e],"_", monthcodes[m], ".*.", ".tif", sep="")), sep=""))
    chelsa <- crop(chelsa, prism)
    if(elements[e]=="Pr") values(chelsa) <- log2(values(chelsa))
    values(chelsa)[!is.finite(values(chelsa))] <- NA
    
    # use a coastal buffer to create a curved boundary for the daymet component of the predictand. 
    daymet.addarea <- project(daymet, dem.akbcus)
    daymet.addarea <- mask(daymet.addarea, buffer, inverse=T)
    temp <- vect(ext(prism)-c(0,0,12,0)); crs(temp) <- crs(prism)
    daymet.addarea <- mask(daymet.addarea, temp)
    # plot(prism)
    # plot(daymet.addarea, add=T)
    
    # compile prism and daymet into a training raster
    clim.training <- cover(prism, daymet.addarea)
    clim.training <- crop(clim.training, ext(prism)-c(18,0,5,0)) #remove unnecessary area (for training)
    # plot(clim.training)
    
    # extract the training data
    points.prism <- crds(crop(prism, ext(prism)-c(18,0,5,0)), na.rm = T)
    # points.prism <- points.prism[sample(1:dim(points.prism)[1], 2500000),] # about half of the total points
    points.daymet <- crds(daymet.addarea, na.rm = T)
    points.daymet <- points.daymet[sample(1:dim(points.daymet)[1], 100000),] # about 10% of total points (because there is very little topographic relief)
    points <- rbind(points.prism, points.daymet)
    predictand.points <- extract(clim.training,points)
    elev.points <- extract(dem.akbcus,points)
    daymet.points <- extract(daymet,points)
    worldclim.points <- extract(worldclim,points)
    chelsa.points <- extract(chelsa,points)
    coastal.points <- extract(coastal,points)
    data.train <- cbind(predictand.points, points, elev.points, daymet.points, worldclim.points, chelsa.points, coastal.points)
    colnames(data.train) <- c("predictand", "x","y","elev","daymet", "worldclim", "chelsa", "coastal")
    data.train <- data.train[!is.na(data.train$elev),]
    data.train <- data.train[!is.na(data.train$worldclim),]
    data.train <- data.train[!is.na(data.train$daymet),]
    # str(data.train)
    
    # train the model
    model.akbcus <- ranger(y=data.train[,1], x=data.train[,-1])
    
    # points for unmapped area to be filled
    dem.fillarea <- dem.akbcus
    dem.fillarea[!is.na(cover(prism, daymet.addarea))] <- NA
    # plot(dem.fillarea)
    points.fillarea <- crds(dem.fillarea, na.rm = F)
    elev.fillarea <- extract(dem.fillarea,points.fillarea)
    daymet.fillarea <- extract(daymet,points.fillarea)
    worldclim.fillarea <- extract(worldclim,points.fillarea)
    chelsa.fillarea <- extract(chelsa,points.fillarea)
    coastal.fillarea <- extract(coastal,points.fillarea)
    data.fillarea <- cbind(points.fillarea,elev.fillarea,daymet.fillarea, worldclim.fillarea, chelsa.fillarea, coastal.fillarea)
    colnames(data.fillarea) <- c("x","y","elev","daymet","worldclim","chelsa", "coastal")
    data.fillarea$id <- seq_along(data.fillarea$x)
    data.fillarea <- data.fillarea[!is.na(data.fillarea$elev),]
    data.fillarea <- data.fillarea[!is.na(data.fillarea$worldclim),]
    data.fillarea <- data.fillarea[!is.na(data.fillarea$daymet),]
    
    # predict to the unmapped area
    pred.fillarea <- predict(model.akbcus, data = data.fillarea)
    rm(model.akbcus)
    
    # compile the composite raster
    clim.pred.fillarea <- dem.fillarea
    values(clim.pred.fillarea) <- NA
    clim.pred.fillarea[data.fillarea$id] <- pred.fillarea$predictions
    clim.pred <- cover(prism, clim.pred.fillarea)
    clim.pred <- cover(clim.pred, daymet.addarea)
    # plot(clim.pred)
    
    # write the raster
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Projects/2024_intercomparison/outputs/Composite_multiPredictor_climatologies/"
    writeRaster(clim.pred, paste(dir, "composite_WNA_1981_2010_", elements[e], monthcodes[m], ".tif", sep=""), overwrite=T)
    
    print(monthcodes[m])  
  }
  print(elements[e])
}



