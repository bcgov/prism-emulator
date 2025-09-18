
# compare climatological maps from various source using the BC PRISM DEM domain
# Colin Mahony 

library(terra)
library(climr)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

# function for preparing data
prep <- function(x, studyarea, element, breaks){
  x <- crop(x, studyarea)
  if(element=="Pr") values(x) <- log2(values(x))
  values(x)[!is.finite(values(x))] <- NA
  values(x)[values(x)>max(breaks)] <- max(breaks)
  values(x)[values(x)<min(breaks)] <- min(breaks)
  return(x)
}

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "Pr")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

e <- 2
m <- 7

# load the source STATION data for the BC prism
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
stn.info <- stn.info[-which(El_Flag=="@"),]
stn.data <- stn.info[,get(month.abb[m])]
stn.data <- if(elements[e]=="Pr") log2(stn.data) else stn.data/10
stn.info <- stn.info[is.finite(stn.data),]
stn.data <- stn.data[is.finite(stn.data)]

#PRISM DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))

# load the BC PRISM  data for the variable
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/", sep="")
file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
prism.bc <- rast(paste(dir, file, sep=""))
if(elements[e]=="Pr") values(prism.bc) <- log2(values(prism.bc))
values(prism.bc)[!is.finite(values(prism.bc))] <- NA

# color scheme
combined <- c(values(prism.bc), stn.data)
combined <- combined[is.finite(combined)]
inc=diff(range(combined))/500
breaks=seq(quantile(combined, 0.001)-inc, quantile(combined, .999)+inc, inc)
ColScheme <- colorRampPalette(if(elements[e]=="Pr") brewer.pal(9, "YlGnBu") else rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColPal <- colorBin(ColScheme, bins=breaks, na.color = "white")
ColPal.raster <- colorBin(ColScheme, bins=breaks, na.color = "transparent")

# load the western canada 2km PRISM data for the variable
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_canw/PRISM_",c("tmin", "tmax", "ppt")[e],"_canw_1961-1990_normal_2kmM1_", monthcodes[m], "_asc/", sep="")
file <- list.files(dir)
prism.canw <- rast(paste(dir, file, sep=""))/100
prism.canw <- prep(prism.canw, studyarea=dem.bc, element=elements[e], breaks=breaks)
delta.e <- rast(paste0("//objectstore2.nrs.bcgov/ffec/TransferAnomalies/delta.from.1961_1990.to.1981_2010.", elements[e], ".tif"))
delta <- project(delta.e[[m]], prism.canw)
prism.canw <- if(elements[e]=="Pr") prism.canw/delta else prism.canw+delta

# load the worldclim data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/WorldClim/"
worldclim <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tmin", "tmax", "prec")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
worldclim <- prep(worldclim, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the CHELSA data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/CHELSA/"
chelsa <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tasmin", "tasmax", "pr")[e],"_", monthcodes[m], ".*.", ".tif", sep="")), sep=""))
chelsa <- prep(chelsa, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the daymet data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/daymet/"
daymet <- rast(paste(dir, list.files(dir, pattern=paste(".*.", c("tmin", "tmax", "pr")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
daymet <- prep(daymet, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the ClimateNA data for the variable
my_grid <- as.data.frame(dem.bc, cells = TRUE, xy = TRUE)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects
var <- paste(c("Tmin", "Tmax", "PPT")[e], monthcodes[m], sep="_")
climr <- downscale(xyz = my_grid, which_refmap = "refmap_climatena", vars = var)
climatena <- rast(dem.bc) # use the DEM as a template raster
climatena[climr[, id]] <- climr[,get(var)]
climatena <- prep(climatena, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the ANUSPLIN data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/ANUSPLIN/mly60arcsecond_1981-2010/"
anusplin <- rast(paste(dir, list.files(dir, pattern=paste(c("mint60", "maxt60", "pcp60")[e],"_", monthcodes[m], ".tif", sep="")), sep=""))
anusplin <- prep(anusplin, studyarea=dem.bc, element=elements[e], breaks=breaks)

# load the USask WRF data for the variable
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/"
dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
wrfUsask <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
wrfUsask <- prep(wrfUsask, studyarea=dem.bc, element=elements[e], breaks=breaks)


# leaflet map
labels <- paste(stn.info$Name, "(El. ", stn.info$Elevation, "m)", sep="")
map <- leaflet(stn.info) %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  # addRasterImage(dem, colors =terrain.colors(99), opacity = 1, maxBytes = 6 * 1024 * 1024, group = "elevation") %>%
  # addRasterImage(dem.usask, colors =terrain.colors(99), opacity = 1, maxBytes = 6 * 1024 * 1024, group = "elevation") %>%
  addRasterImage(prism.bc, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "BC PRISM") %>%
  addRasterImage(prism.canw, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "W. Can. PRISM") %>%
  addRasterImage(daymet, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "Daymet") %>%
  addRasterImage(chelsa, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "CHELSA") %>%
  addRasterImage(worldclim, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WorldClim (1971-2000)") %>%
  addRasterImage(anusplin, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "ANUSPLIN") %>%
  addRasterImage(climatena, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "ClimateNA (climr)") %>%
  addRasterImage(wrfUsask, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WRF (USask)") %>%
  addCircleMarkers(lng = ~Long, lat = ~Lat, color="black", fillColor = ~ ColPal(stn.data), opacity = 1, fillOpacity = 1, popup = labels, radius=6, weight=2, group = "Stations") %>%
  addLayersControl(
    baseGroups = c("basemap", "sat photo"),
    overlayGroups = c("WRF (USask)", "ClimateNA (climr)", "ANUSPLIN", "WorldClim (1971-2000)", "CHELSA", "Daymet", "W. Can. PRISM", "BC PRISM", "Stations"),
    options = layersControlOptions(collapsed = FALSE)
  )
map
# 
# leaflet map (basics)
labels <- paste(stn.info$Name, "(El. ", stn.info$Elevation, "m)", sep="")
map <- leaflet(stn.info) %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  addRasterImage(prism.bc, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "BC PRISM") %>%
  # addRasterImage(wrfUsask, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WRF (USask)") %>%
  addCircleMarkers(lng = ~Long, lat = ~Lat, color="black", fillColor = ~ ColPal(stn.data), opacity = 1, fillOpacity = 1, popup = labels, radius=6, weight=2, group = "Stations") %>%
  addLayersControl(
    baseGroups = c("basemap", "sat photo"),
    # overlayGroups = c("WRF (USask)", "BC PRISM", "Stations"),
    overlayGroups = c("BC PRISM", "Stations"),
    options = layersControlOptions(collapsed = FALSE)
  )
map
# 
# # leaflet map (elevation)
# map <- leaflet() %>%
#   addTiles(group = "basemap") %>%
#   addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
#   addRasterImage(dem.bc, colors =terrain.colors(99), opacity = 1, maxBytes = 6 * 1024 * 1024, group = "BC PRISM") %>%
#   addRasterImage(dem.usask, colors =terrain.colors(99), opacity = 1, maxBytes = 6 * 1024 * 1024, group = "WRF (USask)") %>%
#   addLayersControl(
#     baseGroups = c("basemap", "sat photo"),
#     overlayGroups = c("WRF (USask)", "BC PRISM"),
#     options = layersControlOptions(collapsed = FALSE)
#   )
# map


