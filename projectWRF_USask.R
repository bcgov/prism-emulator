## determining correct projection of the USask WRF simulations
## Colin Mahony colin.mahony@gov.bc.ca
## May 24, 2024
## NOTE: i had to use raster because terra was acting buggy on the projection

# =======================================
# Step 1: determine the extent of the WRF projection

library(raster)
library(sp)
library(ncdf4)

filename <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/HGT.nc"

P4S.LCC <- CRS("+proj=lcc +lat_1=55 +lat_2=65 +lat_0=60.0000038146973 +lon_0=-117 +units=m +datum=WGS84 +no_defs")
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

# Grab the lat and lon from the data
lat <- raster(filename, varname='XLAT')
lon <- raster(filename, varname='XLONG')

# Convert to points and match the lat and lons
plat <- rasterToPoints(lat)
plon <- rasterToPoints(lon)
lonlat <- cbind(plon[,3], plat[,3])

# Specify the lonlat as spatial points with projection as long/lat
lonlat <- SpatialPoints(lonlat, proj4string = P4S.latlon)
plonlat <- spTransform(lonlat, CRSobj = P4S.LCC)

extent.usask <- extent(plonlat) # determine the extent

# =======================================
# step 2: Example implementation of the projection to lonlat

library(raster)

crs.usask <- "+proj=lcc +lat_1=55 +lat_2=65 +lat_0=60.0000038146973 +lon_0=-117 +units=m +datum=WGS84 +no_defs"
extent.usask <- c(-1281382.77433409, 1281382.29708051, -1397506.2674158, 1400375.6158319)

# for the WRF terrain height
dem <- raster("//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/HGT.nc")
crs(dem) <- crs.usask
extent(dem) <- extent.usask
plot(dem)
dem.latlon <- projectRaster(dem, crs="+proj=longlat +datum=WGS84", res = 0.02) # set a higher resolution to avoid losing resolution in the south. 
plot(dem.latlon)
writeRaster(dem.latlon, "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/HGT.latlon.tif", overwrite=TRUE)

# for a WRF 2D field (arbitrarily using EPSG instead)
test <- raster("//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/wrf_tmin_sample_dt.nc")
crs(test) <- crs.usask
extent(test) <- extent.usask
test.latlon <- projectRaster(test, crs="EPSG:4326", res = 0.02)
plot(test.latlon)


# =======================================
# step 3: test that the projection is correct (it's not)

dem.prism <- rast("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/PRISM_dem.asc")
dem.wrf <- rast("//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/HGT.latlon.tif")

ext <- ext(-124, -120, 50, 52)
dem.prism <- crop(dem.prism, ext)
dem.wrf <- crop(dem.wrf, ext)
dem.wrf.proj <- project(dem.wrf, dem.prism)

range <- range(c(values(dem.prism), values(dem.wrf)), na.rm=T)
breaks <- seq(range[1], range[2], 100)
plot(dem.prism, ext=extent(-124, -120, 50, 52), breaks=breaks)
plot(dem.wrf, add=T, breaks=breaks)

     