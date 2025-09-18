library(terra)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

studyarea <- ext(c(-137, -119, 56, 64))
trainarea <- ext(c(-134, -128, 57, 60))

# load the AK PRISM  data for the variable
# prism.ak <- rast('//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_AK/ak_tmax_1981_2010.03.asc')
# plot(prism.ak)

# load the GAN Prediction
yukon_pred <- rast('//objectstore2.nrs.bcgov/ffec/Mosaic_Yukon/GAN_Preds_June18_stand.tif')
# yukon_pred <- project(yukon_pred, prism.bc)
# plot(yukon_pred)

# load the BC PRISM  data for the variable
prism.bc <- rast('//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/tmax_monClim_PRISM_historical_198101-201012_3.tif')
prism.bc <- crop(prism.bc, yukon_pred)
plot(prism.bc)
plot(testarea, add=T)

# load the Random Forest Prediction
yukon_RF <- rast('C:/Users/CMAHONY/OneDrive - Government of BC/Projects/2024_intercomparison/outputs/Composite_worldclim/composite_WNA_1981_2010_Tmax03.tif')
yukon_RF <- crop(yukon_RF, yukon_pred)

# color scheme
combined <- c(values(prism.bc), values(yukon_pred), values(yukon_RF))
combined <- combined[is.finite(combined)]
inc=diff(range(combined))/500
breaks=seq(quantile(combined, 0)-inc, quantile(combined, 1)+inc, inc)
ColScheme <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColPal <- colorBin(ColScheme, bins=breaks, na.color = "white")
ColPal.raster <- colorBin(ColScheme, bins=breaks, na.color = "transparent")


# leaflet map
map <- leaflet() %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  # addRasterImage(prism.ak, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "AK PRISM") %>%
  addRasterImage(yukon_pred, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "GAN") %>%
  addRasterImage(yukon_RF, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "RF") %>%
  addRasterImage(prism.bc, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "BC PRISM") %>%
  addLayersControl(
    overlayGroups = c("GAN", "RF", "BC PRISM"),
    options = layersControlOptions(collapsed = FALSE)
  )
map
