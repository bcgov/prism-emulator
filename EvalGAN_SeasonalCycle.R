# Evaluation of GAN climatological model in generating a credible annual cycle

library(terra)
library(data.table)


monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("tmin", "tmax", "prec")
element.names <- c("mean\ndaily minimum temperature (\u00b0C)", "mean\ndaily maximum temperature (\u00b0C)", "precipitation (mm)")

studyarea <- ext(c(-137, -120, 59, 64.5))

dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/climr_mosaic/", sep="")
dem <- rast(paste(dir, "climr_mosaic_dem_800m.tif", sep=""))
dem <- crop(dem, studyarea)

# -----------------------------------------
# Assemble the monthly files into a 12-month raster stack
for(e in 2:length(elements)){
  
  for(m in 1:12){
    
    dir <- paste0("//objectstore2.nrs.bcgov/ffec/Mosaic_Yukon/operational/WorldClim/", elements[e], "/", tolower(month.abb[m]), "/Predictions/")
    temp <- rast(paste(dir, list.files(dir, pattern=".*.nc"), sep=""))
    assign(paste(elements[e], "v2024", sep="."), if(m==1) temp else c(get(paste(elements[e], "v2024", sep=".")), temp))
    
    
    temp <- rast(paste("//objectstore2.nrs.bcgov/ffec/Mosaic_Yukon/Tirion/Results/foundational_model/", elements[e], "/Model1/", tolower(month.abb[m]), "/", tolower(month.abb[m]), "_fullregion_masked.nc", sep=""))
    assign(paste(elements[e], "v2025", sep="."), if(m==1) temp else c(get(paste(elements[e], "v2025", sep=".")), temp))
    
    print(monthcodes[m])  
  }
  print(elements[e])
}

# -----------------------------------------
# load the source STATION data for the BC prism
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
stn.info <- stn.info[-which(El_Flag=="@"),]
stn.info <- stn.info[complete.cases(stn.info[, ..month.abb])]
stn.info <- stn.info[Lat>60 & Lat<60 & Long > -114]


# -----------------------------------------
# plot comparison of station data
stns <- stn.info$St_ID
for(stn in stns){
  par(mfrow=c(2,1), mar=c(3,4,2,1), mgp=c(2, 0.25, 0))
  
  plot(dem)
  pt <- stn.info[St_ID == stn, c("Long", "Lat")]
  points(pt, cex=1, pch=16)
  text(pt, stn, cex=0.95, pos=4)
  
  stn.data <- as.vector(unlist(stn.info[St_ID == stn, ..month.abb]))
  v2024 <- as.vector(unlist(extract(get(paste(elements[e], "v2024", sep=".")), stn.info[St_ID == stn, c("Long", "Lat"), with = FALSE])[-1]))
  v2025 <- as.vector(unlist(extract(get(paste(elements[e], "v2025", sep=".")), stn.info[St_ID == stn, c("Long", "Lat"), with = FALSE])[-1]))
  
  plot(stn.data, type="l", ylim=range(c(stn.data, v2024, v2025)), lwd=2, xlab="", xaxt="n", ylab=element.names[e], main=paste("Station:", stn))
  axis(1, at=1:12, labels = month.abb, las=2, tck=0)
  lines(v2025, col="dodgerblue")
  lines(v2024, col="gray", lty=2)
  legend("topleft", legend=c("Station data", "GAN2025", "GAN2024"), lwd=c(2,1,1), lty=c(1,1,2), col=c("black", "dodgerblue", "gray"), bty="n")

  print(stn)
  }

par(mfrow=c(1,1))
r <- get(paste(elements[e], "v2025", sep="."))
plot(r[[1]])
points(stn.info$Long, stn.info$Lat, cex=0.9)


dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/climr_mosaic/", sep="")
dem <- rast(paste(dir, "climr_mosaic_dem_800m.tif", sep=""))
dem <- crop(dem, r)
plot(dem)
points(stn.info$Long, stn.info$Lat, cex=0.9)

evalarea <- ext(c(-141, -120, 60, 65))
plot(evalarea, add=T)
