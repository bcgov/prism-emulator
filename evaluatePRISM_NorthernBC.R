
# EDA plots comparing BC PRISM Temperature to PRISM station data. 
# Colin Mahony 
# May 31, 2024

library(terra)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "Pr")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


regions <- c("NorthernBC", "SouthernBC")

region <- regions[1]
for(region in regions){
  plotarea <- if(region == "NorthernBC") ext(c(-132, -126, 56, 59)) else ext(c(-124, -118, 49, 52)) 
  
  dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
  dem.prism <- rast(paste(dir, "PRISM_dem.asc", sep=""))
  
  # function for rounding to the nearest odd integer. 
  round_to_odd <- function(x) {
    if (round(x) %% 2 == 0) {
      return(floor(x/2) * 2 + 1)
    } else {
      return(round(x))
    }
  }
  
  sources <- c("prism", "usask")
  source <- "prism"
  for(source in sources){
    
    if(source=="prism"){
      # PRISM dem
      dem <- dem.prism
      dem <- crop(dem, plotarea)
      # plot(dem)
    } else if(source=="usask"){
      # WRF dem
      dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/HGT/", sep="")
      dem <- rast(paste(dir, "HGT_regrid.nc", sep=""))
      dem <- dem[[1]]
      dem <- project(dem, dem.prism)
      dem <- crop(dem, plotarea)
      # plot(dem)
    }
    
    #---------------------------------
    # Climate data
    #---------------------------------
    
    e <- 1
    for(e in 1:3){
      
      m <- 1
      
      for(m in c(1,3,7,10)){
        
        # load the source STATION data for the BC prism
        dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
        stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
        for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
        stn.info <- stn.info[-which(El_Flag=="@"),]
        stn.data <- stn.info[,get(month.abb[m])]
        stn.info <- stn.info[is.finite(stn.data),]
        stn.info <- stn.info[-which(stn.info$Elevation == -9999),]
        stn.pts <- vect(stn.info, geom=c("Long", "Lat"), keepgeom=T)
        stn.pts <- crop(stn.pts, dem)
        # stn.dem <- extract(dem, stn.pts, method="bilinear")[,-1]
        # stn.pts$el.diff <- stn.dem-stn.pts$Elevation
        # stn.pts <- stn.pts[is.finite(stn.pts$el.diff)]
        # hist(stn.pts$el.diff)
        # cutoff <- max(abs(quantile(stn.pts$el.diff, c(0.05, 0.95)))) # tolerance (in meters) for the difference between the DEM and the station elevation
        # stn.pts <- stn.pts[0-cutoff < stn.pts$el.diff & stn.pts$el.diff < cutoff] # prune stations that have an elevation difference beyond the tolerance
        
        if(source=="prism"){
          # BC prism
          dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
          clim <- rast(paste(dir, list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],".*._", m, ".tif", sep="")), sep=""))
          clim <- crop(clim, plotarea) #crop to a smaller test area
          values(clim)[!is.finite(values(clim))] <- NA
        } else if(source=="usask"){
          # USask wrf
          dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/"
          clim <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
          clim <- project(clim, dem.prism)
          clim <- crop(clim, plotarea) #crop to a smaller test area
          if(elements[e] == "Pr") clim <- clim*monthdays[m] #convert to monthly precip
          values(clim)[!is.finite(values(clim))] <- NA
        }
        
        
        #============================================
        ## smoothed and residual elevation
        latFactor <- cos(mean(ext(dem)[3:4])*(pi/180)) # latitudinal correction for longitudinal length of cell
        w <- 77
        
        dem.d <- focal(dem, w=c(w, round_to_odd(w*latFactor)), fun=mean, na.rm=T)
        dem.o <- dem
        dem.r <- dem.o-dem.d
        
        #============================================
        ## smoothed and residual climate data
        if(elements[e] == "Pr") clim <- log2(clim)
        clim.d <- focal(clim, w=c(w, round_to_odd(w*latFactor)), fun=mean, na.rm=T)
        clim.o <- clim
        clim.r <- clim.o-clim.d
        clim.r <- if(elements[e] == "Pr") 100*(2^(clim.o-clim.d)-1) else clim.o-clim.d
        # breaks=seq(0-max(values(abs(clim.d), na.rm=T)), max(values(abs(clim.d), na.rm=T)), diff(range(as.vector(values(clim.d)), na.rm=T))/200)
        # plot(clim.d, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaks)), breaks=breaks)
        # breaks=seq(0-max(values(abs(clim.r), na.rm=T)), max(values(abs(clim.r), na.rm=T)), diff(range(as.vector(values(clim.r)), na.rm=T))/200)
        # plot(clim.r, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaks)), breaks=breaks)
        
        
        #============================================
        # scatterplots for selected areas (using residual values) 
        #============================================
        
        
        # window centroids
        point.look <- matrix(c(plotarea[1]+1.5, plotarea[4]-0.75, plotarea[2]-1.5, plotarea[4]-0.75, plotarea[1]+1.5, plotarea[3]+0.75, plotarea[2]-1.5, plotarea[3]+0.75), nrow = 4, byrow = T)
        box.size <- res(dem)*w
        
        pal <- rev(hcl.colors(99, "RdBu")) 
        if(elements[e] != "Pr") title <- paste(month.name[m], elements[e], "residuals (Kelvins above regional average temperature)") else title <- paste(month.name[m], elements[e], "residuals (% above regional average precipitation)")
        
        par(mfrow=c(1,1))
        png(filename=paste("results/", region, "_lapserates_residuals_", source, "_", elements[e], "_", monthcodes[m], ".png", sep=""), type="cairo", units="in", width=8, height=9, pointsize=14, res=600)
        
        # mat <- cbind(c(1,2,4), c(1,3,5)); layout(mat)
        mat <- rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(2,3,4,5)); layout(mat)
        par(mar=c(2,1,1,0.5), mgp=c(1.75, 0.25, 0))
        
        X <- clim.r
        lim <- max(values(clim.r), na.rm=T)
        X[X>lim] <- lim
        X[X < 0-lim] <- 0-lim
        plot(crop(X, plotarea), legend=TRUE, type="continuous", col=pal, 
             main=title)
        points(stn.pts)
        for(i in 1:dim(point.look)[1]){
          box.look <- ext(c(point.look[i,1]+box.size*c(-1, 1)/latFactor, point.look[i,2]+box.size*c(-1, 1)))
          plot(box.look, add=T)
          text(point.look[i,1]+box.size*0.8/latFactor, point.look[i,2]+box.size*0.8, i, cex=3, font=2)
        }
        
        par(mar=c(4,4,1,1), mgp=c(1.75, 0.25, 0))
        for(i in 1:dim(point.look)[1]){
          box.look <- ext(c(point.look[i,1]+box.size*c(-1, 1)/latFactor, point.look[i,2]+box.size*c(-1, 1)))
          dem.r.crop <- crop(dem.r, box.look)
          clim.r.crop <- crop(clim.r, dem.r.crop) # had to do it this way because for some reason this was cropping with one row less than the dem
          dem.look <- values(dem.r.crop)
          clim.look <- values(clim.r.crop)
          
          plot(dem.look, clim.look, yaxt="n", pch=16, cex=0.15, col="gray", xlab="Residual elevation (m)", tck = -0.01, ylab = if(elements[e] == "Pr") "Precipitation residual (%)" else "Residual temperature (K)")
          axis(2, at=pretty(clim.look), labels=pretty(clim.look), tck=-0.01, las=2)
          lm <- lm(clim.look~dem.look)
          if(is.finite(lm$coefficients[2])) abline(lm)
          mtext(i, adj=0.95, line = -2.5, side=3, cex=1.5, font=2)
          unit <- if(elements[e] != "Pr") "K" else "%"
          mtext(paste0("lapse = ", round(lm$coefficients[2]*1000, 1), " ", unit, "/km"), adj=0.05, line = -1.5, side=1, cex = 0.6)
        }
        
        dev.off()
        
        #============================================
        # scatterplots for selected areas (using raw values) 
        #============================================

        par(mfrow=c(1,1))
        png(filename=paste("results/", region, "_lapserates_raw_", source, "_", elements[e], "_", monthcodes[m], ".png", sep=""), type="cairo", units="in", width=8, height=9, pointsize=14, res=600)
        
        # mat <- cbind(c(1,2,4), c(1,3,5)); layout(mat)
        mat <- rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(2,3,4,5)); layout(mat)
        par(mar=c(2,1,1,0.5), mgp=c(1.75, 0.25, 0))
        
        if(elements[e] != "Pr") pal <- rev(hcl.colors(99, "RdBu")) else pal <- rev(hcl.colors(99, "YlGnBu")) 
        if(elements[e] != "Pr") title <- paste(month.name[m], elements[e], "(degrees Celsius)", sep=" ") else title <- paste(month.name[m], elements[e], "(log2(mm))", sep=" ")
        
        plot(crop(clim.o, plotarea), legend=TRUE, type="continuous", col=pal, 
             main=title)
        points(stn.pts)
        for(i in 1:dim(point.look)[1]){
          box.look <- ext(c(point.look[i,1]+box.size*c(-1, 1)/latFactor, point.look[i,2]+box.size*c(-1, 1)))
          plot(box.look, add=T)
          text(point.look[i,1]+box.size*0.8/latFactor, point.look[i,2]+box.size*0.8, i, cex=3, font=2)
        }
        
        par(mar=c(4,3,1,1), mgp=c(1.75, 0.25, 0))
        for(i in 1:dim(point.look)[1]){
          box.look <- ext(c(point.look[i,1]+box.size*c(-1, 1)/latFactor, point.look[i,2]+box.size*c(-1, 1)))
          stn.pts.look <- crop(stn.pts, box.look)
          el.look <- as.data.frame(stn.pts.look[,which(names(stn.pts.look)=="Elevation")])[,1] 
          stn.look <- as.data.frame(stn.pts.look[,which(names(stn.pts.look)==month.abb[m])])[,1] # vector of the climate variable of interest
          stn.look <- if(elements[e]=="Pr") log2(stn.look) else stn.look/10
          dem.o.crop <- crop(dem.o, box.look)
          clim.o.crop <- crop(clim.o, dem.o.crop) # had to do it this way because for some reason this was cropping with one row less than the dem
          dem.look <- values(dem.o.crop)
          clim.look <- values(clim.o.crop)
          
          plot(dem.look, clim.look, yaxt="n", xlim=range(c(el.look, dem.look)), ylim=range(c(stn.look, clim.look)), pch=16, cex=0.15, col="gray", xlab="Elevation (m)", tck = -0.01, ylab = if(elements[e] == "Pr") "Precipitation (mm)" else "Temperature (K)")
          if(elements[e] == "Pr"){
            axis(2, at=log2(pretty(2^clim.look)), labels=pretty(2^clim.look), tck=-0.01, las=2)
          } else axis(2, at=pretty(clim.look), labels=pretty(clim.look), tck=-0.01, las=2)
          points(el.look, stn.look, pch=16, cex=1.5, col="blue")
          lm <- lm(stn.look~el.look)
          if(is.finite(lm$coefficients[2])) abline(lm, col="blue")
          unit <- if(elements[e] != "Pr") "K" else "%"
          mtext(paste0("lapse =", round(if(elements[e] == "Pr") 100*(2^(lm$coefficients[2])-1)*1000 else lm$coefficients[2]*1000, 1), " ", unit, "/km"), col="blue", adj=0.05, line = -1.5, side=1, cex = 0.6)
          lm <- lm(clim.look~dem.look)
          if(is.finite(lm$coefficients[2])) abline(lm)
          mtext(paste0("lapse =", round(if(elements[e] == "Pr") 100*(2^(lm$coefficients[2])-1)*1000 else lm$coefficients[2]*1000, 1), " ", unit, "/km"), adj=0.05, line = -2.5, side=1, cex = 0.6)
          mtext(i, adj=0.95, line = -2.5, side=3, cex=1.5, font=2)
        }
        
        dev.off()
        
        print(m)
      }
      print(e)
    }
    print(source)
  }
  print(region)
}



