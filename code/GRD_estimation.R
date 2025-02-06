
######################## Model Estimation #######################

# 0 Load packages and directories ----
library(tidyverse)
library(terra)
library(rdrobust)

# Parameters
b.width <- 10000
timeperiod <- c(2000:2022)
targ.res <- 500

# Load data
setwd("/mnt/SCRATCH/maltela/PAsDRC")
s1aTMFDef<-rast("s1aTMFdef.tif")
s1bTMFDegr<-rast("s1bTMFdegr.tif")
s2PAs <- vect("s2PAs_joined.shp")

s2PAs <- project(s2PAs,crs(s1aTMFDef))
s2PAs$PA <- c(1:nrow(s2PAs))
bpts <- vect("bpts15000.shp")

#Load functions
  DefPanelPA <- function(defRast,PAs,bandwidth = 10000,years = c(2000:2022),aggregate = 1) {
  
  #Define raster CRS as projection standard
  projCRS <- crs(defRast)
  
  #Generate buffer around PA and crop raster accordingly 
  bufferOut <- buffer(PAs, width = bandwidth)
  bufferIn <- buffer(PAs, width = -bandwidth)
  buffer <- bufferOut - bufferIn
  
  defRast <- crop(defRast,buffer)
  
  if (aggregate !=1) {
    treatRast <- aggregate(defRast,fact=aggregate,fun="mean")
    #treatRast <- rast(ext(treatRast),vals=0,crs=treatRast,nrows=nrow(treatRast),ncols = ncol(treatRast))
    #  treatRast <- rast(ext(origRaster),vals=0,crs=origRaster,nrows=nrow(origRaster),ncols = ncol(origRaster))
  } else {
    treatRast <- defRast
  }
  
  treatRast <-rasterize(PAs,treatRast,field = "PA",background = 0,touches = F, input = T,cover = F)
  treatRast <- mask(treatRast,buffer)
  treatDF <- as.data.frame(treatRast,xy=T,cells=T)
  colnames(treatDF) <- c("cell","x","y","PA")
  
  # Calculate annual deforestation
  for (y in years) {
    print(y)
    defRast.y <- terra::ifel(defRast==y,1,0)
    if (aggregate!=1) {
      defRast.y <- aggregate(defRast.y,fact=aggregate,fun="mean")
    }
    defRast.y <- mask(defRast.y,buffer)
    tempDefDF <- as.data.frame(defRast.y,cells=T)
    colnames(tempDefDF) <- c("cell","def")
    tempDefDF$year <- y
    if (y > min(years)) {
      defDF <- merge(defDF,tempDefDF,by = c("cell","year","def"),all.x=T,all.y=T)
    } else {
      defDF <- tempDefDF
      #panelDF <- merge(panelDF,tempDefDF, by= c("cell"))
    }
    if (exists("distanceDF")) {
      defDF <- merge(defDF,distanceDF,by=c("cell"))
    } 
  }
  
  panelDF <- merge(treatDF,defDF,by=c("cell"))
  return(panelDF)
}

  estGRD <- function (depvar="def",data, cutoff.points.v, maxbandwidth = 10000, years = 2000:2020, minobs = 10, bwfix_m = NA, 
                      sparse.exclusion = FALSE, ...) 
  {
    colnames(data)[which(colnames(data)==depvar)] <- "outcome"
    n <- nrow(data)
    np <- nrow(cutoff.points.v)
    bwfix <- F
    if (!is.na(bwfix_m)) {
      bwfix <- TRUE
    }
    
    if (sparse.exclusion == TRUE) {
      Kt = function(u) {
        (1 - abs(u)) * (abs(u) <= 1)
      }
    }
    columnames <- c("bpoint","year","Estimate", "SE_Conv", "SE_Rob", 
                    "p_Conv", "p_Rob", "Ntr", "Nco", "bw_l", "bw_r", "CI_Conv_l", 
                    "CI_Conv_u", "CI_Rob_l", "CI_Rob_u","Est_l","Est_r")
    results <- data.frame(matrix(ncol = length(columnames), nrow = np*length(years)))
    colnames(results) <- columnames
    #data[["dist2cutoff"]] <- 0
    row.ind <- 1
    #cutoff.points.v <- terra::vect(cutoff.points)
    cells <- data %>%
      group_by(cell) %>%
      summarise(PA = mean(PA),x = mean(x),y = mean(y))
    
    cells <- terra::vect(cells, geom=c("x", "y"),crs=cutoff.points.v)
    for (i in 1:np) { # Loop over border points
      if (bwfix == TRUE) {
        bw <- bwfix_m
      }
      PAnr <- cutoff.points.v$PA[i]
      
      buff.cutoff<-terra::buffer(cutoff.points.v[i],maxbandwidth)
      cells.in.buff <- cells %>% 
        terra::crop(buff.cutoff) %>%
        tidyterra::filter(PA %in% c(0,PAnr))
      
      cells.in.buff$dist2cutoff <- distance(cells.in.buff,cutoff.points.v[i,])
      cells.in.buff <- as.data.frame(cells.in.buff)
      dataCrop <- data %>%
        filter(cell %in% cells.in.buff$cell,year %in% years) %>%
        left_join(cells.in.buff,by=c("cell","PA")) %>%
        tidyterra:: mutate(treated = if_else(PA == PAnr,1,0),dist2cutoff = if_else(treated==0,dist2cutoff*-1,dist2cutoff))
      
      #dataCrop.v <- st_as_sf(dataCrop.v)
      #dataCrop.v[["dist2cutoff"]] <- distance(t1[1],t1) # SWITH DISTANCE CALCULATION TO TERRA
      #dataCrop.v[["dist2cutoff"]][dataCrop.v[[treated]] == 0] <- dataCrop.v[["dist2cutoff"]][dataCrop.v[[treated]] == 
      #                                                                                        0] * (-1)
      
      for (yr in years) {
        data.yr <- dataCrop %>% filter(year == yr)
        err <- FALSE
        if (bwfix == TRUE) {
          rdrob_bwflex <- tryCatch(rdrobust::rdrobust(y=data.yr$outcome, 
                                                      x = data.yr$dist2cutoff, c = 0, h = bw, ...), 
                                   error = function(e) {
                                     err <<- TRUE
                                   })
        }
        else {
          rdrob_bwflex <- tryCatch(rdrobust::rdrobust(y=data.yr$outcome, 
                                                      x = data.yr$dist2cutoff, c = 0, ...), error = function(e) {
                                                        err <<- TRUE
                                                      })
        }
        if (err) {
          message("Skipped one boundary point due to error, possibly not enough observations in local neighbourhood. Check vignette for FAQ!")
          next
        }
        results[row.ind, "bpoint"] <- i
        results[row.ind, "year"] <- yr
        results[row.ind, "Estimate"] <- round(rdrob_bwflex[["coef"]][["Conventional", 
                                                                      1]], 2)
        results[row.ind, "SE_Conv"] <- round(rdrob_bwflex[["se"]][["Conventional", 
                                                                   1]], 2)
        results[row.ind, "SE_Rob"] <- round(rdrob_bwflex[["se"]][["Robust", 
                                                                  1]], 2)
        results[row.ind, "p_Conv"] <- round(rdrob_bwflex[["pv"]][["Conventional", 
                                                                  1]], 2)
        results[row.ind, "p_Rob"] <- round(rdrob_bwflex[["pv"]][["Robust", 
                                                                 1]], 2)
        results[row.ind, "Ntr"] <- rdrob_bwflex[["N_h"]][[1]]
        results[row.ind, "Nco"] <- rdrob_bwflex[["N_h"]][[2]]
        results[row.ind, "bw_l"] <- round(rdrob_bwflex[["bws"]][["h", 
                                                                 1]]/1000, 1)
        results[row.ind, "bw_r"] <- round(rdrob_bwflex[["bws"]][["h", 
                                                                 2]]/1000, 1)
        results[row.ind, "CI_Conv_l"] <- round(rdrob_bwflex$ci["Conventional", 
                                                               1], 2)
        results[row.ind, "CI_Conv_u"] <- round(rdrob_bwflex$ci["Conventional", 
                                                               2], 2)
        results[row.ind, "CI_Rob_l"] <- round(rdrob_bwflex$ci["Robust", 
                                                              1], 2)
        results[row.ind, "CI_Rob_u"] <- round(rdrob_bwflex$ci["Robust", 
                                                              2], 2)
        results[row.ind, "Est_l"] <- round(rdrob_bwflex[["tau_bc"]][[1]],2)
        results[row.ind, "Est_r"] <- round(rdrob_bwflex[["tau_bc"]][[2]],2)
        row.ind <- row.ind + 1
      }
    }
    results[["Ntr"]] <- as.integer(results[["Ntr"]])
    results[["Nco"]] <- as.integer(results[["Nco"]])

    results %>% dplyr::filter(.data$Nco > minobs & .data$Ntr > 
                                minobs)
  }

# Calculate deforestation panel from raster
aggr.fact <- targ.res/res(s1aTMFDef)[[1]]
i1aDefPanel <- DefPanelPA(defRast = s1aTMFDef,
                          PAs = s2PAs,
                          bandwidth = b.width,
                          years = timeperiod,
                          aggregate = aggr.fact)

save(i1aDefPanel,file = "i1aDefPanel500.Rdata")

i1bDegrPanel <- DefPanelPA(defRast = s1bTMFDegr,
                          PAs = s2PAs,
                          bandwidth = b.width,
                          years = timeperiod,
                          aggregate = aggr.fact)

save(i1bDegrPanel,file = "i1bDegrPanel500.Rdata")

r1adef <- estGRD(depvar = "def",
                data = i1aDefPanel,
                cutoff.points = bpts,
                minobs = 0,
                years = 2000:2022)

save(r1adef,file="r1adef500_15000.Rdata")

r1bdegr <- estGRD(depvar = "def",
                data = i1bDegrPanel,
                cutoff.points = bpts,
                minobs = 0,
                years = 2000:2022)

save(r1def,file="r1bdegr500_15000.Rdata")
