#################GRDD Functions and plotting of results#########################

# Plot 

# - translate distance() function from sf to raster

# Calculate Distance from and to PA-#######################

DefPanelPA <- function(defRast,PAs,bandwidth = 10000,years = c(2000:2022),aggregate = 1,add_distances = T) {
  
  # Test PA:
  #selectedPA <- 64
  #origRaster <- s1aTMFDef
  
  #Function to calculate the distance of each raster cell to the closest PA
  #Cells inside are coded with negative distance, outside with positive distance
  
  #rasgrid:   raster layer with forest cells
  #PAs:       Shapefile with PA polygons
  #PADF:      Data frame with PA characteristics
  #bandwidth: Bandwith around polygon in within which to calculate distance
  #selectedPA: Number indicator of PA to be calculated
  
  #Returns raster stack with the following layers:
  # 1) Distance to closest PA boundary
  # 2) Binary deforestation indicator after PA establishment
  # 3) Binary deforestation indicator before PA establishment (if available)
  # 4) PA number
  
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
  #treatRast <- mask(treatRast,PAs,updatevalue=1,inverse=T,touches=F)
  treatRast <- mask(treatRast,buffer)
  treatDF <- as.data.frame(treatRast,xy=T,cells=T)
  colnames(treatDF) <- c("cell","x","y","PA")
  
if (add_distances == T) {
  distanceOut<-distance(treatRast,PAs,rasterize=T) #->Distance of NA cells to closest non-NA cell
  
  #Now calculate the same for cells inside
  #Create a polygon of cropped raster
  PAinv<-vect(ext(treatRast),crs=projection)
  
  #Generate polygon of inverse PA shape
  PAinv<-symdif(PAinv,PAs)
  #Calculate distance
  distanceIn<-distance(treatRast,PAinv,rasterize=T)
  distanceIn<-distanceIn*(-1)
  
  #Write distances into one 'distances' grid
  distanceIn<-ifel(distanceIn==0,NA,distanceIn)
  distances<-cover(distanceIn,distanceOut)
  distances <- mask(distances,buffer)
  distanceDF <- as.data.frame(distances,cells=T)
  colnames(distanceDF) <- c("cell","distance")
}
  #distanceDF <- as.data.frame(distances,xy=T,cells=T)
  #colnames(distanceDF) <- c("cell","x","y","distance")


  #panelDF <- data.frame(cell = rep(distanceDF$cell,each = length(years)),year = rep(years,times=nrow(distanceDF)))

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

AnnualRastersToPanel <- function(years,file.pattern,folder.directory,PAs,bandwidth = 15000) {

  for (y in years) {
    print(y)
    source.rast <- rast(here(folder.directory,glue::glue("{file.pattern}{y}.tif")))
    if (y == years[[1]]) {
      target.crs <- crs(source.rast)
      PAs <- project(PAs,target.crs)
      #Generate buffer around PA and crop raster accordingly 
      buffer.out <- buffer(PAs, width = bandwidth)
      buffer.in <- buffer(PAs, width = -bandwidth)
      buffered.PAs <- buffer.out - buffer.in
      
      #Create raster to indicate the cells that lie within PAs
      treat.rast <-rasterize(PAs,source.rast,field = "PA",background = 0,touches = F, input = T,cover = F)
      treat.rast <- mask(treat.rast,buffered.PAs)
      panel.df <- as.data.frame(treat.rast,xy=T,cells=T)
      colnames(panel.df) <- c("cell","x","y","PA")
    }
    
    #
    target.rast <- crop(source.rast,buffered.PAs)
    target.rast <- mask(target.rast,buffered.PAs)
    target.df <- as.data.frame(target.rast,cells=T)
    colnames(target.df) <- c("cell",paste("fcover",y,sep=""))
    panel.df <- rbind(panel.df,target.df)
  }
  return(panelDF)
  }

# Covariate balance---------------------------------------------------------------------

#Function to test balance of covariates

balanceTest <- function (distgrid = cells, covariates, polynomial = 0, wmin=100, level = 0.15) {
  
  if (!("dist2cutoff" %in% names(distgrid))) {
    grid.df <- grid.df %>% group_by(cell) %>% summarise(x = mean(x), y = mean(y), PA = mean(PA)) %>% ungroup()
    cells <- terra::vect(grid.df, geom=c("x", "y"),crs=s1cTMFCov)
    boundaries <- s2PAs %>% as.lines() %>% aggregate()
    
    cells$dist2cutoff <- distance(cells,boundaries)
    cells %>% tidyterra:: mutate(dist2cutoff = if_else(PA>0,dist2cutoff*-1,dist2cutoff))
    
    for (rast in covariates) {
      covrast <- rast(here(dataInt,glue::glue("{rast}_500.tif")))
      covrast <- project(covrast,crs(cells))
      print(glue::glue("Extract values for {rast}"))
      cells[[rast]] <- extract(covrast,cells)[,2]
    }
  }
  
  distances <- cells$dist2cutoff
  covMX <- as.matrix(cells[,covariates])
  balanceTable <- rdwinselect(R = distances,X = covMx,p = polynomial, wmin = wmin, level = level) 
}
  
# Run GRD---------------------------------------------------------------------

#Function to run GRD

estGRD <- function (depvar="def",data, cutoff.points, maxbandwidth = 10000, years = 2000:2020, minobs = 50, bwfix_m = NA, 
                    sparse.exclusion = FALSE, ...) 
{
  colnames(data)[which(colnames(data)==depvar)] <- "outcome"
  n <- nrow(data)
  np <- nrow(cutoff.points)
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
  row.ind <- 1
  cutoff.points.v <- terra::vect(cutoff.points)
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
  message("Skipped")
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

# Calcualte overall ATT by PA and bootstrap CIs ####################

AvgAnnualATT <- function(res.df,years = 2000:2022,blockbstrap=T) {
  
  colnames <- c("year","ATT","c.low","c.high","est_l","est_r","l.c.low","l.c.high","r.c.low","r.c.high")
  ATTs.annual <- data.frame(matrix(ncol = length(colnames), nrow = length(years)))
  colnames(ATTs.annual) <- colnames
  rowind <- 1
  for (y in years) {
    t1 <- res.df %>% filter(year==y)   
    ATT <- mean(t1$est_rob,na.rm=T)
    #ATT.ci.boot <- ATE.boot(data=t1,y="Estimate")
    #boot.ci <- ATEboot.cluster(data=t1,y="Estimate")
    #ATT_se <- sum((t1$SE_Rob^2))/length(t1$SE_Rob)
    est_l <- mean(t1$est_l,na.rm=T)
    est_r <- mean(t1$est_r,na.rm=T)
    if (blockbstrap==T) {
      boot.ci <- block.bstrap(data=t1,estimate="est_rob")
      boot.ci.l <- block.bstrap(data=t1,estimate="est_l")
      boot.ci.r <- block.bstrap(data=t1,estimate="est_r")
    } else {
      boot.ci <- bstrap(data=t1,estimate="est_rob")
      boot.ci.l <- bstrap(data=t1,estimate="est_l")
      boot.ci.r <- bstrap(data=t1,estimate="est_r")
    }
    ATTs.annual[rowind,"year"] <- y
    ATTs.annual[rowind,"ATT"] <- ATT
    ATTs.annual[rowind,"c.low"] <- boot.ci$bt.cl
    ATTs.annual[rowind,"c.high"] <-boot.ci$bt.ch
    ATTs.annual[rowind,"est_l"] <- est_l
    ATTs.annual[rowind,"est_r"] <- est_r
    ATTs.annual[rowind,"l.c.low"] <- boot.ci.l$bt.cl
    ATTs.annual[rowind,"l.c.high"] <-boot.ci.l$bt.ch
    ATTs.annual[rowind,"r.c.low"] <- boot.ci.r$bt.cl
    ATTs.annual[rowind,"r.c.high"] <-boot.ci.r$bt.ch
    rowind <- rowind+1
  }
return(ATTs.annual)
}

# Calcualte ATT by PA and bootstrap CIs ####################

ATTbyPA <- function(res.df,years = 2000:2022, PAs=NULL) {
  
  if(is.null(PAs)) PAs <- unique(res.df$PA)
  
#Discontinuities by year and PA
ATT.PAs.df <- data.frame()
for (PA.ind in PAs) {
rowind <- 1
print(glue::glue("PA {PA.ind}"))
colnames <- c("year","ATT","c.low","c.high","est_l","est_r","l.c.low","l.c.high","r.c.low","r.c.high","PA.ind")
temp.df <- data.frame(matrix(ncol = length(colnames), nrow = length(years)))
colnames(temp.df) <- colnames
 for (y in years) {
#t1 <- r1.df[[3]] %>% filter(PA_year >2004,year==y,PA==PA.ind)
tempres.df <- res.df %>% filter(year==y,PA==PA.ind)   
ATT <- mean(tempres.df$est_rob,na.rm=T)
#ATT.ci.boot <- ATE.boot(data=res.df,y="Estimate")
#boot.ci <- ATEboot.cluster(data=res.df,y="Estimate")
boot.ci <- bstrap(data=tempres.df,estimate="est_rob")
#ATT_se <- sum((res.df$SE_Rob^2))/length(res.df$SE_Rob)
est_l <- mean(tempres.df$est_l,na.rm=T)
est_r <- mean(tempres.df$est_r,na.rm=T)
boot.ci.l <- bstrap(data=tempres.df,estimate="est_l")
boot.ci.r <- bstrap(data=tempres.df,estimate="est_r")
temp.df[rowind,"year"] <- y
temp.df[rowind,"ATT"] <- ATT
temp.df[rowind,"c.low"] <- boot.ci$bt.cl
temp.df[rowind,"c.high"] <-boot.ci$bt.ch
temp.df[rowind,"est_l"] <- est_l
temp.df[rowind,"est_r"] <- est_r
temp.df[rowind,"l.c.low"] <- boot.ci.l$bt.cl
temp.df[rowind,"l.c.high"] <-boot.ci.l$bt.ch
temp.df[rowind,"r.c.low"] <- boot.ci.r$bt.cl
temp.df[rowind,"r.c.high"] <-boot.ci.r$bt.ch
temp.df[rowind,"PA.ind"] <- PA.ind
rowind <- rowind+1
 }
ATT.PAs.df <- rbind(ATT.PAs.df,temp.df)
}
return(ATT.PAs.df)
}

# Calcualte ATT before and after DBF and bootstrap CIs ####################

AvgATTBefAft <- function(res.df,blockbstrap=T) {
  
  colnames <- c("time","ATT","c.low","c.high","est_l","est_r","l.c.low","l.c.high","r.c.low","r.c.high")
  ATTs.annual <- data.frame(matrix(ncol = length(colnames), nrow = length(years)))
  colnames(ATTs.annual) <- colnames
  rowind <- 1
  res.df <- res.df %>% group_by(bpoint) %>% filter(rel_time %in% c(-1,max(rel_time)),PA_year<year) %>% ungroup()
  for (y in c("pre","post")) {
    if (y =="pre") {
      t1 <- res.df %>% group_by(bpoint) %>% filter(rel_time == -1,PA_year<year) %>% ungroup()
    } else {
      t1 <- res.df %>% group_by(bpoint) %>% filter(rel_time == max(rel_time),PA_year<year) %>% ungroup()
    }
    ATT <- mean(t1$est_rob)
    #ATT.ci.boot <- ATE.boot(data=t1,y="Estimate")
    #boot.ci <- ATEboot.cluster(data=t1,y="Estimate")
    #ATT_se <- sum((t1$SE_Rob^2))/length(t1$SE_Rob)
    est_l <- mean(t1$est_l)
    est_r <- mean(t1$est_r)
    if (blockbstrap==T) {
      boot.ci <- block.bstrap(data=t1,estimate="est_rob")
      boot.ci.l <- block.bstrap(data=t1,estimate="est_l")
      boot.ci.r <- block.bstrap(data=t1,estimate="est_r")
    } else {
      boot.ci <- bstrap(data=t1,estimate="est_rob")
      boot.ci.l <- bstrap(data=t1,estimate="est_l")
      boot.ci.r <- bstrap(data=t1,estimate="est_r")
    }
    ATTs.annual[rowind,"time"] <- y
    ATTs.annual[rowind,"ATT"] <- ATT
    ATTs.annual[rowind,"c.low"] <- boot.ci$bt.cl
    ATTs.annual[rowind,"c.high"] <-boot.ci$bt.ch
    ATTs.annual[rowind,"est_l"] <- est_l
    ATTs.annual[rowind,"est_r"] <- est_r
    ATTs.annual[rowind,"l.c.low"] <- boot.ci.l$bt.cl
    ATTs.annual[rowind,"l.c.high"] <-boot.ci.l$bt.ch
    ATTs.annual[rowind,"r.c.low"] <- boot.ci.r$bt.cl
    ATTs.annual[rowind,"r.c.high"] <-boot.ci.r$bt.ch
    rowind <- rowind+1
  }
  return(ATTs.annual)
}

# Bootstrap ATE standard errors ####################

bstrap <- function(estimate="Estimate", data, N=500, seq=NULL, cil=0.025, cih=0.975, ...){ 
  
  #store Y and X data
  data <- data[estimate]
  
  #number of rows in data
  datasize <- dim(data)[1]
  
  #create matrix to hold loess predictions
  predmat <- matrix(0, N, 1)
  
  #loop for bootstraps
  for(i in 1:N){
    
    #take sample of rows of length(X)
    xind <- sample(1:datasize, datasize, replace=T)
    #create X,Y data based on random sample
    x=data[xind,]
    #loess on random sample
    predmat[i,] <- mean(x[[1]],na.rm=T)
    
    #dim confidence interval vectors
    cih.lo <- vector()
    cil.lo <- vector()
    warn <-vector()
    j <- 1
    
    #store CI data foreach prediction of seq
    cih.lo <- quantile(predmat, cih, na.rm=T)
    cil.lo <- quantile(predmat, cil, na.rm=T)
  }#end for
  
  #store output
  b.res <- list(bt.estimates=predmat, bt.cl=cil.lo, bt.ch=cih.lo)
  return(b.res)
}

# Block bootstrap ATE standard errors ####################

block.bstrap <- function(data,estimate="Estimate",cluster = "PA",N=500) {
  cl <- unique(data[[cluster]])
  print(glue::glue("bootstrap from {length(cl)} clusters"))
  
  stat <- function(x, i) {
    # select the observations to subset based on the cluster var
    block_obs <- unlist(lapply(i, function(n) which(x[n] == data[[cluster]])))
    mean(data[block_obs,][[estimate]])
  }
  bt <- boot::boot(cl, stat, N)
  b.res <- list()
  b.res[["bt.estimates"]] <- bt$t
  b.res[["bt.cl"]] <- quantile(bt$t,0.025,na.rm=T)
  b.res[["bt.ch"]] <- quantile(bt$t,0.975,na.rm=T)
  return(b.res)
}


# Plot results and data#######################

PlotAnnualATTs <- function(ATT.df, left.right = T,years = c(2000:2022),PAs = NULL) {
  if (is.null(PAs)) {
    plots <- ggplot(data = ATTs.annual,aes(x=year, y=ATT)) + 
      geom_ribbon(aes(ymin=c.low, ymax=c.high),alpha=0.18) +#,position=position_dodge(width=0.5),color=NA) + 
      geom_ribbon(aes(ymin=l.c.low, ymax=l.c.high),alpha=0.18) +#,position=position_dodge(width=0.5),color=NA) + 
      geom_ribbon(aes(ymin=r.c.low, ymax=r.c.high),alpha=0.18) +#,position=position_dodge(width=0.5),color=NA) + 
      geom_line(position=position_dodge(width=0.5),linewidth=.8) +
      geom_line(aes(x=year,y=Est_l,color = "Est_l")) +
      geom_line(aes(x=year,y=Est_r,color = "Est_r")) +
      geom_point(size=2.3,position=position_dodge(width=0.5)) +
      geom_hline(yintercept=0, linetype="solid", color = "black") +
      theme_bw() +
      #theme(legend.position="bottom") +
      scale_color_manual(values = c("Est_l" = "#766297","Est_r" = "#64ad6b")) +
      #scale_color_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      #scale_fill_viridis(option = "turbo",discrete = T,limits=c("All land use","Farming","Settlements","Mining","Other")) +
      scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      scale_shape_manual(values = c(15,16,17,18,8)) +
      #scale_fill_continuous(guide = "colourbar") +
      xlim(c(2000,2022)) +
      xlab("Year") +
      ylab("Outcome") +
      labs(x = "Year",y="Outcome",title="All PAs <2000")
  } else {
    
  }
  for (PA in PAs) {
    plots <- list()
    plots[[p.ind]] <- ggplot(data = ATTs.annual[ATTs.annual$PA.ind==PA,],aes(x=year, y=ATT)) + 
      geom_ribbon(aes(ymin=c.low, ymax=c.high),alpha=0.18) +#,position=position_dodge(width=0.5),color=NA) + 
      geom_line(position=position_dodge(width=0.5),linewidth=.8) +
      geom_line(aes(x=year,y=Est_l,color = "Est_l")) +
      geom_line(aes(x=year,y=Est_r,color = "Est_r")) +
      geom_point(size=2.3,position=position_dodge(width=0.5)) +
      geom_hline(yintercept=0, linetype="solid", color = "black") +
      geom_vline(xintercept = year.of.est, linetype = "dotdash") +
      theme_bw() +
      #theme(legend.position="bottom") +
      scale_color_manual(values = c("Est_l" = "#766297","Est_r" = "#64ad6b")) +
      #scale_color_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      #scale_fill_viridis(option = "turbo",discrete = T,limits=c("All land use","Farming","Settlements","Mining","Other")) +
      scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      scale_shape_manual(values = c(15,16,17,18,8)) +
      #scale_fill_continuous(guide = "colourbar") +
      xlim(c(2000,2022)) +
      xlab("Year") +
      ylab("Outcome") +
      labs(x = "Year",y="Outcome",title=PA.ind,)
    
    p.ind <- p.ind + 1
  }
  return(plots)
}

####################################################

# Local randomization aproach

estLocRand <- function (depvar,data, border.segments, calculate.windows = T, bandwidth = 10000, level = 0.15, years = 2000:2020, df.crs, covariates = c("s3aAccess","s3bAgrSuit","s3cSlope","s3dAltitude","s3ePop2000"), locrand = F, ...) 
{
  colnames(data)[which(colnames(data)==depvar)] <- "outcome"
  np <- nrow(border.segments)
  
  columnames <- c("border_seg","PA","cov_ind","year","est","p_value","se_bc","se_rob","est_l","est_r","window_l","window_r","segment_length","rel_time")
  results <- data.frame(matrix(ncol = length(columnames), nrow = np*length(years)))
  colnames(results) <- columnames
  row.ind <- 0
  cells <- data %>%
    group_by(cell) %>%
    summarise(PA = mean(PA),x = mean(x),y = mean(y))
  
  cells <- terra::vect(cells, geom=c("x", "y"),crs=df.crs)
  #cells <- project(cells,border.segments)
  
  for (rast in covariates) {
    covrast <- rast(here(dataInt,glue::glue("{rast}_500.tif")))
    covrast <- project(covrast,crs(cells))
    print(glue::glue("Extract values for {rast}"))
    cells[[rast]] <- extract(covrast,cells)[,2]
  }
  
  
  for (i in 1:np) { # Loop over border segments with double frontier
    print(glue::glue("segment {i}/{np}"))
    PAnr <- border.segments$PA[i]
    yr.mine <- border.segments$year.mine0[i]
    buffed.segments <- buffer(border.segments[i],bandwidth)
    cells.in.buff <- cells %>%
        project(border.segments) %>%
        terra::crop(buffed.segments)
      
    cells.in.buff$dist2cutoff <- distance(cells.in.buff,border.segments[i,])

    # for (rast in covariates) {
    #   covrast <- rast(here(dataInt,glue::glue("{rast}_500.tif")))
    #   covrast <- project(covrast,crs(cells.in.buff))
    #   cells.in.buff[[rast]] <- extract(covrast,cells.in.buff)[,2]
    # }
    #   
    cells.in.buff <- as.data.frame(cells.in.buff)
    
    dataCrop <- data %>%
      filter(cell %in% cells.in.buff$cell,year %in% years) %>%
      left_join(cells.in.buff,by=c("cell","PA")) %>%
      tidyterra:: mutate(treated = if_else(PA == PAnr,1,0),
                         dist2cutoff = if_else(treated==0,dist2cutoff*-1,dist2cutoff))
      
    base_out <- dataCrop$outcome[dataCrop$year==(yr.mine-1)]
    #dataCrop <- dataCrop %>% mutate(diff_out = outcome - base_out)
      
    for (yr in years) {
      row.ind <- row.ind + 1
      data.yr <- dataCrop %>% filter(year == yr) %>% mutate(diff_out = outcome - base_out)
      err <- FALSE
      if(locrand == T) {
      rd.est <- tryCatch(rdlocrand::rdrandinf(Y=data.yr$diff_out,
                                              R = data.yr$dist2cutoff, 
                                              cutoff = 0,
                                              plot=T,
                                              level = level,
                                              covariates = as.matrix(data.yr[,covariates]),quietly = T, ...), 
                                              error = function(e) {
                                                      err <<- TRUE
                                                    })
      if (err) {
        next
      }
      
      results[row.ind, "border_seg"] <- i
      results[row.ind, "PA"] <- border.segments$PA[i]
      results[row.ind, "cov_ind"] <- border.segments$Ind[i]
      results[row.ind, "year"] <- yr
      results[row.ind, "est"] <- rd.est$obs.stat
      results[row.ind, "p_value"] <- rd.est$p.value
      results[row.ind, "se_bc"] <- round(rd.est[["se"]][2], 4)
      results[row.ind, "se_rob"] <- round(rd.est[["se"]][3], 4)
      results[row.ind, "est_l"] <- rd.est$sumstats[3,1]
      results[row.ind, "est_r"] <- rd.est$sumstats[3,2]
      results[row.ind, "window_l"] <- rd.est$window[1]
      results[row.ind, "window_r"] <- rd.est$window[2]
      results[row.ind, "segment_length"] <- border.segments$length.intersection[i]
      results[row.ind, "rel_time"] <- yr - yr.mine
      } else {
        rd.est <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, 
                                                    x = data.yr$dist2cutoff, c = 0, ...), error = function(e) {
                                                      err <<- TRUE
                                                    })
        results[row.ind, "border_seg"] <- i
        results[row.ind, "PA"] <- border.segments$PA[i]
        results[row.ind, "cov_ind"] <- border.segments$Ind[i]
        results[row.ind, "year"] <- yr
        results[row.ind, "segment_length"] <- border.segments$length.intersection[i]
        results[row.ind, "rel_time"] <- yr - yr.mine
        
        if (err) {
          next
        }
        
        results[row.ind, "est"] <- round(rd.est[["coef"]][["Conventional", 
                                                                        1]], 2)
        results[row.ind, "se_bc"] <- round(rd.est[["se"]][2], 4)
        results[row.ind, "se_rob"] <- round(rd.est[["se"]][3], 4)
        results[row.ind, "p_value"] <- round(rd.est[["pv"]][["Robust", 
                                                                   1]], 2)
        results[row.ind, "est_l"] <- round(rd.est[["tau_bc"]][[1]],2)
        results[row.ind, "est_r"] <- round(rd.est[["tau_bc"]][[1]],2)
        results[row.ind, "window_l"] <- round(rd.est[["bws"]][["h", 
                                                                     1]]/1000, 1)
        results[row.ind, "window_r"] <- round(rd.est[["bws"]][["h", 
                                                               2]]/1000, 1)
      }
      
      #c("bpoint","year","est_conv","est_rob" "se_conv", "se_rob", 
      #  "p_conv", "p_rob", "Ntr", "Nco", "bw_l", "bw_r","CI_rob_l", "CI_rob_u","est_l","est_r")
    }
  }
  return(results)  
}

######################### Diff in Dis

DiffInDisc <- function (depvar,data, doubleFronts,treatvar, calculate.windows = T, bandwidth = 10000, level = 0.15, years = 2000:2020, df.crs, covariates = NULL, locrand = F, ...) 
{
  colnames(data)[which(colnames(data)==depvar)] <- "outcome"
  names(doubleFronts)[which(names(doubleFronts)==treatvar)] <- "tg"
  groups = unique(doubleFronts$tg)
  columnames <- c("group","year","est","p_value","se_bc","se_rob","est_l","est_r","window_l","window_r","segment_length","rel_time","n_segments","n_l","n_r")
  results <- data.frame(matrix(ncol = length(columnames), nrow = length(groups)*length(years)))
  colnames(results) <- columnames
  row.ind <- 0
  cells <- data %>%
    group_by(cell) %>%
    summarise(PA = mean(PA),x = mean(x),y = mean(y))
  
  cells <- terra::vect(cells, geom=c("x", "y"),crs=df.crs)

  for (rast in covariates) {
    if (rast %in% c("x","y","objectid")) next
    covrast <- rast(here(dataInt,glue::glue("{rast}_500.tif")))
    covrast <- project(covrast,crs(cells))
    print(glue::glue("Extract values for {rast}"))
    cells[[rast]] <- extract(covrast,cells)[,2]
  }
  
  for(g in groups) {
    print(g)
    border.segments <- doubleFronts %>% tidyterra::filter(tg == g)
    border.segments <- aggregate(border.segments)
    buffed.segments <- buffer(border.segments,bandwidth)
    cells.in.buff <- cells %>%
      project(border.segments) %>%
      terra::crop(buffed.segments)
    
    cells.in.buff$dist2cutoff <- distance(cells.in.buff,border.segments)
    cells.in.buff$objectid <- nearby(cells.in.buff,doubleFronts,k=1)[,2]
    cells.in.buff <- as.data.frame(cells.in.buff)
    
    dataCrop <- data %>%
      filter(cell %in% cells.in.buff$cell,year %in% years) %>%
      left_join(cells.in.buff,by=c("cell","PA")) %>%
      tidyterra:: mutate(treated = if_else(PA != 0,1,0),
                         dist2cutoff = if_else(treated==0,dist2cutoff*-1,dist2cutoff))
    #dataCrop <- dataCrop %>% group_by(cell) %>% mutate(outcome = if_else(lag(outcome) != 0,(dplyr::lag(outcome)-outcome)/dplyr::lag(outcome),0)) %>% ungroup()
    base_out <- dataCrop$outcome[dataCrop$year==(g-1)]
    #dataCrop <- dataCrop %>% mutate(diff_out = outcome - base_out)
    
    for (yr in years) {
      row.ind <- row.ind + 1
      data.yr <- dataCrop %>% filter(year == yr) %>% mutate(diff_out = outcome - base_out)
      err <- FALSE
      if(locrand == T) {
        rd.est <- tryCatch(rdlocrand::rdrandinf(Y=data.yr$diff_out,
                                                R = data.yr$dist2cutoff, 
                                                cutoff = 0,
                                                plot=T,
                                                level = level,
                                                covariates = as.matrix(data.yr[,covariates]),quietly = T, ...), 
                           error = function(e) {
                             err <<- TRUE
                           })
        if (err) {
          next
        }
        
        results[row.ind, "group"] <- g
        results[row.ind, "n"] <- length(border.segments)
        results[row.ind, "year"] <- yr
        results[row.ind, "est"] <- rd.est$obs.stat
        results[row.ind, "p_value"] <- rd.est$p.value
        results[row.ind, "est_l"] <- rd.est$sumstats[3,1]
        results[row.ind, "est_r"] <- rd.est$sumstats[3,2]
        results[row.ind, "window_l"] <- rd.est$window[1]
        results[row.ind, "window_r"] <- rd.est$window[2]
        results[row.ind, "segment_length"] <- sum(border.segments$length.intersection)
        results[row.ind, "rel_time"] <- yr - g
      } else {
        #rd.est <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, 
         #                                     x = data.yr$dist2cutoff, c = 0, covs = as.matrix(data.yr[,covariates]), ...), error = function(e) {
        #                                        err <<- TRUE
        #                                      })
        rd.est <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, 
                                              x = data.yr$dist2cutoff, c = 0, covs = as.matrix(data.yr[,covariates]), ...), error = function(e) {
                                                err <<- TRUE
                                              })
        
        rd.est.add <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, 
                                              x = data.yr$dist2cutoff, c = 0, ...), error = function(e) {
                                                err <<- TRUE
                                              })
        results[row.ind, "group"] <- g
        results[row.ind, "year"] <- yr
        results[row.ind, "segment_length"] <- sum(border.segments$length.intersection)
        results[row.ind, "rel_time"] <- yr - g
        results[row.ind,"n_segments"] <- length(border.segments)
        
        if (err) {
          next
        }
        
        results[row.ind, "est"] <- rd.est$coef[2]
        results[row.ind, "se_bc"] <- rd.est$se[2]
        results[row.ind, "se_rob"] <- rd.est$se[3]
        results[row.ind, "p_value"] <- rd.est$pv[2]
        results[row.ind, "est_l"] <- rd.est.add$tau_bc[1]
        results[row.ind, "est_r"] <- rd.est.add$tau_bc[2]
        results[row.ind, "window_l"] <- round(rd.est[["bws"]][["h", 
                                                               1]]/1000, 1)
        results[row.ind, "window_r"] <- round(rd.est[["bws"]][["h", 
                                                               2]]/1000, 1)
        results[row.ind, "n_l"] <- rd.est$N_h[1]
        results[row.ind, "n_r"] <- rd.est$N_h[2]
      }
  

    }
  }
  return(results)  
}

## Parametric specifications ------------------------------------------------

ParamRD <- function()

## Plot binned averages over threshold---------------------------------------

PlotBins <-function(df,covar,boxplot = FALSE) {
  
  plots <- list()
  plotpos <- 1
  for (var in covar) {
  
  tags <- c("[-5,-4)","[-4,-3)", "[-3,-2)", "[-2,-1)", "[-1,0)", "[0,1)","[1,2)", "[2,3)","[3,4)", "[4,5)")
  to_plot <- df %>% select(dist2cutof ,all_of(var)) %>%
    mutate(dist2cutof = dist2cutof*200000)
  to_plot <- as_tibble(to_plot) %>% 
    mutate(tag = case_when(
      dist2cutof  >= -5000 &dist2cutof  < -4000 ~ tags[1],
      dist2cutof  >= -4000 & dist2cutof  < -3000 ~ tags[2],
      dist2cutof  >= -3000 & dist2cutof  < -2000 ~ tags[3],
      dist2cutof  >= -2000 & dist2cutof  < -1000 ~ tags[4],
      dist2cutof  >= -1000 & dist2cutof  < 0 ~ tags[5],
      dist2cutof  >= 0 & dist2cutof  < 1000 ~ tags[6],
      dist2cutof  >= 1000 & dist2cutof  < 2000 ~ tags[7],
      dist2cutof  >= 2000 & dist2cutof  < 3000 ~ tags[8],
      dist2cutof  >= 3000 & dist2cutof  < 4000 ~ tags[9],
      dist2cutof  >= 4000 & dist2cutof  <= 5000 ~ tags[10],
    ))
  colnames(to_plot)<-c("distance","var","tag")
  
  #Generate mean and standard errors
  # test_m <- do.call(rbind, lapply(split(test, test$tag), function(d) {
  #   data.frame(mean = mean(d$var), sd = sd(d$var), tag = d$tag)
  # }))
  # #https://ggplot2.tidyverse.org/reference/ggplot.html
  # 
  # test<-test[!is.na(test$tag),]
  # ggplot(data = test, mapping = aes(x=tag,y=var)) + 
  #   geom_boxplot(fill="bisque",color="black",alpha=0.3,outlier.shape = NA) + 
  #   #scale_y_continuous(limits = c(0,0.03)) +
  #   labs(x='distance') +
  #   geom_point(data = test_m, aes(tag, mean), colour = 'red', size = 3) +
  #   guides(color="none") +
  #   theme_minimal() +
  #   coord_cartesian(ylim = quantile(test$var, c(0.1, 0.9)))
  # 

  to_plot<-to_plot[!is.na(to_plot$tag),]
  to_plot$tag <- factor(to_plot$tag, c("[-5,-4)","[-4,-3)", "[-3,-2)", "[-2,-1)", "[-1,0)", "[0,1)","[1,2)", "[2,3)","[3,4)", "[4,5]"))
  #vline location:
  myLoc <- (which(levels(to_plot$tag) == "[-1,0)") + which(levels(to_plot$tag) == "[0,1)")) / 2
  
  
  plots[[plotpos]]<-ggplot(to_plot, aes(x=tag, y=var, fill=tag)) +
    xlab("Distance to PA (km)") + ylab("Forest cover") +
    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = myLoc,linetype = "dotted") +
    theme(legend.position="none") +
    scale_fill_brewer(palette="RdBu") +
    coord_cartesian(ylim = quantile(to_plot$var, c(0.1, 0.9),na.rm=T))

  plotpos <- plotpos+1
  }
  
  return(plots)
}


## Function to plot estimates from spatial RD and confidence intervals by border point ####

PlotGRD <-function(results,PA,mines,cover = FALSE) {
#results:     Results sf object from SpatialRDmod command with 
#mines:       SF object with mining concessions
#PA:          Number of PA (PAind value)
#cover:       Plot line for cover left and right of border point

  mat = st_intersects(results, mines, sparse = FALSE)
  #mergedres_def$mine = st_covered_by(mergedres_def,GFWmines,sparse = T)
  results$mine<-as.integer(apply(mat, 1, any))
  
  if (cover == T) {
  GRDplot<-results[results$PAind==PA,] %>% mutate(Color = ifelse(mine > 0, "green", "red")) %>%
    ggplot(aes(x=Point,y=Estimate, color=Color)) +
    geom_point(shape=20, size=3) +
    geom_errorbar(aes(ymin = CI_Rob_l,ymax = CI_Rob_u)) +
    geom_hline(yintercept = 0,linetype = "dotted") + 
    geom_line(aes(x = Point, y = cover_l,color="blue"))+
    geom_line(aes(x = Point, y = cover_r,color="black"))+
    scale_color_identity(labels = c("Forest cover outside","Forest cover inside","In concession","Outside concession"),guide = "legend")
    
  } else {
    GRDplot<-results[results$PAind==PA,] %>% mutate(Color = ifelse(mine > 0, "green", "red")) %>%
      ggplot(aes(x=Point,y=Estimate, color=Color)) +
      geom_point(shape=20, size=3) +
      geom_errorbar(aes(ymin = CI_Rob_l,ymax = CI_Rob_u)) +
      geom_hline(yintercept = 0,linetype = "dotted")
  }  
  return(GRDplot)
}

## Map estimates by PA ####
MapGRD <-function(results,PA,cover = FALSE) {
  #Results:     Result dataframe with estimates form SpatialRD command
  #PA:          PA to be mapped
GRDmap<-ggplot(data = PAs[PA,]) +
  geom_sf() +
  geom_sf(data = ConcessionsIn, fill="red") +
  geom_sf_label(data = results[results$PAind==PA,],aes(label = Estimate),size=2)

return(GRDmap)
}

## Plot protection mechanism ####

## Plot protection mechanism ####
plotSankey <- function(Def = r1aDef,Cov = r1cCov, metrics = "jenks",years = c(2002,2007,2012,2017,2022)) {
  #  years = c(2002,2007,2012,2017,2022)
  #  metrics = "hosonuma"
  resDF <- Def %>% 
    rename(left_def = est_l, right_def = est_r, disc_def = est_rob) %>%
    filter(!is.na(bpoint),on.border!=1,outside.tmf!=1)
  
  resDF <- Cov %>% 
    select(bpoint, year,left_cov = est_l, right_cov = est_r, disc_cov = est_rob) %>% 
    right_join(resDF, by = c("bpoint","year")) #%>%
  #mutate(no_forest = if_else(year==2000 & left_cov <0.1 & right_cov < 0.1,1,0))
  
  #resDF <- resDF %>% filter(!(bpoint %in% resDF$bpoint[resDF$no_forest==1]))
  
  y <- 2000
  df <- data.frame()
  # Metrics "own"
  for (year in years) {
    first_y = y +1
    y = year
    if (metrics == "jenks") {
      
      df <- resDF %>% arrange(bpoint,year) %>% 
        filter(!is.na(right_cov*right_def*left_cov*left_def),
               year %in% c(first_y:y)) %>%
        group_by(bpoint) %>%
        summarise(left_cov = min(left_cov),
                  right_cov = min(right_cov),
                  left_def = mean(left_def),
                  right_def = mean(right_def),
                  year = y
        ) %>%
        ungroup() %>% # Jenks for year 2000 (average left + right)
mutate(mechanism_cov = case_when( #from most specific to most general
  right_def <=.005 & left_def >=.005 & right_cov > 0.335 ~ "Contained",
  right_def >= 0.005 ~ "Sprawling",
  right_cov <=0.335 ~ "Exhausted",
  left_cov > .79 & right_cov > .79 & left_def <= 0.005 & right_def <= 0.005 ~ "Dormant",
  right_cov > 0.335 & right_def < 0.005 ~ "Consolidated"
)) %>%
  bind_rows(df)
    }
    
    
    if (metrics == "desy") {
      
      df <- resDF %>% arrange(bpoint,year) %>% 
        filter(!is.na(right_cov*right_def*left_cov*left_def),
               year %in% c(first_y:y)) %>%
        group_by(bpoint) %>%
        summarise(left_cov = min(left_cov),
                  right_cov = min(right_cov),
                  left_def = mean(left_def),
                  right_def = mean(right_def),
                  year = y
        ) %>%
        ungroup() %>% 
        mutate(mechanism_cov = case_when( #from most specific to most general
          right_def <=.0037 & left_def >=.0037 & right_cov > 0.15 ~ "Contained",
          right_def >= 0.0037 ~ "Sprawling",
          right_cov <=0.15 ~ "Exhausted",
          left_cov > .5 & right_cov > .5 & left_def <= 0.0037 & right_def <= 0.0037 ~ "Dormant",
          right_cov > 0.015 & right_def < 0.0037 ~ "Consolidated"
        )) %>%
        bind_rows(df)
    }
    
    # Metrics "Buchadas et al."
    if (metrics == "buchadas") {
      df <- resDF %>% arrange(bpoint,year) %>% 
        filter(!is.na(right_cov*right_def*left_cov*left_def),
               year %in% c(first_y:y)) %>%
        group_by(bpoint) %>%
        summarise(left_cov = min(left_cov),
                  right_cov = min(right_cov),
                  left_def = mean(left_def),
                  right_def = mean(right_def),
                  year = y
        ) %>%
        ungroup() %>%
        mutate(mechanism_cov = case_when( #from most specific to most general
          right_def <=.006 & left_def >.006 & right_cov > 0.1 ~ "Contained",
          right_def > 0.006 ~ "Sprawling",
          right_cov <=0.1 ~ "Exhausted",
          left_cov > .55 & right_cov > .55 & left_def <= 0.006 & right_def <= 0.006 ~ "Dormant",
          right_cov > 0.01 & right_def < 0.006 ~ "Consolidated"
        )) %>%
        bind_rows(df)
    }
    
    # Metrics "Hosonuma"
    if(metrics =="hosonuma") {
      df <- resDF %>% arrange(bpoint,year) %>% 
        filter(!is.na(right_cov*right_def*left_cov*left_def),
               year %in% c(first_y:y)) %>%
        group_by(bpoint) %>%
        summarise(left_cov = min(left_cov),
                  right_cov = min(right_cov),
                  left_def = mean(left_def),
                  right_def = mean(right_def),
                  year = y
        ) %>%
        ungroup() %>%
        mutate(mechanism_cov = case_when( #from most specific to most general
          right_cov > 0.5 & right_def/right_cov <.0025 & left_def/left_cov >0.0025 ~ "Contained",
          left_cov > .5 & right_cov > .5 & left_def/left_cov < 0.0025 & right_def/right_cov < 0.0025 ~ "Dormant",
          right_def/right_cov > 0.0025 ~ "Sprawling",
          right_cov <0.5 & right_def/right_cov <0.0025 ~ "Exhausted",
        )) %>%
        bind_rows(df)
      
    }
  }
  
  plot <- df %>% arrange(bpoint,year) %>% pivot_wider(id_cols = "bpoint",
                                                      names_from = "year",
                                                      names_prefix = "",
                                                      values_from = "mechanism_cov") %>%
    drop_na() %>%
    ggsankey::make_long(starts_with("20")) %>% 
    filter(!(is.na(node))) %>% 
    ggplot(aes(x=x,
               node=node,
               next_x = next_x,
               next_node =  next_node,
               fill = node,
               label = node)) +
    ggsankey::geom_sankey(flow.alpha = .6,
                          node.color = "gray30") + 
    ggsankey::geom_sankey_label(size = 3, color = "white", fill = "gray30",alpha = 0.6) +
    scale_fill_viridis_d(option = "viridis",drop = FALSE) +
    ggsankey::theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5))
  
  return(plot)
}

###################
# 
# Dbf <- DbfMin %>% filter(exploitat0 == 1)
# plots <- list()
# p = 1
# for (b in unique(Dbf$border_seg)) {
#   plots[[p]] <- r2.1LocRandCov.Mine5000.tmf %>% filter(border_seg == b) %>%
#     mutate(est = if_else(rel_time==-1,0,est)) %>%
#     ggplot(aes(x=rel_time,y=est)) +
#     geom_point() +
#     geom_line() +
#     geom_errorbar(aes(ymin = est - se_rob*1.96,ymax = est + se_rob*1.96))
#   p = p+1
# }
# plots[[1]]
# OLD STUFF ------------------------------------------------------------------


# Generate distance raster-################################

DistanceRaster <- function(bandwidth = 5000,PAs = s2PAs,TMFrast,PAcount,PADF) {
  
  #Function calculates a raster with distance to PA as cell values for a given bandwidth around PAs
  #bandwidth:   Radius in m for distance calculation around PA (determines observations for non-parametric estimation)
  #PAcount:     list with row numbers referring to the PAs in dataframe to be used for analysis
  #TMFrast:     Raster layer with the same extent and crs as the deforestation raster. Will serve as a template to write distances
  #PAs:         Vector with polygons of PAs
  #PADF:        Dataframe with 
  
  for (x in PAcount) {
    TMFdistTemp <-distToPA(TMFrast,PAs,PADF,bandwidth,x)
    if (exists("TMFdist")) {
      TMFdist<-mosaic(TMFdistTemp[[1]],TMFdist,fun="min")
      TMFPAs<-extend(TMFPAs,TMFdist)
      TMFyr<-extend(TMFyr,TMFdist)
      TMFdistTemp<-extend(TMFdistTemp,TMFdist)
      #TMFPAs<-cover(TMFPAs,TMFdistTemp[[2]])
      TMFPAs<-ifel((TMFdist>TMFdistTemp[[1]])|is.na(TMFPAs),TMFdistTemp[[2]],TMFPAs)
      TMFyr<-ifel((TMFdist>TMFdistTemp[[1]])|is.na(TMFyr),TMFdistTemp[[3]],TMFyr)
      TMFPAs<-extend(TMFPAs,TMFdist)
    } else {
      TMFdist<-TMFdistTemp[[1]]
      TMFPAs<-TMFdistTemp[[2]]
      TMFyr<-TMFdistTemp[[3]]
    }
  }
  
  TMFdist<-c(TMFdist,TMFPAs,TMFyr)
  return(TMFdist)
}

#Add TMF data
# DefYear<-list.files(path="Data/TMF data/Deforestation", pattern='TMF', full.names = TRUE)
# for (tile in DefYear) {
#   temptile<-rast(tile)
#   if (exists("TMFDef")) {
#   TMFDef<-merge(TMFDef,temptile)
#   } else {
#     TMFDef<-temptile
#   }
# }
# writeRaster(TMFDef,"Data/TMF data/Deforestation/TMFDef.tif",overwrite=T)
### Repeat with Degradation and Forest Cover! ###

#Merge DF with TMF data for each year
# TMFyears<-list.files(path="Data/TMF data/annual", pattern='TMF', full.names = TRUE)
# for (x in TMFyears) {
#   TMFx<-rast(x)
#   TMFxDF<-as.data.frame(TMFx,cells=T)
#   colnames(TMFxDF)<-c("cell",paste("cover",x,sep=""),paste("degr",x,sep=""),paste("def",x,sep=""),paste("reg",x,sep=""))
#   mergedDF<-merge.data.frame(mergedDF,TMFxDF,all.x=T,by="cell")
# }
# save(mergedDF,file="Data/R Outputs/mergedTS_FA.RData",overwrite=T)
# 
# BaseColNr<-ncol(mergedDF)
# 
# mergedDF$defBef<-NA
# mergedDF$defAft<-NA
# 
# for (x in PAcount) {
#   year<-mergedDF$STATUS_YR[mergedDF$PAind==x][[1]]
#   mergedDF$defBef[mergedDF$PAind==x]<-NA
#   mergedDF$defAft[mergedDF$PAind==x]<-NA
#   


#   #Disaggregate data by outcome
#   #Column names
#   newnames<-c(2020:2000)
#   outcomes<-c("cover","degr","def","reg")
#   
#   fcover<-mergedDF[,grep(cover, names(mergedDF), value=TRUE)]
#   colnames(fcover)<-newnames
#   
#   substrRight <- function(x, n){
#     substr(x, nchar(x)-n+1, nchar(x))
#   }
#   fcover<-cbind(mergedDF[,BaseColNr],fcover)
# }
# save(cover,file="Data/R Outputs/fcover.RData",overwrite=T)
# save(degr,file="Data/R Outputs/fdegr.RData",overwrite=T)
# save(def,file="Data/R Outputs/fdef.RData",overwrite=T)
# save(reg,file="Data/R Outputs/freg.RData",overwrite=T)



#To be turned into function:
#Routine to plot border points over years:
# pd <- position_dodge(0.5)
# #resPre[resPre$PAind==PA,] %>% mutate(Color = ifelse(mine > 0, "green", "red")) %>%
#   ggplot(resPre[resPre$PAind==PA,],aes(x=Defy,y=Estimate, color=Point, group=Point)) +
#   geom_point(shape=20, size=3,position=pd) +
#   geom_errorbar(aes(ymin = CI_Rob_l,ymax = CI_Rob_u), width=.2, position=pd) +
#   geom_hline(yintercept = 0,linetype = "dotted")