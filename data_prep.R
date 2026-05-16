#######################
# Data preparation
# Note: due to data storage limitations, source data to run this script will be provided on request.
#######################



load(file = here(dataInt,"def.y.df.Rdata"))
load(file = here(scrap,"def.y.df.s.Rdata"))

s2PAs <- vect(here(dataPrep,"s2PAs_joined.shp"))
def.y <- rast(here(dataPrep,"undisturbed",glue::glue("TMF_undisturbed_DRCDec1995.tif")))
s2PAs <- project(s2PAs,def.y)
treatRast <- rasterize(s2PAs,def.y,field = 1, background = 0)
treat.df <- as.data.frame(treatRast,xy=T,cells=F)
names(treat.df)[[3]] <- "treat"

undist.panel <- data.frame()

for (y in c(1990:1999)) {
  undist.y <- rast(here(dataPrep,"undisturbed",glue::glue("TMF_undisturbed_DRCDec{y}.tif")))
  undist.y.df <- as.data.frame(undist.y,xy=T,cells = F)
  colnames(undist.y.df)[[3]] <- "undisturbed"
  undist.y.df$year <- y
  undist.y.df <- undist.y.df %>% left_join(treat.df)
  undist.panel <- rbind(undist.panel,undist.y.df)
  print(y)
}

save(undist.panel,file = here(dataInt,"undist.panel.Rdata"))

load(file = here(dataInt,"def.y.df.Rdata"))


bp.coords <- bpoints.df %>% mutate(x = round(x,4),y=round(y,4))
bpoints.v <- st_read(here(dataPrep,"bpoints_pts_clean.shp")) 
bpoints.v <- st_transform(bpoints.v,crs(def.y))
#  bpoints.v <- bpoints.v[!st_is_empty(bpoints.v),]
  st_write(bpoints.v,here(dataPrep,"bpoints_tmf_pts.shp"),append=F)
coords.v <- as.data.frame(st_coordinates(bpoints.v)) %>%
  mutate(revised = 1,x=round(X,4),y=round(Y,4)) %>% right_join(bpoints.df)
bpts.df <- as.data.frame(bpoints.v,geom = "XY")
st_transform(bpoints.v,)

drop2 <- c(1,5,15:20,39:43,85:92,98,99,103,107:109,111,112,166,167,240,306,442,443,444,445,446,447,448,456,457,458,459,461,462,463,464,465,577,578,581,582,583,584,670,671,740,793,794,828,829,830,831,832,833,834,850,851,859,860)

c(1,5,15,16,17,18,19,20,39:43,74:93,98,99,103,106:109,111,112,166:171,239,240,244,305,306,351,359:366,391,442:448,456:465,515:517,577:587,670:671,724,725,740,756,785:802,804,810,811,821,827:839,841,847,849:852,859,860,862)

#Required packages
# Package library ############################################################
library(pacman)
#Basic packages
p_load("tidyverse","dplyr","here","ggplot2","terra","SpatialRDD","sf","stars","terra",install=T)

# Directories ####

projectfolder<-here::here()
#source(here("code","0afunctionSheet.R"))
dataPrep <- here(projectfolder,"dataPrep")
dataOut <- here(projectfolder,"dataOut")
dataInt <- here(projectfolder,"dataInt")

figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
tables <- here(projectfolder,"tables")

################
# (1) Generate panel data
################
s2PAs <- vect(here(dataPrep,"s2PAs_joined.shp"))

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

i1cCoverPanel <- AnnualRastersToPanel(years = 2000:2022,file.pattern = "fcover",folder.directory = here(dataInt,"aggregated.fcover500"),PAs = s2PAs)


#Mining concession polygons
s3MineCon<-vect(here(dataPrep,"permis_miniers_actif.shp"))
s3aMineExpl <- vect(dataPrep)

#Concession types to keep:
#6210301 -> Autorisation d’Exploitation de Carrières Permanente
#6210100 -> Permis d'Exploitation (PE)
#6210101 -> Permis d'Exploitation de Petite Mine (PEPM)
#6210201 -> Permis d'Exploitation des Rejets (PER)

################
# (2) Covariates
################

s3aMineExpl <- vect(here(dataPrep,"permis_miniers.shp"))

s3aMineExpl <- s3aMineExpl %>% tidyterra::mutate(exploitation = if_else(type_ %in% c(6210301,6210100,6210101,6210201),1,0)) %>%
  tidyterra::mutate(year.mine.start = as.numeric(strftime(as.Date(date_2,format ="%d/%m/%Y"),"%Y")),year.mine.end = as.numeric(strftime(as.Date(date_3,format ="%d/%m/%Y"),"%Y")),year.mine.appl = as.numeric(strftime(as.Date(date_1,format ="%d/%m/%Y"),"%Y"))) %>%
  tidyterra::select(year.mine.start,year.mine.end,year.mine.appl,type = type_,desc = desc_type,name = nom_ste,ressource,status = statu_perm,exploitation,province)

writeVector(s3aMineExpl,here(dataPrep,"s3aMineClean.shp"),overwrite=T)

bpts.diff.v <- vect(here(dataPrep,"bpoints_pts_clean.shp"))
covariates <- c("s3aAccess","s3bAgrSuit","s3cSlope","s3dAltitude")

for (rast in covariates) {
  covrast <- rast(here(dataPrep,glue::glue("{rast}.tif")))
  covrast <- project(covrast,crs(bpts.diff.v))
  bpts.diff.v[[rast]] <- extract(covrast,bpts.diff.v)[,2]
  names(bpts.diff.v)[length(names(bpts.diff.v))] <- rast
}
bpts.diff <- st_as_sf(bpts.diff.v)
bpts.diff$bpoint <- 1:nrow(bpts.diff)

roads <- vect(here(dataPrep,"roads_FA.shp"))
roads <- roads %>% filter(desc_type %in% c("Route nationale","Route provinciale","Route locale")) %>% aggregate() %>% project(bpts.diff.v)
bpts.diff$road.dist.major <- distance(bpts.diff.v,roads)[,1]/1000

grip_roads <- vect(here(dataPrep,"GRIP4.shp"))
grip_roads <- project(grip_roads,bpts.diff.v)
grip_roads_agg <- aggregate(grip_roads)
bpts.diff$GRIP_road <- distance(bpts.diff.v,grip_roads_agg)[,1]/1000

loggroad_old <- vect(here(dataPrep,"logg_roads_old.shp"))
loggroad_old <- project(loggroad_old,bpts.diff.v)
loggroad_old_agg <- aggregate(loggroad_old)
bpts.diff$loggroads_old <- distance(bpts.diff.v,loggroad_old_agg)[,1]/1000

loggroad <- vect(here(dataPrep,"logg_roads.shp"))
loggroad <- project(loggroad,bpts.diff.v)
loggroad_agg <- aggregate(loggroad)
bpts.diff$loggroads <-  distance(bpts.diff.v,loggroad_agg)[,1]/1000

s3conflict <- vect(here(dataPrep,"s3gConflict.shp"))
conflicts <- as.data.frame(s3conflict) %>%
  mutate(id.y = 1:nrow(s3conflict)) %>%
  rename(yr.conflict = year)

buff.bpoints <-  buffer(bpts.diff.v,7500)
buff.bpoints$bpoint <- 1:nrow(buff.bpoints)
conflicts.bpts <- extract(buff.bpoints,s3conflict) %>%
  left_join(conflicts,by="id.y") %>%
  filter(yr.conflict >= 2000,!is.na(bpoint)) %>%
  group_by(bpoint) %>%
  summarise(conflicts = n())

bpts.diff<- bpts.diff %>% 
  left_join(conflicts.bpts) %>%
  mutate(conflicts = if_else(is.na(conflicts),0,conflicts))

forestroads <- vect(here(dataPrep,"forestroads.shp"))
forestroads <- project(forestroads,bpts.diff.v)
forestroads <- aggregate(forestroads)
bpts.diff$new_forestroads <- distance(bpts.diff.v,forestroads)[,1]/1000

bpts.cov <- st_transform(bpts.diff,crs(def.y))
to_merge <- bpoints.df %>% mutate(x = round(x,4),y = round(y,4))
coords.v <- as.data.frame(st_coordinates(bpts.diff)) %>%
  mutate(wrong_bpoint = 1:nrow(bpts.diff)) %>%
  mutate(revised = 1,x=round(X,4),y=round(Y,4)) %>% right_join(to_merge)
bpts.cov <- bpts.diff %>% rename(wrong_bpoint = bpoint) %>% right_join(coords.v)
st_write(bpts.cov,here(dataInt,"bpoints_cov.shp"),append=F)
bpts.cov <- bpts.cov %>% st_drop_geometry()
bpts.cov <- bpts.cov %>% rename(PA.year = year,iucn= cat_ucn)
save(bpts.cov,file = here(dataInt,"bpts.clean.cov.Rdata"))


# Spillover/leakage ------------------------------------------------------------

# Buffer rings to estimate outcomes
buffer_rings.v <- vect(here(dataInt,"robust_buffers_0(2)10.shp"))
tmf_rast <- rast(here(dataPrep,"deforested",glue::glue("TMF_DRC_deforestedDec2024.tif")))
buffer_rings.v <- project(buffer_rings.v,tmf_rast)
buffer_rast <- rasterize(buffer_rings.v,tmf_rast,field = "distance")
buff.df <- as.data.frame(buffer_rast,xy=T,cells=F)
names(treat.df)[[3]] <- "dist"

load(file = here(dataInt,"def.panel.Rdata"))
reg.dat <- def.panel %>% right_join(buff.df) 

# Allocate cells to PA for clustered SEs
buffer10.v <- vect(here(dataInt,"robust_buffers_10.shp"))
buffer10.v <- project(buffer10.v,tmf_rast)
PA_rast <- rasterize(buffer10.v,tmf_rast,field = "PA")
PA.df <- as.data.frame(PA_rast,xy=T,cells=F)
reg.dat <- reg.dat %>% left_join(PA.df)

# Covariates to include

#Acess:
access <- rast(here(dataPrep,"accessibility_to_cities_2015_v1.0.tif"))
access <- project(access,tmf_rast)
access <- terra::resample(access,tmf_rast)
access.df <- as.data.frame(access,xy=T,cells =F)
names(access.df)[[3]] <- "accessibility"
reg.dat <- reg.dat %>% left_join(access.df)

#Altitude:
altitude <- rast(here(dataPrep,"s3dAltitude.tif"))
altitude <- project(altitude,tmf_rast)
altitude <- terra::resample(altitude,tmf_rast)
altitude.df <- as.data.frame(altitude,xy=T,cells =F)
names(altitude.df)[[3]] <- "altitude"
reg.dat <- reg.dat %>% left_join(altitude.df)

#Agricultural suitablity:
agr_suitability <- rast(here(dataPrep,"s3bAgrSuit.tif"))
agr_suitability <- project(agr_suitability,tmf_rast)
agr_suitability <- terra::resample(agr_suitability,tmf_rast)
agr_suitability.df <- as.data.frame(agr_suitability,xy=T,cells =F)
names(agr_suitability.df)[[3]] <- "agr_suitability"
reg.dat <- reg.dat %>% left_join(agr_suitability.df)

#Slope:
slope <- rast(here(dataPrep,"s3cSlope.tif"))
slope <- project(slope,tmf_rast)
slope <- terra::resample(slope,tmf_rast)
slope.df <- as.data.frame(slope,xy=T,cells =F)
names(slope.df)[[3]] <- "slope"
reg.dat <- reg.dat %>% left_join(slope.df)



reg.dat <- reg.dat %>% mutate(ring_neg4 = if_else(distance == -4000,1,0),
                              ring_neg2 = if_else(distance == -2000,1,0),
                              ring_2 = if_else(distance == 2000,1,0),
                              ring_4 = if_else(distance == 4000,1,0),
                              ring_6 = if_else(distance == 6000,1,0),
                              ring_8 = if_else(distance == 8000,1,0),
                              ring_10 = if_else(distance == 10000,1,0),
)

save(reg.dat,file = here(dataInt,"robust_spillover_def.Rdata"))
