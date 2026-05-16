
# 0 Load packages and directories ----
library(pacman)
p_load("tidyverse",
       "here",
       "terra",
       "tidyterra",
       "sf",
       "stars",
       install=T)

projectfolder<-here::here()
dataPrep <- here(projectfolder,"dataPrep")
dataInt <- here(projectfolder,"dataInt")
dataOut <- here(projectfolder,"dataOut")
figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
tables <- here(projectfolder,"tables")

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
load(file = here(dataInt,"robust_spillover_fcov.Rdata"))

library(fixest)
reg1 <- feols(def ~ ring_neg4 + ring_2 + ring_4 + ring_6 + ring_8 + ring_10 | year + PA,data=reg.dat, cluster = "PA")
reg2 <- feols(def ~ ring_neg4 + ring_2 + ring_4 + ring_6 + ring_8 + ring_10 + accessibility + altitude + agr_suitability + slope | year + PA ,data=reg.dat, cluster = "PA")
summary(reg1)
summary(reg2)

# library(splm)
# spml(undisturbed ~ x + y + I(x*y) + ring_neg4 + ring_2 + ring_4 + ring_6 + ring_8 + ring_10 | year, 
#      data = reg.dat, 
#      model = "within",
#      effect = "twoways",
#      lag = T)


# Endogenous timing ------------------------------------------------------------

# 
# load(file= here(data,"bpoints.tmf.clean.Rdata"))
# PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
# bpoints.df <- as.data.frame(PAs.v) %>%
#   select(PA,name = nom_orig,PA.year = year,name_full = nom_ap,province, cat_uicn,type_gouv) %>%
#   right_join(bpoints.df) 
# 
# weights <- data.frame(bpoints.df[,10])
# colnames(weights) <- "bpoint"
# # Define weighting vectors for each IUCN category
# 
# for (i in unique(bpoints.df$PA)) {
#   weights[bpoints.df$PA == i,ncol(weights)+1] <- 1/nrow(bpoints.df[bpoints.df$PA == i,])
#   colnames(weights)[ncol(weights)] <- glue::glue("PA{i}")
# }
# weights[is.na(weights)] <- 0
# 
# # Aggregate boundary points by IUCN category
# 
# AATE.PA.cov <- data.frame()
# AATE.PA.def <- data.frame()
# 
# for(i in 1:length(unique(bpoints.df$PA))) {
#   for (y in 1:25) {
#     #Fcover
#     yr <- y+1999
#     fcov_RDs[[y]]
#     AATE.cov.y <- AATE_bpt_est(x=fcov_RDs[[y]],AATE = weights[,i+1])
#     AATE.cov.y <- AATE.cov.y %>% mutate(year = yr,
#                                         PA = unique(bpoints.df$PA)[i])
#     
#     AATE.PA.cov <- rbind(AATE.PA.cov,AATE.cov.y)
#     
#     #Deforestation
#     yr <- y+1999
#     def_RDs[[y]]
#     AATE.def.y <- AATE_bpt_est(x=def_RDs[[y]],AATE = weights[,i+1])
#     AATE.def.y <- AATE.def.y %>% mutate(year = yr,
#                                         PA = unique(bpoints.df$PA)[i])
#     
#     AATE.PA.def <- rbind(AATE.PA.def,AATE.def.y)
#   }
# }
# 
# # Plot results
# PAs.df <- as.data.frame(PAs.v) %>%
#   select(name =nom_orig,PA.year = year,iucn = cat_uicn,PA)
# AATE.PA.cov <- AATE.PA.cov %>%
#   left_join(PAs.df,relationship = "many-to-one") %>%
#   mutate(name = if_else(PA == 27,"Rubi-Tele",name),name = if_else(PA == 9,"Salonga I",name),name = if_else(PA == 22,"Salonga II",name),
#          less_than_10_bpts = if_else(PA %in% c(7,32,10,5,14,18,35,3),1,0),
#          PA.year = if_else(PA==30,1939,PA.year))
# 
# AATE.PA.def <- AATE.PA.def %>%
#   left_join(PAs.df,relationship = "many-to-one") %>%
#   mutate(name = if_else(PA == 27,"Rubi-Tele",name),name = if_else(PA == 9,"Salonga I",name),name = if_else(PA == 22,"Salonga II",name),
#          less_than_10_bpts = if_else(PA %in% c(7,32,10,5,14,18,35,3),1,0),
#          PA.year = if_else(PA==30,1939,PA.year))
# 
# #save(AATE.PA.cov,file = here(dataOut,"AATE.PA.cov.Rdata"))
# 
# figA2a <- AATE.PA.cov %>% 
#   mutate(border.PA = if_else(PA %in% c(23,34,15),1,0),
#          rel_time = year - PA.year) %>%
#   filter(rel_time %in% c(-3:-1),border.PA != 1,PA.year>=2001,name != "Kibali-Ituri") %>% #
#   filter(!(PA %in% c(7,32,10,5,14,18,35,3))) %>% # PAs with more than 10 points
#   mutate(PA = factor(PA,levels = unique(bpoints.df$PA)),
#          rel_time = factor(rel_time)) %>%
#   ggplot(aes(y=name,x=est_AATE,color = rel_time, group = rel_time)) +
#   geom_point(position = position_dodge(width = 0.5),size=.4) +
#   geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.35,height = .51,position = position_dodge(width = 0.5)) +
#   geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
#   theme_bw(base_size = 7) +
#   theme(legend.key.width=unit(.2,"cm")) +
#   scale_color_manual("rel. time",values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
#   #scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
#   #scale_shape_manual(values = c(15,16,17,18,8,3)) +
#   xlab("Discontinuity in forest cover") +
#   ylab("PA")
# 
# figA2b <- AATE.PA.def %>% 
#   mutate(border.PA = if_else(PA %in% c(23,34,15),1,0),
#          rel_time = year - PA.year) %>%
#   filter(rel_time %in% c(-3:-1),border.PA != 1,PA.year>=2001,name != "Kibali-Ituri") %>% #
#   filter(!(PA %in% c(7,32,10,5,14,18,35,3))) %>% # PAs with more than 10 points
#   mutate(PA = factor(PA,levels = unique(bpoints.df$PA)),
#          rel_time = factor(rel_time)) %>%
#   ggplot(aes(y=name,x=est_AATE,color = rel_time, group = rel_time)) +
#   geom_point(position = position_dodge(width = 0.5),size=.6) +
#   geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.6,height = .51,position = position_dodge(width = 0.5)) +
#   geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
#   theme_bw(base_size = 7) +
#   theme(legend.key.width=unit(.2,"cm")) +
#   scale_color_manual("rel. time",values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
#   #scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
#   #scale_shape_manual(values = c(15,16,17,18,8,3)) +
#   xlab("Discontinuity in deforestation") +
#   ylab("PA")
# 
# figA2a + figA2b
# 
# 
# 


  # Cohort pre-establishment  -----------------------------------------------------------
# merge with info on PAs:
load(file= here(data,"bpoints.tmf.Rdata"))
bpoints.df$bpoint <- c(1:nrow(bpoints.df))
PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
bpoints.df <- as.data.frame(PAs.v) %>%
  select(PA,name = nom_orig, PA.year=year, name_full = nom_ap,province, cat_uicn,type_gouv) %>%
  right_join(bpoints.df) %>%
  filter(!bpoint %in% drop)

weights <- data.frame(bpoints.df[,10])
colnames(weights) <- "bpoint"

# Define weighting vectors for each IUCN category

#cat1_bpoints <- bpoints.df$bpoint[bpoints.df$cat_uicn %in% c("Ia","Ib")]
y16_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2016]# & bpoints.df$PA.year<=2000]
y12_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2012]#& bpoints.df$PA.year<=2000]
y9_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2009]
y7_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2007]
y6_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2006 & !(bpoints.df$name == "Itombwe")]
y4_bpoints <- bpoints.df$bpoint[bpoints.df$PA.year == 2004]

weights <- weights %>% mutate(y16 = if_else(bpoint %in% y16_bpoints,1/length(y16_bpoints),0),
                              y12 = if_else(bpoint %in% y12_bpoints,1/length(y12_bpoints),0),
                              y9 = if_else(bpoint %in% y9_bpoints,1/length(y9_bpoints),0),
                              y7 = if_else(bpoint %in% y7_bpoints,1/length(y7_bpoints),0),
                              y6 = if_else(bpoint %in% y6_bpoints,1/length(y6_bpoints),0),
                              y4 = if_else(bpoint %in% y4_bpoints,1/length(y4_bpoints),0))

# Aggregate boundary points by IUCN category

AATE.cohort.cov <- data.frame()
AATE.cohort.def <- data.frame()
AATE.cohort.degr <- data.frame()

for(i in 1:6) {
  for (y in 1:25) {
    yr <- y+1999
    
    #forest cover results
    fcov_RDs[[y]]
    AATE.cov.y <- AATE_bpt_est(x=fcov_RDs[[y]],AATE = weights[,i+1])
    AATE.cov.y <- AATE.cov.y %>% mutate(year = yr,
                                        y = case_when(
                                          i == 1 ~ 2016,
                                          i == 2 ~ 2012,
                                          i == 3 ~ 2009,
                                          i == 4 ~ 2007,
                                          i == 5 ~ 2006,
                                          i == 6 ~ 2004
                                        ))
    
    AATE.cohort.cov <- rbind(AATE.cohort.cov,AATE.cov.y)
    
    #deforestation results
    AATE.def.y <- AATE_bpt_est(x=def_RDs[[y]],AATE = weights[,i+1])
    AATE.def.y <- AATE.def.y %>% mutate(year = yr,
                                        y = case_when(
                                          i == 1 ~ 2016,
                                          i == 2 ~ 2012,
                                          i == 3 ~ 2009,
                                          i == 4 ~ 2007,
                                          i == 5 ~ 2006,
                                          i == 6 ~ 2004
                                        ))
    
    AATE.cohort.def <- rbind(AATE.cohort.def,AATE.def.y)
    
    #degradation results
    AATE.degr.y <- AATE_bpt_est(x=degr_RDs[[y]],AATE = weights[,i+1])
    AATE.degr.y <- AATE.degr.y %>% mutate(year = yr,
                                          y = case_when(
                                            i == 1 ~ 2016,
                                            i == 2 ~ 2012,
                                            i == 3 ~ 2009,
                                            i == 4 ~ 2007,
                                            i == 5 ~ 2006,
                                            i == 6 ~ 2004
                                          ))
    
    AATE.cohort.degr <- rbind(AATE.cohort.degr,AATE.degr.y)
  }
}

figA2a <- AATE.cohort.cov %>% mutate(rel.time = year-y) %>%
  filter(rel.time %in% c(-1,-5,-10)) %>% #
  mutate(cohort = factor(y),
         rel_time = factor(rel.time)) %>%
  ggplot(aes(y=cohort,x=est_AATE,color = rel_time, group = rel_time)) +
  geom_point(position = position_dodge(width = 0.5),size=.6) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.6,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
  theme_bw(base_size = 7) +
  theme(legend.key.width=unit(.2,"cm")) +
  scale_color_manual("rel. time",values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Discontinuity in deforestation") +
  ylab("Cohort of PA establishment")

figA2b <- AATE.cohort.def %>% mutate(rel.time = year-y) %>%
  filter(rel.time %in% c(-3:-1)) %>% #
  mutate(cohort = factor(y),
         rel_time = factor(rel.time)) %>%
  ggplot(aes(y=cohort,x=est_AATE,color = rel_time, group = rel_time)) +
  geom_point(position = position_dodge(width = 0.5),size=.6) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.6,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
  theme_bw(base_size = 7) +
  theme(legend.key.width=unit(.2,"cm")) +
  scale_color_manual("rel. time",values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Discontinuity in deforestation") +
  ylab("Cohort of PA establishment")

figA2c <- AATE.cohort.degr %>% mutate(rel.time = year-y) %>%
  filter(rel.time %in% c(-3:-1)) %>% #
  mutate(cohort = factor(y),
         rel_time = factor(rel.time)) %>%
  ggplot(aes(y=cohort,x=est_AATE,color = rel_time, group = rel_time)) +
  geom_point(position = position_dodge(width = 0.5),size=.6) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.6,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
  theme_bw(base_size = 7) +
  theme(legend.key.width=unit(.2,"cm")) +
  scale_color_manual("rel. time",values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  #scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Discontinuity in deforestation") +
  ylab("Cohort of PA establishment")

figA2 <- figA2a + figA2b + figA2c + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")
ggsave(here(figures,"figA2.pdf"),plot=figA2,
       width = 140, height = 120, units = "mm", dpi = 300)