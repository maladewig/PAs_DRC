
######################## Model Estimation #######################

# 0 Load packages and directories ----
library(pacman)
p_load("tidyverse",
       "here",
       "terra",
       "tidyterra",
       "sf",
       "stars",
       "kableExtra",
       "rd2d",
       "patchwork",
       "ggsankey",
       "viridis",
       install=T)

projectfolder<-here::here()
data <- here(projectfolder,"data")
figures <- here(projectfolder,"figures")
tables <- here(projectfolder,"tables")

load(file=here(data,"bpoints.tmf.clean.Rdata"))

PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
bpoints.df <- as.data.frame(PAs.v) %>%
  select(PA,name = nom_orig,PA.year = year,name_full = nom_ap,province, cat_uicn,type_gouv) %>%
  right_join(bpoints.df)

AATE_bpt_est <- function(x,AATE = NULL) {
  
  if (is.null(AATE)){
    AATE <- rep(1/nrow(x$results),nrow(x$results))
  } else {
    AATE <- AATE / sum(AATE)
  }
  
  ATEdf <- data.frame(sum(AATE * x$results$Est.q))
  colnames(ATEdf) <- "est_AATE"
  ATEdf$se_AATE <- sqrt(AATE %*% x$cov.q %*% AATE)[1,1]
  ATEdf$zvalue_AATE <- sum(AATE * x$results$Est.q)/ATEdf$se_AATE
  ATEdf$pvalue_AATE <- 2 * pnorm(abs(ATEdf$zvalue_AATE),lower.tail = FALSE)
  
  if (x$opt$side == "two"){
    zval <- qnorm((x$opt$level + 100)/ 200)
    ATEdf$CI.lower_AATE <- sum(AATE * x$results$Est.q) - zval * ATEdf$se_AATE
    ATEdf$CI.upper_AATE <- sum(AATE * x$results$Est.q) + zval * ATEdf$se_AATE
  }
  if (x$opt$side == "left"){
    zval <- qnorm(x$opt$level / 100)
    ATEdf$CI.upper_AATE <- x$results$Est.p + zval * ATEdf$se_AATE
    ATEdf$CI.lower_AATE <- rep(-Inf, length(CI.upper_AATE))
  }
  if (x$opt$side == "right"){
    zval <- qnorm(x$opt$level / 100)
    ATEdf$CI.lower_AATE <- x$results$Est.p - zval * ATEdf$se_AATE
    ATEdf$CI.upper_AATE <- rep(Inf, length(CI.lower_AATE))
  }
  
  return(ATEdf)
}

## Load RD estimates (computed on HPC) -----------------------------------------
def_RDs <-readRDS(here(data,"est_def_clean.Rdata"))
fcov_RDs <-readRDS(here(data,"est_fcov_clean.Rdata"))
degr_RDs <- readRDS(here(data,"est_degr_clean.Rdata"))

## Results by bpoint -----------------------------------------------------------
res.list <- fcov_RDs
r1.fcov.bpoints.df <- data.frame()

for (i in 01:25) {
year <- i+1999
r1.bpoints.df <- as.data.frame(fcov_RDs[[i]]$results$Est.q)
colnames(r1.bpoints.df) <- c("est")
r1.bpoints.df$p_value <- round(fcov_RDs[[i]]$results$`P>|z|`,3)
r1.bpoints.df$se <- fcov_RDs[[i]]$results$Se.q 
r1.bpoints.df$CI_l <- fcov_RDs[[i]]$results$CI.lower
r1.bpoints.df$CI_u <- fcov_RDs[[i]]$results$CI.upper
r1.bpoints.df$CB_l <- fcov_RDs[[i]]$results$CB.lower
r1.bpoints.df$CB_u <- fcov_RDs[[i]]$results$CB.upper
r1.bpoints.df$Nh_0 <- fcov_RDs[[i]]$results$Nh0
r1.bpoints.df$Nh_1 <- fcov_RDs[[i]]$results$Nh1
r1.bpoints.df$est_l <- fcov_RDs[[i]]$results.A0$Est.q
r1.bpoints.df$est_r <- fcov_RDs[[i]]$results.A1$Est.q
r1.bpoints.df$se_l <- fcov_RDs[[i]]$results.A0$Se.q
r1.bpoints.df$se_r <- fcov_RDs[[i]]$results.A0$Se.q
r1.bpoints.df$year <- year
r1.bpoints.df$bpoint <- bpoints.df$bpoint
r1.fcov.bpoints.df <- rbind(r1.fcov.bpoints.df,r1.bpoints.df)
}

res.list <- def_RDs

r2.def.bpoints.df <- data.frame()
for (i in 01:25) {
  year <- i+1999
  r1.bpoints.df <- as.data.frame(def_RDs[[i]]$results$Est.q)
  colnames(r1.bpoints.df) <- c("est")
  r1.bpoints.df$p_value <- round(def_RDs[[i]]$results$`P>|z|`,3)
  r1.bpoints.df$se <- def_RDs[[i]]$results$Se.q 
  r1.bpoints.df$CI_l <- def_RDs[[i]]$results$CI.lower
  r1.bpoints.df$CI_u <- def_RDs[[i]]$results$CI.upper
  r1.bpoints.df$CB_l <- def_RDs[[i]]$results$CB.lower
  r1.bpoints.df$CB_u <- def_RDs[[i]]$results$CB.upper
  r1.bpoints.df$Nh_0 <- def_RDs[[i]]$results$Nh0
  r1.bpoints.df$Nh_1 <- def_RDs[[i]]$results$Nh1
  r1.bpoints.df$est_l <- def_RDs[[i]]$results.A0$Est.q
  r1.bpoints.df$est_r <- def_RDs[[i]]$results.A1$Est.q
  r1.bpoints.df$se_l <- def_RDs[[i]]$results.A0$Se.q
  r1.bpoints.df$se_r <- def_RDs[[i]]$results.A0$Se.q
  r1.bpoints.df$year <- year
  r1.bpoints.df$bpoint <- bpoints.df$bpoint
  r2.def.bpoints.df <- rbind(r1.bpoints.df,r2.def.bpoints.df)
}

## Drop invalid points ---------------------------------------------------------

# (i) Points on national boundaries
admin1.v <- vect(here(data,"s0admin1.shp"))
bpoints.v <- vect(here(data,"bpoints_tmf_pts.shp"))
bpoints.v <- project(bpoints.v,admin1.v)
admin1_shift <- buffer(admin1.v,-5000) # to exclude points within 5km from border
nat.border <- relate(bpoints.v,admin1_shift,relation="within")

# Omit 36 points on nat. borders:
#bpoints.df <- bpoints.df[nat.border,]

# (ii) Points without initial forest cover
 fcov_90RDs <- readRDS(here(data,"est_fcov90-93.Rdata"))
 View(as.data.frame(fcov_90RDs[[1]]$results.A1$Est.p)) > 0.1

 # Kibali-Ituri (not yet established)
 kibali_ituri <- bpoints.df$bpoint[bpoints.df$name == "Kibali-Ituri"]
 kibali_ituri <- which(bpoints.df$bpoint %in% kibali_ituri)

# Overall AATE -------------------------------------------------------------------
drop_estab_aft_2010 <- c(61:69, 379:434)
drop_estab_aft_2007 <- which(bpoints.df$PA.year>2007)

AATE.cov.est <- data.frame()
AATE.def.est <- data.frame()
AATE.degr.est <- data.frame()

drop_bpoints <- unique(c(drop_estab_aft_2007,kibali_ituri))
for (yr in 1:25) {
  
  y <- yr+1999
  
## Forest cover ####
  
AATEcov = c(rep(1/(nrow(fcov_RDs[[yr]]$results[-drop_bpoints,])),nrow(fcov_RDs[[yr]]$results)))
AATEcov[drop_bpoints] <- 0
AATE.cov.y <- AATE_bpt_est(x=fcov_RDs[[yr]],AATE = AATEcov)
AATE.cov.y$year <- y
AATE.cov.y$est_l <- sum(fcov_RDs[[yr]]$results.A0$Est.q * AATEcov)
AATE.cov.y$est_r <- sum(fcov_RDs[[yr]]$results.A1$Est.q * AATEcov)
AATE.cov.est <- rbind(AATE.cov.est,AATE.cov.y)

## Deforestation
AATEdef = c(rep(1/(nrow(def_RDs[[yr]]$results[-drop_bpoints,])),nrow(def_RDs[[yr]]$results)))
AATEdef[drop_bpoints] <- 0
AATE.def.y <- AATE_bpt_est(x=def_RDs[[yr]],AATE = AATEdef)
AATE.def.y$year <- y
AATE.def.y$est_l <- sum(def_RDs[[yr]]$results.A0$Est.q* AATEdef)
AATE.def.y$est_r <- sum(def_RDs[[yr]]$results.A1$Est.q * AATEdef)
AATE.def.est <- rbind(AATE.def.est,AATE.def.y)

## Forest degradation
AATEdegr = c(rep(1/(nrow(degr_RDs[[yr]]$results[-drop_bpoints,])),nrow(degr_RDs[[yr]]$results)))
AATEdegr[drop_bpoints] <- 0
AATE.degr.y <- AATE_bpt_est(x=degr_RDs[[yr]],AATE = AATEdegr)
AATE.degr.y$year <- y
AATE.degr.y$est_l <- sum(degr_RDs[[yr]]$results.A0$Est.q * AATEdegr)
AATE.degr.y$est_r <- sum(degr_RDs[[yr]]$results.A1$Est.q * AATEdegr)
AATE.degr.est <- rbind(AATE.degr.est,AATE.degr.y)
}
AATE.def.est$disturbance <- "deforestation"
AATE.degr.est$disturbance <- "degradation"

fig3a <- AATE.cov.est %>%
  filter(year>=2007) %>%
  ggplot(aes(x=year, y=est_AATE)) + 
  geom_line(position=position_dodge(width=0.5),linewidth=.6,color = "#c7842a") +
  geom_point(position=position_dodge(width=0.5),color = "#c7842a",size=.8) +
  geom_ribbon(aes(ymin=CI.lower_AATE, ymax=CI.upper_AATE),alpha=0.18,position=position_dodge(width=0.5),color = "#c7842a",fill = "#c7842a",linewidth=.2) +
  theme_bw(base_size = 7) +
  scale_color_manual(values = c("#c7842a","#d8c656","#7d6756")) +
  scale_fill_manual(values = c("#c7842a","#d8c656","#7d6756")) +
  scale_shape_manual(values = c(15,16,17,18,8)) +
  xlab("Year") +
  ylab("Discontinuity in undisturbed forest") +
  guides(fill=guide_legend(title=" "),color=guide_legend(title=" "))


fig3b <- AATE.def.est %>% bind_rows(AATE.degr.est) %>%
  filter(year>=2007) %>%
  ggplot(aes(x=year, y=est_AATE,group = disturbance,color = disturbance,shape = disturbance,fill = disturbance)) + 
  geom_line(linewidth=.6)+
  geom_point(size=.8) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",linewidth = .6) +
  geom_ribbon(aes(ymin=CI.lower_AATE, ymax=CI.upper_AATE),alpha=0.18,linewidth=.2) +#,position=position_dodge(width=0.5)) +
  theme_bw(base_size = 7) +
  scale_color_manual(values = c("#d8c656","#440154")) +
  scale_fill_manual(values = c("#d8c656","#440154")) +
  scale_shape_manual(values = c(15,16,17,18,8)) +
  xlab("Year") +
  ylab("Discontinuity in forest disturbances") +
  theme(legend.title=element_blank())

fig3 <- fig3a + fig3b + plot_annotation(tag_levels = "a")
ggsave(here(figures,"fig3.pdf"), plot = fig3,
       width = 190, height = 100, units = "mm", dpi = 300)

# PA heterogeneities---------------------------------------------------------------------
drop <- c(1,5,15,16,17,18,19,20,39:43,74:93,98,99,103,106:109,111,112,166:171,239,240,244,305,306,351,359:366,391,442:448,456:465,515:517,577:587,670:671,724,725,740,756,785:802,804,810,811,821,827:839,841,847,849:852,859,860,862)

load(file= here(data,"bpoints.tmf.Rdata"))
bpoints.df$bpoint <- c(1:nrow(bpoints.df))
PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
bpoints.df <- as.data.frame(PAs.v) %>%
  select(PA,name = nom_orig,PA.year = year,name_full = nom_ap,province, cat_uicn,type_gouv) %>%
  right_join(bpoints.df) %>%
  filter(!(bpoint %in% c(drop))) 

weights <- data.frame(bpoints.df[,10])
colnames(weights) <- "bpoint"

# Define weighting vectors for each IUCN category

for (i in unique(bpoints.df$PA)) {
  weights[bpoints.df$PA == i,ncol(weights)+1] <- 1/nrow(bpoints.df[bpoints.df$PA == i,])
  colnames(weights)[ncol(weights)] <- glue::glue("PA{i}")
}
weights[is.na(weights)] <- 0

# Aggregate boundary points by IUCN category

AATE.PA.cov <- data.frame()
AATE.PA.def <- data.frame()

for(i in 1:length(unique(bpoints.df$PA))) {
  for (y in 1:25) {
    yr <- y+1999
    AATE.cov.y <- AATE_bpt_est(x=fcov_RDs[[y]],AATE = weights[,i+1])
    AATE.cov.y <- AATE.cov.y %>% mutate(year = yr,
                                        PA = unique(bpoints.df$PA)[i])
    
    AATE.PA.cov <- rbind(AATE.PA.cov,AATE.cov.y)
  }
}

# Plot results

PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
PAs.df <- as.data.frame(PAs.v) %>%
  select(name =nom_orig,PA.year = year,iucn = cat_uicn,PA)
AATE.PA.cov <- AATE.PA.cov %>%
  left_join(PAs.df,relationship = "many-to-one") %>%
  mutate(name = if_else(PA == 27,"Rubi-Tele",name),
         name = if_else(PA == 9,"Salonga I",name),
         name = if_else(PA == 22,"Salonga II",name),
         less_than_10_bpts = if_else(PA %in% c(7,32,10,5,14,18,35,3),1,0),
         PA.year = if_else(PA==30,1939,PA.year)) %>%
  filter(!(name=="Kibali-Ituri"))

#save(AATE.PA.cov,file = here(dataOut,"AATE.PA.cov.Rdata"))

fig4a <- AATE.PA.cov %>% 
  mutate(border.PA = if_else(PA %in% c(23,34,15),1,0)) %>%
  filter(year %in% c(2024,2014,2004),border.PA != 1,PA.year<=2004,name != "Kibali-Ituri") %>% #
  filter(!(PA %in% c(7,32,10,5,14,18,35,3))) %>% # PAs with more than 10 points
  mutate(year = factor(year,levels = c(2024,2014,2004)),PA = factor(PA,levels = unique(bpoints.df$PA))) %>%
  ggplot(aes(y=name,x=est_AATE,color = year, group = year,shape=year)) +
  geom_point(position = position_dodge(width = 0.5),size=.4) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.35,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
  theme_bw(base_size = 7) +
  theme(legend.key.width=unit(.2,"cm")) +
  scale_color_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Discontinuity in forest cover") +
  ylab("PA")

# IUCN heterogeneity -----------------------------------------------------------
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

cat1_bpoints <- bpoints.df$bpoint[bpoints.df$cat_uicn %in% c("Ia","Ib",NA) & !(bpoints.df$name == "Kibali-Ituri")]# & bpoints.df$PA.year<=2000]
cat1_bpoints <- cat1_bpoints[!(cat1_bpoints %in% drop_estab_aft_2007)]
cat1_bpoints <- cat1_bpoints[!(cat1_bpoints %in% drop_estab_aft_2007)]
cat2_bpoints <- bpoints.df$bpoint[bpoints.df$cat_uicn %in% c("II")]
cat2_bpoints <- cat2_bpoints[!(cat2_bpoints %in% drop_estab_aft_2007)]
cat4_bpoints <- bpoints.df$bpoint[bpoints.df$cat_uicn %in% c("IV")]
cat4_bpoints <- cat4_bpoints[!(cat4_bpoints %in% drop_estab_aft_2007)]
cat6_bpoints <- bpoints.df$bpoint[bpoints.df$cat_uicn %in% c("VI")]
cat6_bpoints <- cat1_bpoints[!(cat6_bpoints %in% drop_estab_aft_2007)]

weights <- weights %>% mutate(iucn1 = if_else(bpoint %in% cat1_bpoints,1/length(cat1_bpoints),0),
                              iucn2 = if_else(bpoint %in% cat2_bpoints,1/length(cat2_bpoints),0),
                              iucn4 = if_else(bpoint %in% cat4_bpoints,1/length(cat4_bpoints),0),
                              iucn6 = if_else(bpoint %in% cat6_bpoints,1/length(cat6_bpoints),0),)

# Aggregate boundary points by IUCN category

AATE.iucn.cov <- data.frame()
AATE.iucn.def <- data.frame()
AATE.iucn.degr <- data.frame()

for(i in 1:4) {
  for (y in 1:25) {
    yr <- y+1999
    
    #forest cover results
    fcov_RDs[[y]]
    AATE.cov.y <- AATE_bpt_est(x=fcov_RDs[[y]],AATE = weights[,i+1])
    AATE.cov.y <- AATE.cov.y %>% mutate(year = yr,
                                        iucn = case_when(
                                          i == 1 ~ "other",
                                          i == 2 ~ "II",
                                          i == 3 ~ "IV",
                                          i == 4 ~ "VI"
                                        ))
    
    AATE.iucn.cov <- rbind(AATE.iucn.cov,AATE.cov.y)
    
    #deforestation results
    AATE.def.y <- AATE_bpt_est(x=def_RDs[[y]],AATE = weights[,i+1])
    AATE.def.y <- AATE.def.y %>% mutate(year = yr,
                                        iucn = case_when(
                                          i == 1 ~ "other",
                                          i == 2 ~ "II",
                                          i == 3 ~ "IV",
                                          i == 4 ~ "VI"
                                        ))
    
    AATE.iucn.def <- rbind(AATE.iucn.def,AATE.def.y)
    
    #degradation results
    AATE.degr.y <- AATE_bpt_est(x=degr_RDs[[y]],AATE = weights[,i+1])
    AATE.degr.y <- AATE.degr.y %>% mutate(year = yr,
                                          iucn = case_when(
                                            i == 1 ~ "other",
                                            i == 2 ~ "II",
                                            i == 3 ~ "IV",
                                            i == 4 ~ "VI"
                                          ))
    
    AATE.iucn.degr <- rbind(AATE.iucn.degr,AATE.degr.y)
  }
}

# Plot results

fig4b <- AATE.iucn.cov %>% 
  filter(year %in% c(2007:2024)) %>%
  mutate(iucn = if_else(iucn %in% c("I",NA),"other",iucn)) %>%
  mutate(iucn = factor(iucn,levels = c("other","VI","IV","II"))) %>%
  ggplot(aes(y=iucn,x=est_AATE,color = year, group = year)) +
  geom_point(position = position_dodge(width = 1),size=.4) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.35,height = .51,position = position_dodge(width = 1)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth=.2) +
  geom_hline(yintercept=1.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=3.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  theme_bw(base_size = 7) +
  scale_color_viridis(discrete = F) +
  xlab("Forest cover discontinuity") +
  ylab("IUCN category")

fig4c <- AATE.iucn.def %>% filter(year %in% c(2007:2024)) %>%
  mutate(iucn = if_else(iucn %in% c("I",NA),"other",iucn)) %>%
  mutate(iucn = factor(iucn,levels = c("other","VI","IV","II"))) %>%
  ggplot(aes(y=iucn,x=est_AATE,color = year, group = year)) +
  geom_point(position = position_dodge(width = 1),size=.4) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.35,height = .51,position = position_dodge(width = 1)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth=.2) +
  geom_hline(yintercept=1.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=3.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  theme_bw(base_size = 7) +
  theme(axis.title.y=element_blank())+
  scale_color_viridis(discrete = F) +
  scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Defor. discontinuity") +
  ylab("IUCN category") +
  xlim(c(-0.00225,0.0014))

fig4d <- AATE.iucn.degr %>% filter(year %in% c(2007:2024)) %>%
  mutate(iucn = if_else(iucn %in% c("I",NA),"other",iucn)) %>%
  mutate(iucn = factor(iucn,levels = c("other","VI","IV","II"))) %>%
  ggplot(aes(y=iucn,x=est_AATE,color = year, group = year)) +
  geom_point(position = position_dodge(width = 1),size=.4) +
  geom_errorbarh(aes(xmin=CI.lower_AATE, xmax=CI.upper_AATE),linewidth=.35,height = .51,position = position_dodge(width = 1)) +
  geom_vline(xintercept=0, linetype="solid", color = "black",linewidth = .2) +
  geom_hline(yintercept=1.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  geom_hline(yintercept=3.5, linetype="dashed", color = "darkgrey",linewidth=.2) +
  theme_bw(base_size = 7) +
  theme(axis.title.y=element_blank()) +
  scale_color_viridis(discrete = F) +
  scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Forest degr. discontinuity") +
  ylab("IUCN category") +
  xlim(c(-0.00225,0.0014))

fig4bcd <- fig4b + fig4c + fig4d + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.key.width=unit(.2,"cm")) 

fig4 <- fig4a / fig4bcd + plot_annotation(tag_levels = "a")

ggsave(here(figures,"fig4.pdf"),plot=fig4,
       width = 140, height = 140, units = "mm", dpi = 300)

# Typology -------------------------------------------------------------------

# Jenks clustering
library(classInt)
levels.fcov.df <- data.frame()
levels.def.df <- data.frame()
for (y in c(2000:2024)) {
  #Forest cover thresholds
  fcover.y <- fcov_RDs[[y-1999]]
  intervals_fcov <- levels(classify_intervals((fcover.y$results.A0$Est.q + fcover.y$results.A0$Est.q) /2,n=3,style="jenks"))[1:2]
  levels.fcov.df <- rbind(levels.fcov.df,t(data.frame(intervals_fcov)))
  
  #Deforestation thresholds
  def.y <- def_RDs[[y-1999]]
  intervals_def <- levels(classify_intervals((def.y$results.A0$Est.q + def.y$results.A0$Est.q) /2,n=3,style="jenks"))[1:2]
  levels.def.df <- rbind(levels.def.df,t(data.frame(intervals_def)))
}

thresholds.fcov <- levels.fcov.df %>% mutate(threshold1 = as.numeric(str_extract(
  string = levels.fcov.df$V1,
  pattern = "(?<=,)[\\s]*-?([0-9]+\\.?[0-9]*)"
)),
threshold2 = as.numeric(str_extract(
  string = levels.fcov.df$V2,
  pattern = "(?<=,)[\\s]*-?([0-9]+\\.?[0-9]*)"
)))

threshold1.fcov <- mean(thresholds.fcov$threshold1)
threshold2.fcov <- mean(thresholds.fcov$threshold2)

thresholds.def <- levels.def.df %>% mutate(threshold1 = as.numeric(str_extract(
  string = levels.def.df$V1,
  pattern = "(?<=,)[\\s]*-?([0-9]+\\.?[0-9]*)"
)),
threshold2 = as.numeric(str_extract(
  string = levels.def.df$V2,
  pattern = "(?<=,)[\\s]*-?([0-9]+\\.?[0-9]*)"
)))

threshold1.def <- mean(thresholds.def$threshold1)
threshold2.def <- mean(thresholds.def$threshold2)

res.topology <- r1.fcov.bpoints.df %>% 
  select(left_cov = est_l,right_cov = est_r,bpoint,year)
res.topology <- r2.def.bpoints.df %>% 
  select(left_def = est_l,right_def = est_r,bpoint,year) %>%
  inner_join(res.topology)

res.topology <- res.topology %>%
  mutate(category = case_when( # Note: Hierarchical in order (not exclusive by definition)
    right_def >= threshold1.def ~ "Sprawling", 
    right_cov <=threshold1.fcov ~ "Exhausted",
    left_def >=threshold1.def ~ "Contained", 
    right_cov > threshold2.fcov & left_cov > threshold2.fcov ~ "Dormant",
    left_def < threshold1.def & right_def < threshold1.def ~ "Consolidated")
  )

drop_estab_aft_2007 <- bpoints.df$bpoint[bpoints.df$PA.year>2007]
fig5 <- res.topology %>%
  filter(!(bpoint %in% drop_estab_aft_2007)) %>%
  pivot_wider(id_cols = "bpoint",
              names_from = "year",
              names_prefix = "",
              values_from = "category") %>%
  select(bpoint,`2007`,`2010`,`2015`,`2020`,`2024`) %>% 
  ggsankey::make_long(starts_with("20")) %>% 
  filter(!(is.na(node))) %>%
  ggplot(aes(x=x,
             node=node, 
             next_x = next_x, 
             next_node =  next_node, 
             fill = node, 
             label = node)) +
  ggsankey::geom_sankey(flow.alpha = .6,
                        node.color = "gray30",
                        linewidth=.2) + 
  ggsankey::geom_sankey_label(size = 1.5, color = "white", fill = "gray30",alpha = 0.6,linewidth=.2) +
  scale_fill_viridis_d(option = "viridis",drop = FALSE) +
  ggsankey::theme_sankey(base_size = 9) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

ggsave(here(figures,"fig5.pdf"),plot=fig4,
       width = 140, height = 70, units = "mm", dpi = 300)

save(res.topology,file=here(data,"res.typology.Rdata"))

# Boxplots with heterogeneities ------------------------------------------------

load(file=here(data,"bpts.clean.cov.Rdata"))
load(file=here(data,"res.typology.Rdata"))

fig6a <-res.topology %>%
  left_join(bpts.cov) %>%
  mutate(iucn = if_else(iucn %in% c("Ia","Ib",NA),"other",iucn)) %>%
  filter(year == 2024) %>%
  mutate(category = fct_reorder(category, loggrds)) %>%
  mutate(category = factor(category, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("other", "II", "IV", "VI"))) %>%
  ggplot(aes(y=loggrds, x=category,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5,lwd = 0.3,varwidth=T,width = 8/length(unique(bpts.cov$iucn)),outlier.size = 1) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw(base_size = 7)  +
  xlab("") +
  ylab("Dist. to forest road (km)") +
  theme(axis.text.x = element_text(angle = -25, vjust = -.15,size = 7),
        axis.title.y = element_text(size = 7))

fig6b <- res.topology %>%
  left_join(bpts.cov) %>%
  mutate(iucn = if_else(iucn %in% c("Ia","Ib",NA),"other",iucn)) %>%
  filter(year == 2024) %>%
  mutate(cnflcts = if_else(is.na(cnflcts),0,cnflcts))  %>%
  mutate(category = fct_reorder(category, cnflcts)) %>%
  mutate(category = factor(category, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("other", "II", "IV", "VI"))) %>%
  ggplot(aes(y=cnflcts, x=category,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5,lwd = 0.3,varwidth = T,width = 8/length(unique(bpts.cov$iucn)),outlier.size = 1) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw(base_size = 7)  +
  xlab("") +
  ylab("Conflict incidences (#)") +
  annotate("point",x=3.88,y=90,color = "darkgrey",size=1) +
  annotate("point",x=3.92,y=90,color = "darkgrey",size=1) +
  annotate("text",x=3.83,y=90,label = "(",color = "darkgrey",size=1.2) +
  annotate("text",x=3.97,y=90,label = ")",color = "darkgrey",size=1.2) +
  annotate("point",x=3.7,y=90,color = "darkgrey",size=1) +
  annotate("text",x=3.65,y=90,label = "(",color = "darkgrey",size=1.2) +
  annotate("text",x=3.75,y=90,label = ")",color = "darkgrey",size=1.2) +
  annotate("point",x=5.1,y=90,color = "darkgrey",size=1) +
  annotate("text",x=5.05,y=90,label = "(",color = "darkgrey",size=1.2) +
  annotate("text",x=5.15,y=90,label = ")",color = "darkgrey",size=1.2) +
  ylim(0,90) +
  theme(axis.text.x = element_text(angle = -25, vjust = -.15,size=7),
        axis.title.y = element_text(size = 7))


fig6 <- fig6a / fig6b + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(here(figures,"fig6.pdf"),plot=fig6,
       width = 90, height = 90, units = "mm", dpi = 600)

# Table with overlap between mining/logging and PA -----------------------------

typology24 <- res.topology %>% filter(year == 2024) %>% select(bpoint,category24 = category)

## Mining ####
bpoints.v <- vect(here(data,"bpoints_cov.shp"))

  ### Expl. permit ####
  mining_dbf <- vect(here(data,"s3aMineClean.shp"))
  mining_dbf <- mining_dbf %>%
    project(bpoints.v) %>% buffer(width=50)
  
  bpts_mining_expl <- bpoints.v %>%
    intersect(mining_dbf) %>%
    filter(type %in% c(6210301,6210100,6210101,6210201)) %>% 
    as.data.frame() %>%
    left_join(typology24) %>%
    count(category24) %>%
    pivot_wider(
      names_from = category24,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
    
  all_bpts_df <- data.frame()
  for(type in colnames(bpts_mining_expl)) {
    if (!(type %in% c("sumBpts","Frontier"))) {
      all_bpts_df[1,type] <- paste(round(bpts_mining_expl[[type]]/bpts_mining_expl$sumBpts*100,2),"%",sep="")
    }
  }
  all_bpts_df[1,"N"] <- bpts_mining_expl[["sumBpts"]]
  all_bpts_df[1,"Frontier"] <- "Mining concession (exploitation permit)"

  ### Research permit ####
  bpts_mining_rech <- bpoints.v %>%
    intersect(mining_dbf) %>%
    filter(type %in% c(6210200)) %>%
    as.data.frame() %>%
    left_join(typology24) %>%
    count(category24) %>%
    pivot_wider(
      names_from = category24,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  for(type in colnames(bpts_mining_rech)) {
    if (type == "sumBpts") next
    all_bpts_df[2,type] <- paste(round(bpts_mining_rech[[type]]/bpts_mining_rech$sumBpts*100,2),"%",sep="")
  }
  all_bpts_df[2,"N"] <- bpts_mining_rech[["sumBpts"]]
  all_bpts_df[2,"Frontier"] <- "Mining concession (research permit)"

## Logging ####

#logg_dbf <- vect(here(dataInt,"Heterogeneities","loggconcessions_dbf.shp"))
logg_conc <- vect(here(data,"Forest_concession_agreements.shp"))
logg_conc <- project(logg_conc,bpoints.v)

  ### Active ####
  bpts_logg_amenage <- bpoints.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "amenage")
  
  logg_amenage <- bpoints.v %>% filter(bpoint %in% bpts_logg_amenage$bpoint) %>%
    as.data.frame() %>%
    left_join(typology24) %>%
    count(category24) %>%
    pivot_wider(
      names_from = category24,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  for(type in colnames(logg_amenage)) {
    if (type == "sumBpts") next
    all_bpts_df[3,type] <- paste(round(logg_amenage[[type]]/logg_amenage$sumBpts*100,2),"%",sep="")
  }
  all_bpts_df[3,"N"] <- logg_amenage[["sumBpts"]]
  all_bpts_df[3,"Frontier"] <- "Logging concession (active)"
  
  ### Valid ####
  bpts_logg_valid <- bpoints.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "valide")
  
  logg_valid <- bpoints.v %>% filter(bpoint %in% bpts_logg_valid$bpoint) %>%
    as.data.frame() %>%
    left_join(typology24) %>%
    count(category24) %>%
    pivot_wider(
      names_from = category24,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  for(type in colnames(logg_valid)) {
    if (type == "sumBpts") next
    all_bpts_df[4,type] <- paste(round(logg_valid[[type]]/logg_valid$sumBpts*100,2),"%",sep="")
  }
  all_bpts_df[4,"N"] <- logg_valid[["sumBpts"]]
  all_bpts_df[4,"Frontier"] <- "Logging concession (valid)"
  
  ### En cours ####
  bpts_logg_encours <- bpoints.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "en cours")
  
  logg_encours <- bpoints.v %>% filter(bpoint %in% bpts_logg_encours$bpoint) %>%
    as.data.frame() %>%
    left_join(typology24) %>%
    count(category24) %>%
    pivot_wider(
      names_from = category24,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  for(type in colnames(logg_encours)) {
    if (type == "sumBpts") next
    all_bpts_df[5,type] <- paste(round(logg_encours[[type]]/logg_encours$sumBpts*100,2),"%",sep="")
  }
  all_bpts_df[5,"N"] <- logg_encours[["sumBpts"]]
  all_bpts_df[5,"Frontier"] <- "Logging concession (in process)"
  
  all_bpts_df <- all_bpts_df %>% select(5,2,6,3,1,7,4)
  all_bpts_df[is.na(all_bpts_df)] <- 0
  as.data.frame(all_bpts_df) %>% 
    knitr::kable(format = "latex",booktabs = T)
  
# Diff-in-Disc Twangiza mine ---------------------------------------------------
 source(file=here("DifInDisc_function.R"))
  
## Identify double frontiers and select twangiza concession
bpoints.v <- vect(here(dataInt,"bpoints_cov.shp"))
PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
PAlines <- as.lines(PAs.v)
PAlines <- project(PAlines,bpoints.v)

s3aMineExpl <- vect(here(data,"s3aMineClean.shp"))
s3aMineExpl$Ind <- c(1:nrow(s3aMineExpl))
s3aMineExpl <- project(s3aMineExpl,PAlines)
s3aMineExpl <- buffer(s3aMineExpl,-1)
s3aMineExpl$Ind <- c(1:nrow(s3aMineExpl))
i2aDblFrontMin <- intersect(PAlines,s3aMineExpl)
i2aDblFrontMin$length.intersection <- perim(i2aDblFrontMin)

MineFront <- i2aDblFrontMin %>% filter(Ind %in% c(215,973)) %>% mutate(year.mine0 = 2012)
MineFront <- aggregate(MineFront,by="year.mine0")

# Diff-in-Disc (long runtime!)

load(here(dataInt,"panel.tmf.Rdata"))
DblMin <- DiffInDisc(depvar = "fcover",
                    treat.time = "year.mine0",
                    treat.ind = "treat",
                   data=panel.y.df,
                   doubleFronts = MineFront,
                   years = 2000:2024,
                   covariates = c("x","y"),
                   level = 0.1,
                   bandwidth = 20000)

  # save(DblMin,file = here(data,"DblMin20000_200m.Rdata"))
  # load(file = here(data,"DblMin20000_200m.Rdata"))

  fig8 <- DblMin %>%
    mutate(est = if_else(rel_time==-1,0,est),
           est_l =if_else(rel_time==-1,0,est_l), 
           est_r = if_else(rel_time==-1,0,est_r),
           se_bc = if_else(rel_time==-1,0,se_bc)) %>%
    ggplot(aes(x=year,y=est)) +
    geom_point(color = "#c7842a",size=.6) +
    geom_line(color = "#c7842a",linewidth = 0.6) +
    geom_line(aes(y=est_l),color = "#766297",linewidth = 0.6,linetype = "dashed") +
    geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.6,linetype = "dotdash") +
    annotate("text", x=2019.5, y=-0.07, label= "Forest cover inside",size=2,color = "#64ad6b",angle = -8) +
    annotate("text", x=2018, y=-0.1, label= "Forest cover outside",size=2,color = "#766297",angle = -14) +
    geom_vline(xintercept=2012, linetype="dotted", color = "darkgrey",linewidth = .4) +
    geom_vline(xintercept=2010, linetype="dotted", color = "darkgrey",linewidth = .4) +
    geom_vline(xintercept=2005, linetype="dotted", color = "darkgrey",linewidth = .4) +
    geom_vline(xintercept=2009, linetype="dotted", color = "darkgrey",linewidth = .4) +
    geom_vline(xintercept=2020, linetype="dotted", color = "darkgrey",linewidth = .4) +
    geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4,linewidth = .2) +
    geom_hline(yintercept=0, linetype="solid", color = "black",linewidth = .2) +
    annotate("text", x=2014.5, y=0.15, label= "Production phase",size=1.8) +
    annotate("text", x=2007, y=0.15, label= "Exploration phase",size=1.8) +
    annotate("text", x=2014, y=0.11, label= "Replacement agreement",size=1.8) +
    geom_segment(x = 2005.2, y = 0.14, xend = 2008.8, yend = 0.14,
                 arrow = arrow(length = unit(0.02, "npc"), ends = "both",type = "closed"),size=.4,linewidth = .4) +
    geom_segment(x = 2012.2, y = 0.14, xend = 2019.5, yend = 0.14,
                 arrow = arrow(length = unit(0.02, "npc"), ends = "both",type = "closed"),size=.4,linewidth = .4) +
    geom_curve(x = 2013, y = 0.113, xend = 2010.2, yend = 0.14,
               arrow = arrow(length = unit(0.02, "npc"), ends = "last",type = "closed"),size=.4,
               curvature = 0.3,linewidth = .4) +
    ylab("Difference in forest cover discontinuity") +
    xlab("Year") +
    theme_bw(base_size = 7) +
    theme(text = element_text(size = 7))
  
  ggsave(here(figures,"fig8.pdf"),plot=fig8,
         width = 90, height = 80, units = "mm", dpi = 300)
  
# Diff-in-Disc timber concessions -------------------------------------------------

# Identify Dbfs
bpoints.v <- vect(here(dataInt,"bpoints_cov.shp"))
PAs.v <- vect(here(data,"s2PAsClean_tmf.shp"))
  
PAlines <- as.lines(PAs.v)
PAlines <- project(PAlines,bpoints.v)
s3bTimber <- vect(here(data,"Forest_concession_agreements.shp"))
s3bTimber <- project(s3bTimber,PAlines)
s3bTimber <- buffer(s3bTimber,-1)
i2bDblFrontTimb <- intersect(PAlines,s3bTimber)
i2bDblFrontTimb$length.intersection <- perim(i2bDblFrontTimb)
t1 <- i2bDblFrontTimb %>% 
  filter(objectid %in% c(30,33,59,38),PA!=51) %>%
  mutate(logg_yr = as.numeric(strftime(as.Date(date_attr,format ="%d/%m/%Y"),"%Y"))) %>%
  mutate(logg_yr = if_else(objectid == 44,2020,logg_yr))

# Diff-in-Disc
load(here(dataInt,"panel.tmf.Rdata"))

res.timb.dbf <- data.frame()
for (timb.con in unique(t1$objectid)) {
  print(timb.con)
 dbfs <- t1 %>% filter(objectid == timb.con)
DblTimb <- DiffInDisc(depvar = "fcover",
                                treat.time = "logg_yr",
                                treat.ind = "treat",
                                data=panel.y.df,
                                doubleFronts = dbfs,
                                years = 2000:2024,
                                level = 0.1,
                                bandwidth = 20000)
DblTimb$objectid <- timb.con
res.timb.dbf <- bind_rows(res.timb.dbf,DblTimb)
}
res.timb.dbf <- res.timb.dbf %>% 
  left_join(as.data.frame(t1),by="objectid") %>%
  rename(year = year.x,year_PA = year.y)

# save(res.timb.dbf,file = here(dataOut,"res.timb.dbf.200m.Rdata"))
# load(here(dataOut,"res.timb.dbf.200m.Rdata"))

fig9a <- res.timb.dbf %>% filter(objectid == 38) %>%
  mutate(est = if_else(rel_time==-1,0,est),
         est_l =if_else(rel_time==-1,0,est_l), 
         est_r = if_else(rel_time==-1,0,est_r),
         se_bc = if_else(rel_time==-1,0,se_bc)) %>%
  ggplot(aes(x=year,y=est)) +
  geom_point(color = "#c7842a") +
  geom_line(color = "#c7842a",linewidth = 0.8) +
  geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
  geom_vline(xintercept=2011, linetype="dotted", color = "black") +
  annotate("text", x=2014.3, y=0.15, label= "Year of assignment",size=3) +
  geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
  ylab("Difference in forest cover discontinuity") +
  xlab("Years since production") +
  geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
  geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
  annotate("text", x=2021, y=-0.09, label= "Forest cover inside",size=3,color = "#64ad6b",angle = -8) +
  annotate("text", x=2020.5, y=-0.13, label= "Forest cover outside",size=3,color = "#766297",angle = -18) +
  ylim(c(-0.15,0.2)) +
  theme_bw(base_size = 7)

  fig9b <- res.timb.dbf %>% filter(objectid == 30) %>%
                                     mutate(est = if_else(rel_time==-1,0,est),
                                            est_l =if_else(rel_time==-1,0,est_l), 
                                            est_r = if_else(rel_time==-1,0,est_r),
                                            se_bc = if_else(rel_time==-1,0,se_bc)) %>%
                                     ggplot(aes(x=year,y=est)) +
                                     geom_point(color = "#c7842a") +
                                     geom_line(color = "#c7842a",linewidth = 0.8) +
                                     geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
                                     geom_vline(xintercept=2011, linetype="dotted", color = "black") +
                                     geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
                                     ylab("Difference in forest cover discontinuity") +
                                     xlab("Years since production") +
                                     geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
                                     geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
                                     ylim(c(-0.15,0.2)) +
                                     theme_bw(base_size = 7)

fig9c <- res.timb.dbf %>% filter(objectid == 33) %>%
                                  mutate(est = if_else(rel_time==-1,0,est),
                                         est_l =if_else(rel_time==-1,0,est_l), 
                                         est_r = if_else(rel_time==-1,0,est_r),
                                         se_bc = if_else(rel_time==-1,0,se_bc)) %>%
                                  ggplot(aes(x=year,y=est)) +
                                  geom_point(color = "#c7842a") +
                                  geom_line(color = "#c7842a",linewidth = 0.8) +
                                  geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
                                  geom_vline(xintercept=2011, linetype="dotted", color = "black") +
                                  geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
                                  ylab("Difference in forest cover discontinuity") +
                                  xlab("Years since production") +
                                  geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
                                  geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
                                  ylim(c(-0.15,0.2)) +
                                  theme_bw(base_size = 7)
                                 
fig9d <- res.timb.dbf %>% filter(objectid == 59) %>%
                                  mutate(est = if_else(rel_time==-1,0,est),
                                         est_l =if_else(rel_time==-1,0,est_l), 
                                         est_r = if_else(rel_time==-1,0,est_r),
                                         se_bc = if_else(rel_time==-1,0,se_bc)) %>%
                                  ggplot(aes(x=year,y=est)) +
                                  geom_point(color = "#c7842a") +
                                  geom_line(color = "#c7842a",linewidth = 0.8) +
                                  geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
                                  geom_vline(xintercept=2011, linetype="dotted", color = "black") +
                                  geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
                                  ylab("Difference in forest cover discontinuity") +
                                  xlab("Years since production") +
                                  geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
                                  geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
                                  annotate("text", x=2017.5, y=-0.05, label= "Forest cover inside",size=3,color = "#64ad6b",angle = -10) +
                                  annotate("text", x=2018, y=-0.13, label= "Forest cover outside",size=3,color = "#766297",angle = -20) +
                                  ylim(c(-0.15,0.2)) +
                                  theme_bw(base_size = 7)

fig9 <- fig9a + fig9b + fig9c + fig9c + plot_annotation(tag_levels ="a")

ggsave(here(figures,"fig9.pdf"),plot=fig9,
       width = 7*1.5, height = 6*1.5, units = "in", dpi = 300)

  