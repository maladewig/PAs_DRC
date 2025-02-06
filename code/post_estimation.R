
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
       install=T)

#library(rdmulti)

projectfolder<-here::here()
data <- here(projectfolder,"dataPrep")
results <- here(projectfolder,"results")
figures <- here(projectfolder,"figures")
tables <- here(projectfolder,"tables")

source(here("code","functions.R"))

load(file=here(data,"bpts.diff.cov.Rdata"))
bpts.diff.v <- vect(here(data,"bpts15000DiffTMF.shp"))

s2PAs <- vect(here(data,"s2PAsClean.shp"))
s2PAs <- project(s2PAs,crs(bpts.diff.v))

# Load grd results from Orion and bpts 
load(here(results,"r1aDef500_40000.Rdata"))
load(here(results,"r1bDegr500_40000.Rdata"))
load(here(results,"r1cCov500_40000.Rdata"))

PAs.df <- as.data.frame(vect((here(data,"s2PAsClean.shp"))))

# FIGURE 3 ---------------------------------------------------------------------
bpoints <- bpts.diff$bpoint

res.list <- list(r1aDef,r1bDegr,r1cCov)

AnnATTs <- function(res.df) {
  res.df <- res.df %>% filter(bpoint %in% bpoints)
  AvgAnnualATT(res.df = res.df,years = 2000:2022)
}
AnnATT.list <- lapply(res.list, AnnATTs)

# Merge DFs of different outcomes
r1AnnualATTs.df <- data.frame()
outcomes <- c("deforestation","degradation","forest cover")
for (df in 1:length(AnnATT.list)){
  res.df <- AnnATT.list[[df]] %>% mutate(outcome = outcomes[df])
  r1AnnualATTs.df <- rbind(r1AnnualATTs.df,res.df)
}

save(r1AnnualATTs.df,file=here(results,"r2AnnualATTs.df.Post2000.TMF.Rdata"))

p1 <- r1AnnualATTs.df %>% filter(outcome %in% c("degradation","deforestation")) %>%
  ggplot(aes(x=year, y=ATT, color = outcome,group = outcome,fill = outcome,shape = outcome)) + 
  geom_line(position=position_dodge(width=0.5),linewidth=.8) +
  geom_point(size=2.3,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="solid", color = "black") +
  geom_ribbon(aes(ymin=c.low, ymax=c.high),alpha=0.18,position=position_dodge(width=0.5)) +
  theme_bw() +
  scale_color_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
  scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
  scale_shape_manual(values = c(15,16,17,18,8)) +
  xlab("Year") +
  ylab("Discontinuity in annual disturbances") +
  guides(fill=guide_legend(title=" "),color=guide_legend(title=" "),shape=guide_legend(title=" "))

p2 <- r1AnnualATTs.df %>% filter(outcome %in% c("forest cover")) %>%
  ggplot(aes(x=year, y=ATT, color = outcome,group = outcome,fill = outcome)) + 
  geom_line(position=position_dodge(width=0.5),linewidth=.8) +
  geom_point(size=2.3,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="solid", color = "black") +
  geom_ribbon(aes(ymin=c.low, ymax=c.high),alpha=0.18,position=position_dodge(width=0.5)) +
  theme_bw() +
  scale_color_manual(values = c("#c7842a","#d8c656","#7d6756")) +
  scale_fill_manual(values = c("#c7842a","#d8c656","#7d6756")) +
  scale_shape_manual(values = c(15,16,17,18,8)) +
  xlab("Year") +
  ylab("Discontinuity in undisturbed forest") +
  guides(fill=guide_legend(title=" "),color=guide_legend(title=" "))

p2 + p1 + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')

# FIGURE 4: Sankey ---------------------------------------------------------------------

load(here(results,"r1aDef500_40000.Rdata"))
load(here(results,"r1bDegr500_40000.Rdata"))
load(here(results,"r1cCov500_40000.Rdata"))

Sankeyplot <- plotSankey(r1aDef,r1cCov)

Sankeyplot

# FIGURE 5: Country discontinuities ---------------------------------------------------------------------

bpoints <- bpts.diff$bpoint

  ## Forest cover ####

  est.df = r1cCov %>% filter(bpoint %in% bpoints)
  subsample_cat1 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn %in% c("Ia","Ib"))
  subsample_cat2 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "II")
  subsample_cat4 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "IV")
  subsample_cat6 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "VI")
  subsample_no_cat <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(is.na(cat_uicn))
  
  res.cat1 <- AvgAnnualATT(res.df = subsample_cat1) %>% mutate(cat = "I")
  res.cat2 <- AvgAnnualATT(res.df = subsample_cat2) %>% mutate(cat = "II")
  res.cat4 <- AvgAnnualATT(res.df = subsample_cat4) %>% mutate(cat = "IV")
  res.cat6 <- AvgAnnualATT(res.df = subsample_cat6) %>% mutate(cat = "VI")

  res.by.iucncat <-res.cat1 %>% bind_rows(res.cat2,res.cat4,res.cat6) %>% mutate(outcome = "fcover")
  
  ## Deforestation ####
  est.df = r1aDef %>% filter(bpoint %in% bpoints)
  subsample_cat1 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn %in% c("Ia","Ib"))
  subsample_cat2 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "II")
  subsample_cat4 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "IV")
  subsample_cat6 <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(cat_uicn == "VI")
  subsample_no_cat <- PAs.df %>% select(PA,cat_uicn) %>% right_join(est.df) %>% filter(is.na(cat_uicn))
  
  res.cat1 <- AvgAnnualATT(res.df = subsample_cat1) %>% mutate(cat = "I")
  res.cat2 <- AvgAnnualATT(res.df = subsample_cat2) %>% mutate(cat = "II")
  res.cat4 <- AvgAnnualATT(res.df = subsample_cat4) %>% mutate(cat = "IV")
  res.cat6 <- AvgAnnualATT(res.df = subsample_cat6) %>% mutate(cat = "VI")

  res.by.iucncat <- res.cat1 %>% bind_rows(res.cat2,res.cat4,res.cat6) %>% mutate(outcome = "def") %>% bind_rows(res.by.iucncat)
  save(res.by.iucncat,file = here(results,"r3by.iucncat.Rdata"))
p1 <- res.by.iucncat %>% 
  filter(year %in% c(2002,2012,2022),outcome == "fcover") %>%
  mutate(year = factor(year,levels = c(2022,2012,2002)),cat = factor(cat,levels = c("I","II","IV","VI"))) %>%
  ggplot(aes(y=cat,x=ATT,color = year, group = year,shape=year)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin=c.low, xmax=c.high),linewidth=.8,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black") +
  theme_bw() +
  scale_color_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  scale_fill_manual(values = c("#766297","#c7842a","#d8c656","#7d6756","#440154")) +
  scale_shape_manual(values = c(15,16,17,18,8,3)) +
  xlab("Discontinuity in forest cover") +
  ylab("IUCN category")

p2 <- res.by.iucncat %>% filter(year %in% c(2015:2022),outcome == "def") %>%
  mutate(cat = factor(cat,levels = c("I","II","IV","VI"))) %>%
  ggplot(aes(y=cat,x=ATT,color = year, group = year)) +
  geom_point(position = position_dodge(width = 0.7)) +
  geom_errorbarh(aes(xmin=c.low, xmax=c.high),linewidth=.8,height = .51,position = position_dodge(width = 0.7)) +
  geom_vline(xintercept=0, linetype="solid", color = "black") +
  theme_bw() +
  viridis::scale_color_viridis(discrete=F) +
  xlim(c(-0.007,0.02)) +
  xlab("Discontinuity in deforestation") +
  ylab("IUCN category")

p1 + p2 + plot_annotation(tag_levels = "a")

# FIGURE 6: BOXPLOTS BY IUCN AND PROTECTION MECHANISM ---------------------------------------------------------------------

load(file=here(data,"bpts.diff.cov.Rdata"))

resDF <- r1aDef %>% 
  rename(left_def = est_l, right_def = est_r, disc_def = est_rob) %>%
  filter(!is.na(bpoint),on.border!=1,outside.tmf!=1)

resDF <- r1cCov %>% 
  select(bpoint, year,left_cov = est_l, right_cov = est_r, disc_cov = est_rob) %>% 
  right_join(resDF, by = c("bpoint","year")) 
years = c(2002,2010,2022)
y <- min(years)
df <- data.frame()
for (year in years) {
  first_y = year -5
  bpts.diff <- resDF %>% arrange(bpoint,year) %>% 
      filter(!is.na(right_cov*right_def*left_cov*left_def),
             year %in% c(first_y:y)) %>%
      group_by(bpoint) %>%
      summarise(left_cov = min(left_cov),
                right_cov = min(right_cov),
                left_def = mean(left_def),
                right_def = mean(right_def),
      ) %>%
      ungroup() %>% 
      mutate(mechanism = case_when( # Metrics from Jenks for year 2000 (average left + right)
        right_def <=.005 & left_def >=.005 & right_cov > 0.335 ~ "Contained",
        right_def >= 0.005 ~ "Sprawling",
        right_cov <=0.335 ~ "Exhausted",
        left_cov > .79 & right_cov > .79 & left_def <= 0.005 & right_def <= 0.005 ~ "Dormant",
        right_cov > 0.335 & right_def < 0.005 ~ "Consolidated"
      )) %>%
      select(bpoint,mechanism) %>%
      right_join(bpts.diff)
    colnames(bpts.diff)[which(colnames(bpts.diff) == "mechanism")] <- glue::glue("mechanism{year}")
}

save(bpts.diff,file=here(data,"bpts.diff.cov.Rdata"))

p1 <-bpts.diff %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(mechanism2022),iucn != "Ib") %>%
  mutate(mechanism2022 = fct_reorder(mechanism2022, loggroads)) %>%
  mutate(mechanism2022 = factor(mechanism2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=loggroads, x=mechanism2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5,lwd = 0.3) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Dist. to forest road (km)") +
  theme(axis.text.x = element_text(angle = -25, vjust = -.7,size = 8),
        axis.title.y = element_text(size = 8))

p2 <- bpts.diff %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(mechanism2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(mechanism2022 = fct_reorder(mechanism2022, conflictsAft2010)) %>%
  mutate(mechanism2022 = factor(mechanism2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=conflictsAft2010, x=mechanism2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5,lwd = 0.3) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Conflict incidences (#)") +
  ylim(0,180) +
  theme(axis.text.x = element_text(angle = -25, vjust = -.7,size=8),
        axis.title.y = element_text(size = 8))


p1 / p2 + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# FIGURE 6: HETEROGENEITIES BY PA ---------------------------------------------------------------------

rob.1 <- r1cCov %>% filter(PA_year > 2000)
res.PAs <- data.frame()
for (i in unique(r1cCov$PA)) {
  df1 <- r1cCov %>% filter(PA == i)
  df2 <- AvgAnnualATT(res.df = df1,years = 2000:2022,blockbstrap=F)
  df2$PA <- i
  res.PAs <- bind_rows(res.PAs,df2)
}

res.PAs <- PAs.df %>% 
  rename(PA_year = year) %>% 
  right_join(res1) %>%
  mutate(PA_year = if_else(PA ==22,2011,PA_year)) %>%
  mutate(PA_year = if_else(PA ==23,2006,PA_year)) %>%
  mutate(rel_time = year -PA_year) %>%
  mutate(border.PA = if_else(PA %in% c(12,22,31,34,53),1,0)) %>%
  mutate(nom_orig = if_else(PA == 43,"Rubi-Tele",nom_orig),nom_orig = if_else(PA == 14,"Salonga I",nom_orig),nom_orig = if_else(PA == 32,"Salonga II",nom_orig))

save(res.PAs,file=here(results,"r4by.PA.Rdata"))
load(file=here(results,"r4by.PA.Rdata"))

res.PAs %>% 
  filter(year %in% c(2012,2022),!(nom_orig %in% c("Rutshuru","Mondo-Missa","Luki","Luama-Katanga")),border.PA !=1) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(y=nom_orig,x=ATT,group = year,color = year)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin=c.low, xmax=c.high),linewidth=.7,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black") +
  scale_color_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756","#440154")) +
  xlab("Discontinuity in undisturbed forest cover") +
  ylab("Name of PA") +
  theme_bw()

# TABLE 1: MINING / LOGGING x CONSERVATION ------------------------------------
  
  ## Mining ####
  
  ### Expl. permit ####
  mining_dbf <- vect(here(data,"minconcessions_dbf.shp"))
  mining_dbf <- mining_dbf %>%
    project(bpts.diff.v) %>% buffer(width=50) %>% rename(PAnr = PA)
  
  bpts_mining_expl <- bpts.diff.v %>%
    intersect(mining_dbf) %>%
    filter(type__2 %in% c(6210301,6210100,6210101,6210201))
  
  mine_expl <- bpts.diff %>% filter(bpoint %in% bpts_mining_expl$bpoint) %>%
    filter(year==2022) %>% 
    count(mechanism2022) %>%
    pivot_wider(
      names_from = mechanism2022,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  all_bpts_df <- data.frame()
  for(type in colnames(mine_expl)) {
    if (!(type %in% c("sumBpts","Frontier"))) {
    all_bpts_df[1,type] <- paste(round(mine_expl[[type]]/mine_expl$sumBpts*100,2),"%",sep="")
    }
  }
  all_bpts_df[1,"N"] <- mine_expl[["sumBpts"]]
  all_bpts_df[1,"Frontier"] <- "Mining concession (exploitation permit)"
  
  ### Research permit ####
  bpts_mining_rech <- bpts.diff.v %>%
    intersect(mining_dbf) %>%
    filter(type__2 %in% c(6210200))
  
  mine_rech <- bpts.diff %>% filter(bpoint %in% bpts_mining_rech$bpoint) %>%
    filter(year==2022) %>% 
    count(mechanism2022) %>%
    pivot_wider(
      names_from = mechanism2022,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(sumBpts = rowSums(.))
  
  for(type in colnames(mine_rech)) {
    if (type == "sumBpts") next
    all_bpts_df[2,type] <- paste(round(mine_rech[[type]]/mine_rech$sumBpts*100,2),"%",sep="")
  }
  all_bpts_df[2,"N"] <- mine_rech[["sumBpts"]]
  all_bpts_df[2,"Frontier"] <- "Mining concession (research permit)"
  
  ## Logging ####
  
  logg_dbf <- vect(here(data,"loggconcessions_dbf.shp"))
  logg_conc <- vect(here(data,"Forest_concession_agreements.shp"))
  logg_conc <- project(logg_conc,bpts.diff.v)
  
  ### Active ####
  bpts_logg_amenage <- bpts.diff.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "amenage")
  
  logg_amenage <- bpts.diff %>% filter(bpoint %in% bpts_logg_amenage$bpoint) %>%
    filter(year==2022) %>% 
    count(mechanism2022) %>%
    pivot_wider(
      names_from = mechanism2022,
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
  bpts_logg_valid <- bpts.diff.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "valide")
  
  logg_valid <- bpts.diff %>% filter(bpoint %in% bpts_logg_valid$bpoint) %>%
    filter(year==2022) %>% 
    count(mechanism2022) %>%
    pivot_wider(
      names_from = mechanism2022,
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
  bpts_logg_encours <- bpts.diff.v %>%
    intersect(logg_conc) %>%
    filter(statu_amgt == "en cours")
  
  logg_encours <- bpts.diff %>% filter(bpoint %in% bpts_logg_encours$bpoint) %>%
    filter(year==2022) %>% 
    count(mechanism2022) %>%
    pivot_wider(
      names_from = mechanism2022,
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
  
  as.data.frame(all_bpts_df) %>% 
    knitr::kable(format = "latex",booktabs = T)
  
# FIGURE 8: TWANGIZA MINE DISCONTINUITIES --------------------------------------

## Identify double frontiers and select twangiza concession

PAlines <- as.lines(s2PAs)
PAlines <- project(PAlines,bpts.diff.v)

s3aMineExpl <- vect(here(data,"s3aMineClean.shp"))
s3aMineExpl$Ind <- c(1:nrow(s3aMineExpl))
s3aMineExpl <- project(s3aMineExpl,PAlines)
s3aMineExpl <- buffer(s3aMineExpl,-1)
s3aMineExpl$Ind <- c(1:nrow(s3aMineExpl))
i2aDblFrontMin <- intersect(PAlines,s3aMineExpl)
i2aDblFrontMin$length.intersection <- perim(i2aDblFrontMin)


t1 <- i2aDblFrontMin %>% filter(Ind %in% c(215,973)) %>% mutate(year.mine0 = 2012)
t1 <- aggregate(t1,by="year.mine0")

# Diff-in-Disc

load(here(data,"i1cCovPanel500.Rdata"))

DblMin <- DiffInDisc(depvar = "fcover",
                    treatvar = "year.mine0",
                   data=i1cCovPanel,
                   doubleFronts = t1,
                   df.crs = crs(s1cTMFCov),
                   #covariates = c("x","y"),
                   level = 0.1,
                   bandwidth = 20000)

  
  save(DblMin,file = here(results,"r5DblMin20000.Rdata"))
  load(file = here(results,"r5DblMin20000.Rdata"))

  DblMin %>%
    mutate(est = if_else(rel_time==-1,0,est),
           est_l =if_else(rel_time==-1,0,est_l), 
           est_r = if_else(rel_time==-1,0,est_r),
           se_bc = if_else(rel_time==-1,0,se_bc)) %>%
    ggplot(aes(x=year,y=est)) +
    geom_point(color = "#c7842a") +
    geom_line(color = "#c7842a",linewidth = 0.8) +
    geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
    geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
    geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
    annotate("text", x=2017.5, y=-0.05, label= "Forest cover inside",size=3,color = "#64ad6b",angle = -10) +
    annotate("text", x=2018, y=-0.13, label= "Forest cover outside",size=3,color = "#766297",angle = -20) +
    geom_vline(xintercept=2012, linetype="dotted", color = "black") +
    geom_vline(xintercept=2010, linetype="dotted", color = "black") +
    geom_vline(xintercept=2005, linetype="dotted", color = "black") +
    geom_vline(xintercept=2009, linetype="dotted", color = "black") +
    geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
    annotate("text", x=2014.5, y=0.2, label= "Production phase",size=3) +
    annotate("text", x=2007, y=0.2, label= "Exploration phase",size=3) +
    annotate("text", x=2014, y=0.11, label= "Replacement agreement",size=3) +
    geom_segment(x = 2005.2, y = 0.17, xend = 2008.8, yend = 0.17,
                 arrow = arrow(length = unit(0.03, "npc"), ends = "both")) +
    geom_segment(x = 2012.2, y = 0.17, xend = 2019.5, yend = 0.17,
                 arrow = arrow(length = unit(0.03, "npc"), ends = "both")) +
    geom_curve(x = 2013, y = 0.113, xend = 2010.2, yend = 0.17,
               arrow = arrow(length = unit(0.03, "npc"), ends = "last"),
               curvature = 0.3) +
    ylab("Difference in forest cover discontinuity") +
    xlab("Year") +
    theme_bw()

# FIGURE 9: Timber concessions -------------------------------------------------

# Identify Dbfs

PAlines <- as.lines(s2PAs)
PAlines <- project(PAlines,bpts.diff.v)
s3bTimber <- vect(here(data,"Forest_concession_agreements.shp"))
s3bTimber <- project(s3bTimber,PAlines)
s3bTimber <- buffer(s3bTimber,-1)
#s3bTimber$Ind <- c(1:nrow(s3bTimber))
i2bDblFrontTimb <- intersect(PAlines,s3bTimber)
i2bDblFrontTimb$length.intersection <- perim(i2bDblFrontTimb)
t1 <- i2bDblFrontTimb %>% 
  #filter(objectid %in% c(37,28,44,30,33,59,23,38),PA!=51) %>%
  filter(objectid %in% c(30,33,59,38),PA!=51) %>%
  mutate(logg_yr = as.numeric(strftime(as.Date(date_attr,format ="%d/%m/%Y"),"%Y"))) %>%
  mutate(logg_yr = if_else(objectid == 44,2020,logg_yr))

#i2bDblFrontTimb <- i2bDblFrontTimb[i2bDblFrontTimb$length.intersection > 5000]
#writeVector(i2bDblFrontTimb,here(dataInt,"i2aDblFrontTimb_5000.shp"))

# Diff-in-Disc

load(here(data,"i1cCovPanel500.Rdata"))

res.timb.dbf <- data.frame()
for (timb.con in unique(t1$objectid)) {
  print(timb.con)
 dbfs <- t1 %>% filter(objectid == timb.con)
DblTimb <- DiffInDisc(depvar = "fcover",
                      data=i1cCovPanel,
                      treatvar = "logg_yr",
                      doubleFronts = dbfs,
                      df.crs = crs(s1cTMFCov),
                      covariates = c("x","y"),
                      level = 0.1,
                      bandwidth = 20000)
DblTimb$objectid <- timb.con
res.timb.dbf <- bind_rows(res.timb.dbf,DblTimb)
}
res.timb.dbf <- res.timb.dbf %>% 
  left_join(as.data.frame(t1),by="objectid") %>%
  rename(year = year.x,year_PA = year.y)

save(res.timb.dbf,file = here(results,"r5DblTimb20000.Rdata"))

plots <- list()
p<-1
#for (id in c(37,28,44,23)) {
for (id in c(38,30,33,59)) {
plots[[p]] <- res.timb.dbf %>% filter(objectid == id) %>%
  mutate(est = if_else(rel_time==-1,0,est),
         est_l =if_else(rel_time==-1,0,est_l), 
         est_r = if_else(rel_time==-1,0,est_r),
         se_bc = if_else(rel_time==-1,0,se_bc)) %>%
  ggplot(aes(x=year,y=est)) +
  geom_point(color = "#c7842a") +
  geom_line(color = "#c7842a",linewidth = 0.8) +
  geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
  geom_vline(xintercept=2011, linetype="dotted", color = "black") +
  annotate("text", x=2006, y=0.055, label= "Year of assignment",size=3) +
  geom_curve(x = 2009, y = 0.05, xend = 2010.8, yend = 0.035,
             arrow = arrow(length = unit(0.03, "npc"), ends = "last"),
             curvature = 0.3) +
  geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
  ylab("Difference in forest cover discontinuity") +
  xlab("Years since production") +
  geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
  geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
  annotate("text", x=2017.5, y=-0.05, label= "Forest cover inside",size=3,color = "#64ad6b",angle = -10) +
  annotate("text", x=2018, y=-0.13, label= "Forest cover outside",size=3,color = "#766297",angle = -20) +
  ylim(c(-0.15,0.2)) +
  theme_bw()

p <- p+1
}
library(patchwork)
plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plot_annotations(level ="a")


# APPENDIX ---------------------------------------------------------------------

## FIGURE A1: RDD plots with covariates ####
distCells <- vect(here(dataInt,"distCells.shp"))
dist.df <- data.frame(values(distCells),geom(distCells))
x <- dist.df[["dist2cutof"]]
plots <- list()
p <- 1
covariates <- c("s3bAgrSuit" ,"s3cSlope" ,"s3dAltitude","annual_precipitation","annual_temp","growing_period_days")
for (cov in covariates) {
  print(cov)
  cov_rast <- rast(here(dataPrep,glue::glue("{cov}.tif")))
  if (crs(distCells) != crs(cov_rast)) {
    distCells <- project(distCells,cov_rast)
  }
y <- extract(cov_rast,distCells,ID=F)
plots[[p]] <- rdrobust::rdplot(y=y[[cov]],x=x)
p <- p+1
}
RDplots <- plots
saveRDS(plots, file=here(results,"RDplots.RData"))

p1 <- plots[[1]]$rdplot +
  ggtitle("Agricultural suitability") +
  labs(x="Distance to PA",y="Agricultural suitability score")

p2 <- plots[[2]]$rdplot +
  ggtitle("Slope") +
  labs(x="Distance to PA",y="Slope")

p3 <- plots[[3]]$rdplot +
  ggtitle("Altitude") +
  labs(x="Distance to PA",y="Altitude (m)")

p4 <- plots[[4]]$rdplot +
  ggtitle("Precipitation") +
  labs(x="Distance to PA",y="Precipitation (mm)")

p5 <- plots[[5]]$rdplot +
  ggtitle("Temperature") +
  labs(x="Distance to PA",y="Temperature (C)")

p6 <- plots[[6]]$rdplot +
  ggtitle("Growing period length") +
  labs(x="Distance to PA",y="Growing degree days")

ggsave(filename = here(figures,"RDplot_GDD.pdf"),device = "pdf")

p1 + geom_line(aes(y,x,color="yellow"))

library(patchwork)
(p1 + p2)/(p3 + p4)/(p5 + p6) & plot_annotation(tag_levels = "a")


## FIGURE A2 RDs pre-establishment ####

rob.1 <- r1cCov %>% filter(PA_year > 2000)
res1 <- data.frame()
 for (i in unique(r1cCov$PA)) {
   t1 <- r1cCov %>% filter(PA == i)
   t2 <- AvgAnnualATT(res.df = t1,years = 2000:2022,blockbstrap=F)
   t2$PA <- i
   res1 <- bind_rows(res1,t2)
 }

robust.df <- PAs.df %>% 
  rename(PA_year = year) %>% 
  right_join(res1) %>%
  mutate(PA_year = if_else(PA ==22,2011,PA_year)) %>%
  mutate(PA_year = if_else(PA ==23,2006,PA_year)) %>%
  mutate(rel_time = year -PA_year) %>%
  mutate(border.PA = if_else(PA %in% c(12,22,31,34,53),1,0)) %>%
  mutate(nom_orig = if_else(PA == 43,"Rubi-Tele",nom_orig),nom_orig = if_else(PA == 14,"Salonga I",nom_orig),nom_orig = if_else(PA == 32,"Salonga II",nom_orig))

save(robust.df,file=here(results,"PA.est.Rdata"))
load(file=here(results,"PA.est.Rdata"))

p1 <-robust.df %>% 
  filter(PA_year > 2000,rel_time == -1 | year == 2022,border.PA == 0,!(nom_orig %in% c("Yangambi"))) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(y=nom_orig,x=ATT,group = year,color = year)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin=c.low, xmax=c.high),linewidth=.8,height = .51,position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype="solid", color = "black") +
  theme_bw()

ptest <- r1cCov %>% group_by(PA) %>% filter(year == 2022) %>% summarise(fcover.disc = mean(est_rob))
ptest <- r1aDef %>% group_by(bpoint) %>% 
  filter(year %in% 2018:2022) %>% 
  summarise(def.disc = mean(est_rob),PA = mean(PA)) %>% 
  group_by(PA) %>%
  summarise(def.disc = mean(def.disc)) %>%
  ungroup() %>%
  select(def.disc) %>%
  bind_cols(ptest)

ggplot(ptest,aes(x=fcover.disc,y=def.disc)) +
  geom_point() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

## TABLE A1: PA heterogeneity --------------------------------------------------

PA.cats <- bpts.diff %>%  
  count(PA,mechanism2022) %>%
  filter(!is.na(mechanism2022)) %>%
  pivot_wider(
    names_from = mechanism2022,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(N = Consolidated + Dormant+ Exhausted + Sprawling+ Contained) %>%
  mutate(Dormant = round(Dormant / N,2),
         Exhausted = round(Exhausted / N,2),
         Sprawling = round(Sprawling / N,2),
         Contained = round(Contained / N,2),
         Consolidated = round(Consolidated / N,2)) %>%
  filter(N>10)

s2PAs <- vect(here(dataPrep,"s2PAsClean.shp"))
s2PAs$area <- terra::expanse(s2PAs,unit="ha")

tabA1 <- as.data.frame(s2PAs) %>% 
  select(PA,name = nom_ap,IUCN = cat_uicn,area) %>%
  right_join(PA.cats) %>%
  select(-PA) %>%
  arrange(name) %>%
  knitr::kable(format = "latex",booktabs = T)

## FIGURE A3: Boxplots by IUCN -------------------------------------------------
PAs <-PAs.df %>% select(PA,iucn = cat_uicn,type_gouv)


bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(ptype2022 = fct_reorder(ptype2022, s3bAgrS)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=s3bAgrS, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Agricultural suitability")

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(ptype2022 = fct_reorder(ptype2022, s3Accss)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=s3Accss, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Accessibility") 

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(ptype2022 = fct_reorder(ptype2022, GRIP)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=GRIP, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Distance to road (km)") 

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(ptype2022 = fct_reorder(ptype2022, loggroads_old)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=loggroads_old, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("loggroads_old") 

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(ptype2022 = fct_reorder(ptype2022, loggroads)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=loggroads, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Distance to forest road (km)") 

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(ptype2022 = fct_reorder(ptype2022, s3dAltt)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=s3dAltt, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Altitude") 

bpts.diff %>% st_drop_geometry(bpts.diff) %>%
  mutate(iucn = if_else(is.na(iucn),"--",iucn)) %>%
  filter(!is.na(ptype2022),iucn != "Ib") %>%
  mutate(conflictsAft2010 = if_else(is.na(conflictsAft2010),0,conflictsAft2010))  %>%
  mutate(ptype2022 = fct_reorder(ptype2022, conflictsAft2010)) %>%
  mutate(ptype2022 = factor(ptype2022, levels=c("Dormant", "Contained", "Consolidated", "Sprawling","Exhausted")),
         iucn = factor(iucn, levels=c("Ia", "II", "IV", "VI","--"))) %>%
  ggplot(aes(y=conflictsAft2010, x=ptype2022,fill = iucn)) + 
  geom_boxplot(position="dodge", alpha=0.5) +
  viridis::scale_fill_viridis(discrete=T, name="IUCN") +
  theme_bw()  +
  xlab("") +
  ylab("Conflict events since 2010") +
  ylim(0,180)

## TABLE A2 ####

ptype_IUCN <- bpts.diff %>% count(iucn,ptype2022) 
ptype_IUCN <- ptype_IUCN %>% 
  group_by(iucn) %>% 
  summarise(tot = sum(n)) %>% 
  ungroup() %>% 
  right_join(ptype_IUCN) %>% 
  mutate(share = n/tot) %>%
  pivot_wider(id_cols = c(iucn,tot),names_from = ptype2022, values_from = share) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_iucn.tex"))

ptype_gouv <- bpts.diff %>% 
  count(type_gouv,ptype2022)

ptype_gouv <- ptype_gouv %>% 
  group_by(type_gouv) %>% 
  summarise(tot = sum(n)) %>% 
  ungroup() %>% 
  right_join(ptype_gouv) %>% 
  mutate(share = n/tot) %>%
  pivot_wider(id_cols = c(type_gouv,tot),names_from = ptype2022, values_from = share) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_gouv.tex"))

ptype_newroads <- bpts.diff %>% 
  count(new.loggroad_5km,ptype2022)

ptype_newroads <- ptype_newroads %>% 
  group_by(new.loggroad_5km) %>% 
  summarise(tot = sum(n)) %>% 
  ungroup() %>% 
  right_join(ptype_newroads) %>% 
  mutate(share = n/tot) %>%
  pivot_wider(id_cols = c(new.loggroad_5km,tot),names_from = ptype2022, values_from = share) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_newroads.tex"))

ptype_loggroads <- bpts.diff %>% 
  count(loggroads,ptype2022)

bpts.diff %>% 
  group_by(ptype2022) %>% 
  summarise(loggroads = mean(loggroads)) %>% 
  ungroup() %>% 
  mutate(var = "Logging road distance") %>%
  pivot_wider(id_cols = c(var),names_from = ptype2022, values_from = loggroads) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_loggroads.tex"))


ptype_access <- bpts.diff %>%
  group_by(ptype2022) %>% 
  summarise(access = mean(s3Accss)/60) %>%
  ungroup() %>%
  mutate(var = "Accessibility") %>%
  pivot_wider(id_cols = c(var),names_from = ptype2022, values_from = access) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_access.tex"))

ptype_access <- bpts.diff %>%
  group_by(ptype2022) %>% 
  summarise(altitude = mean(s3dAltt)) %>%
  ungroup() %>%
  mutate(var = "Altitude") %>%
  pivot_wider(id_cols = c(var),names_from = ptype2022, values_from = altitude) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_altitude.tex"))

ptype_road <- bpts.diff %>%
  group_by(ptype2022) %>% 
  summarise(road.dist = mean(road.dist.)/1000) %>%
  ungroup() %>%
  mutate(var = "Road distance (km)") %>%
  pivot_wider(id_cols = c(var),names_from = ptype2022, values_from = road.dist) %>%
  kbl(booktabs = TRUE,format = "latex",digits = 2) %>%
  save_kable(here(tables,"S1_tab_road.tex"))

## FIGURE A4 OLD CONCESSION BOUNDARIES ####

PAlines <- as.lines(vect(here(dataInt,"ITBW.old.shp")))
PAlines <- project(PAlines,bpts.diff.v)

s3aMineExpl <- vect(here(dataPrep,"s3aMineClean.shp"))
s3aMineExpl <- project(s3aMineExpl,PAlines)
s3aMineExpl <- buffer(s3aMineExpl,-1)
s3aMineExpl$Ind <- c(1:nrow(s3aMineExpl))
i2aDblFrontMin <- intersect(PAlines,s3aMineExpl)
i2aDblFrontMin$length.intersection <- perim(i2aDblFrontMin)
i2aDblFrontMin <- i2aDblFrontMin[i2aDblFrontMin$length.intersection > 5000]

t1 <- i2aDblFrontMin %>% filter(Ind %in% c(215)) %>% mutate(year.mine0 = 2012)
t1 <- aggregate(t1,by="year.mine0")

#writeVector(i2aDblFrontMin,here(dataInt,"i2aDblFrontMin_5000.shp"))

# Diff-in-Disc

load(here(dataInt,"i1cCovPanel500.Rdata"))

DblMin <- DiffInDisc(depvar = "fcover",
                     treatvar = "year.mine0",
                     data=i1cCovPanel,
                     doubleFronts = t1,
                     df.crs = crs(s1cTMFCov),
                     #covariates = c("x","y"),
                     level = 0.1,
                     bandwidth = 20000)


save(DblMin,file = here(results,"DblMinRob.Rdata"))
load(file = here(results,"DblMinRob.Rdata"))

DblMin %>%
  mutate(est = if_else(rel_time==-1,0,est),
         est_l =if_else(rel_time==-1,0,est_l), 
         est_r = if_else(rel_time==-1,0,est_r),
         se_bc = if_else(rel_time==-1,0,se_bc)) %>%
  ggplot(aes(x=year,y=est)) +
  geom_point(color = "#c7842a") +
  geom_line(color = "#c7842a",linewidth = 0.8) +
  geom_hline(yintercept=0, linetype="solid", color = "darkgrey") +
  geom_vline(xintercept=2010, linetype="dotted", color = "black") +
  geom_ribbon(aes(ymin = est - se_bc*1.96,ymax = est + se_bc*1.96),fill = "#c7842a",color = "#c7842a",alpha = 0.4) +
  geom_line(aes(y=est_l),color = "#766297",linewidth = 0.8) +
  geom_line(aes(y=est_r),color = "#64ad6b",linewidth = 0.8) +
  #annotate("text", x=2017.5, y=-0.05, label= "Forest cover inside",size=3,color = "#64ad6b",angle = -10) +
  #annotate("text", x=2018, y=-0.13, label= "Forest cover outside",size=3,color = "#766297",angle = -20) +
  ylab("Difference in forest cover discontinuity") +
  xlab("Year") +
  theme_bw()
  