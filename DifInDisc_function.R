#################GRDD Functions and plotting of results#########################

# Plot 

# - translate distance() function from sf to raster

# Calculate Distance from and to PA-#######################



######################### Diff in Disc #########################################

DiffInDisc <- function (depvar= "fcover",data, doubleFronts,treat.ind,treat.time, calculate.windows = T, bandwidth = 20000, level = 0.15, years = 2000:2024, covariates = NULL, locrand = F, ...) 
{
  
  #Define extent of analysis
  doubleFronts$id <- 1:nrow(doubleFronts)
  buff.cutoff <- buffer(doubleFronts,bandwidth)
  cutoff.ext <- ext(buff.cutoff)
  xmin <- cutoff.ext[1]
  xmax <- cutoff.ext[2]
  ymin <- cutoff.ext[3]
  ymax <- cutoff.ext[4]
  data <- data %>% filter(x > xmin, x < xmax,y > ymin,y < ymax, year %in% years)
  colnames(data)[which(colnames(data)==depvar)] <- "outcome"
  colnames(data)[which(colnames(data)==treat.ind)] <- "treat"
  names(doubleFronts)[which(names(doubleFronts)==treat.time)] <- "tg"
  groups = unique(doubleFronts$tg)
  columnames <- c("group","year","est","p_value","se_bc","se_rob","est_l","est_r","window_l","window_r","segment_length","rel_time","n_segments","n_l","n_r")
  results <- data.frame(matrix(ncol = length(columnames), nrow = length(groups)*length(years)))
  colnames(results) <- columnames
  row.ind <- 0
  
  cutoff.sf <- st_as_sf(doubleFronts)
  data.sf <- st_as_sf(data,coords = c("x","y"),crs=st_crs(cutoff.sf))
  
  data$dist2cutoff <- sf::st_distance(data.sf,cutoff.sf) 
  data <- data %>%
    mutate(dist2cutoff = if_else(treat == 1,dist2cutoff,dist2cutoff*-1))
  
  for(g in groups) {
    print(g)
    border.segments <- doubleFronts %>% tidyterra::filter(tg == g)
    border.segments <- aggregate(border.segments)
    buffed.segments <- buffer(border.segments,bandwidth)
    
    base_out <- data$outcome[data$year==(g-1)]

    for (yr in years) {
      row.ind <- row.ind + 1
      data.yr <- data %>% filter(year == yr) %>% mutate(diff_out = outcome - base_out)
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
        rd.est <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, x = data.yr$dist2cutoff, c = 0, ...), error = function(e) {
          err <<- TRUE
        })
        
        
        #rd.est.add <- tryCatch(rdrobust::rdrobust(y=data.yr$diff_out, 
        #                                         x = data.yr$dist2cutoff, c = 0, covs = as.matrix(data.yr[,covariates]), ...) ...), error = function(e) {
        #                                          err <<- TRUE
        #                                       })
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
        results[row.ind, "est_l"] <- rd.est$tau_bc[1]
        results[row.ind, "est_r"] <- rd.est$tau_bc[2]
        results[row.ind, "window_l"] <- round(rd.est[["bws"]][["h", 
                                                               1]]/1000, 1)
        results[row.ind, "window_r"] <- round(rd.est[["bws"]][["h", 
                                                               2]]/1000, 1)
        results[row.ind, "n_l"] <- rd.est$N_h[1]
        results[row.ind, "n_r"] <- rd.est$N_h[2]
      }
      
      print(yr)
    }
  }
  return(results)  
}