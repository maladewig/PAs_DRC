### NOTE! This scritpt has been run on an external server due to its computational demands ###

######################## Estimating Non-parametric GRDs #######################

# 0 Load packages and directories ----
library(tidyverse)
library(stringr)
library(here)
library(Matrix)

load(file = here(dataInt,"panel.tmf.Rdata"))
load(file= here(dataPrep,"bpoints.tmf.Rdata"))


projectfolder<-here::here()
dataPrep <- here(projectfolder,"dataPrep")
dataInt <- here(projectfolder,"dataInt")
dataOut <- here(projectfolder,"dataOut")
figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
tables <- here(projectfolder,"tables")

######################### Functions #############################
rdbw2d <- function(Y, X, t, b, p = 1, deriv = c(0,0), tangvec = NULL,
                   kernel = c("tri","triangular","epa","epanechnikov","uni","uniform","gau","gaussian"),
                   kernel_type = c("prod","rad"),
                   bwselect = c("mserd", "imserd", "msetwo", "imsetwo"),
                   method = c("dpi", "rot"), vce = c("hc1","hc0","hc2","hc3"),
                   bwcheck = 50 + p + 1, masspoints = c("check","adjust","off"),
                   C = NULL, scaleregul = 1, scalebiascrct = 1,
                   stdvars = TRUE){
  
  # Input error handling
  
  bwselect <- match.arg(bwselect)
  kernel <- match.arg(kernel)
  kernel_type <- match.arg(kernel_type)
  method <- match.arg(method)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  verbose <- FALSE
  
  d <- t # renaming the variable
  
  # Check Errors
  
  exit=0
  
  if (length(Y) != length(d) || length(Y) != nrow(X)) {
    print("Y, d, and rows of X must have the same length")
    exit <- 1
  }
  
  if (ncol(X) != 2) {
    print("X must have exactly 2 columns")
    exit <- 1
  }
  
  if (!(is.logical(d) || all(d %in% c(0, 1)))) {
    print("d must be a logical vector or a numeric vector containing only 0 and 1")
    exit <- 1
  }
  
  if (kernel!="gau" & kernel!="gaussian" & kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit <- 1
  }
  
  if (!is.numeric(deriv) | length(deriv) != 2) {
    print("deriv must be a numeric vector of length 2")
    exit <- 1
  } else if (sum(deriv) > p) {
    print("Sum of deriv components must be less than or equal to polynomial order p")
    exit <- 1
  }
  
  if (!is.null(tangvec)) {
    if (!(is.matrix(tangvec) || is.data.frame(tangvec)) ||
        nrow(tangvec) != nrow(b) || ncol(tangvec) != 2) {
      print("tangvec must be a matrix or data frame with the same number of rows as b and exactly 2 columns")
      exit <- 1
    }
  }
  
  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }
  
  if (exit>0) stop()
  
  # Data Cleaning
  
  dat <- cbind(X[,1], X[,2], Y, d)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("x.1", "x.2", "y", "d")
  eval <- as.data.frame(b)
  colnames(eval) <- c("x.1", "x.2")
  neval <- dim(eval)[1]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2) & complete.cases(dat$y) & complete.cases(dat$d)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  N.0 <- dim(dat[dat$d == 0,])[1]
  N.1 <- dim(dat[dat$d == 1,])[1]
  
  min_sample_size <- bwcheck
  if (is.null(bwcheck)) min_sample_size <- 50 + p + 1
  
  if (N < min_sample_size){
    warning("Not enough observations to perform bandwidth calculations.")
    stop()
  }
  
  if (is.null(p))         p <- 1
  kernel   <- tolower(kernel)
  
  e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
  deriv.sum <- deriv[1] + deriv[2]
  if (deriv.sum >= 1){
    e_deriv[,(factorial(deriv.sum+1)/(factorial(deriv.sum-1)* 2)) + deriv[2] + 1] <- 1
  } else {
    e_deriv[,1] <- 1
  }
  
  
  if (!is.null(tangvec)){
    warning("Tangvec provided. Ignore option deriv.")
    e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
    e_deriv[,2] <- tangvec[,1]
    e_deriv[,3] <- tangvec[,2]
    deriv <- c(1,0) # standardization for latter codes
    deriv.sum <- 1
  }
  
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"
  
  # Standardize data if necessary
  
  if (stdvars){
    sd.1 <- sd(dat$x.1)
    sd.2 <- sd(dat$x.2)
    dat$x.1 <- dat$x.1 / sd.1
    dat$x.2 <- dat$x.2 / sd.2
    eval$x.1 <- eval$x.1 / sd.1
    eval$x.2 <- eval$x.2 / sd.2
  } else {
    sd.1 <- 1
    sd.2 <- 1
  }
  
  # Store variance and bias constants for IMSE
  
  bconst <- rep(NA, neval)
  vconst <- rep(NA, neval)
  
  # Check for mass points
  
  M <- N; M.0 <- N.0; M.1 <- N.1
  if (masspoints == "check" | masspoints == "adjust"){
    unique.const <- rd2d_unique(dat)
    unique <- unique.const$unique
    M.0 <- dim(unique[unique$d == 0,])[1]
    M.1 <- dim(unique[unique$d == 1,])[1]
    M <- M.0 + M.1
    mass <- 1 - M / N
    if (mass >= 0.2){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & (masspoints == "check" | masspoints == "adjust")) bwcheck <- 50 + p + 1
    }
  }
  
  # Rule of thumb bandwidth selection
  
  dn <- rdbw2d_rot(dat,kernel.type, M)
  
  # Loop over points of evaluations
  results <- data.frame(matrix(NA, ncol = 16, nrow = neval))
  colnames(results) <- c('b1','b2','h01', 'h02', 'h11', 'h12', 'Nh.0', 'Nh.1',
                         'bias.0', 'bias.1', 'var.0', 'var.1','reg.bias.0','reg.bias.1','reg.var.0','reg.var.1')
  
  for (i in 1:neval){
    
    ev <- eval[i,]
    vec <- e_deriv[i,]
    
    # Center data
    
    dat.centered <- dat[,c("x.1", "x.2", "y", "d")]
    dat.centered$x.1 <- dat.centered$x.1 - ev$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev$x.2
    
    dat.centered$dist <- pmax(abs(dat.centered$x.1), abs(dat.centered$x.2)) # infinity norm
    if (kernel_type == "rad") dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2) # Euclidean norm
    
    if (masspoints == "adjust"){
      unique.centered <- unique
      unique.centered$x.1 <- unique.centered$x.1 - ev$x.1
      unique.centered$x.2 <- unique.centered$x.2 - ev$x.2
      
      unique.centered$dist <- pmax(abs(unique.centered$x.1), abs(unique.centered$x.2)) # infinity norm
      if (kernel_type == "rad") unique.centered$dist <- sqrt(unique.centered$x.1^2 + unique.centered$x.2^2) # Euclidean norm
    }
    
    # Weights
    
    dn.0 <- dn; dn.1 <- dn
    
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      if (masspoints == "adjust"){
        sorted.0 <- sort(unique.centered[unique.centered$d == FALSE,]$dist)
        sorted.1 <- sort(unique.centered[unique.centered$d == TRUE,]$dist)
        bw.min.0   <- sorted.0[min(bwcheck, M.0)]
        bw.min.1   <- sorted.1[min(bwcheck, M.1)]
        bw.max.0   <- sorted.0[length(sorted.0)]
        bw.max.1   <- sorted.1[length(sorted.1)]
      } else{
        sorted.0   <- sort(dat.centered[dat.centered$d == FALSE,]$dist)
        sorted.1   <- sort(dat.centered[dat.centered$d == TRUE,]$dist)
        bw.min.0   <- sorted.0[min(bwcheck, N.0)]
        bw.min.1   <- sorted.1[min(bwcheck, N.1)]
        bw.max.0   <- sorted.0[length(sorted.0)]
        bw.max.1   <- sorted.1[length(sorted.1)]
      }
      dn.0     <- max(dn.0, bw.min.0)
      dn.1     <- max(dn.1, bw.min.1)
      dn.0     <- min(dn.0, bw.max.0)
      dn.1     <- min(dn.1, bw.max.1)
    }
    
    if (kernel_type == "prod"){
      w.0 <- W.fun(dat.centered[dat.centered$d == FALSE,]$x.1/dn.0, kernel) *
        W.fun(dat.centered[dat.centered$d == FALSE,]$x.2/dn.0, kernel) / c(dn.0^2)
      w.1 <- W.fun(dat.centered[dat.centered$d == TRUE,]$x.1/dn.1, kernel) *
        W.fun(dat.centered[dat.centered$d == TRUE,]$x.2/dn.1, kernel) / c(dn.1^2)
    } else{
      w.0   <- W.fun(dat.centered[dat.centered$d == FALSE,]$dist/dn.0, kernel)/c(dn.0^2)
      w.1   <- W.fun(dat.centered[dat.centered$d == TRUE,]$dist/dn.1, kernel)/c(dn.1^2)
    }
    
    eN.0 <- sum(w.0 > 0)
    eN.1 <- sum(w.1 > 0)
    
    vec.q.0 <- get_coeff(dat.centered[dat.centered$d == FALSE,], vec, p, dn.0, kernel, kernel_type)
    vec.q.1 <- get_coeff(dat.centered[dat.centered$d == TRUE,], vec, p, dn.1, kernel, kernel_type)
    
    if (verbose) {print("Coefficients for a linear combination of (p+1)-th derivatives"); print(vec.q.0); print(vec.q.1)}
    
    # Bandwidth for fitting the linear combination of (p+1)-th derivatives using (p+1)-th degree model.
    
    thrshd.0 <- median(dat.centered[dat.centered$d == FALSE,]$dist)
    thrshd.1 <- median(dat.centered[dat.centered$d == TRUE,]$dist)
    
    bn.0 <- thrshd.0 # If method is "rot", use half of control data to estimate (p+1)th derivative of control.
    bn.1 <- thrshd.1 # If method is "rot", use half of treated data to estimate (p+1)th derivative of treated.
    
    if (method == "dpi"){
      
      bn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p + 1, vec.q.0, dn.0, thrshd.0, NULL, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      bn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p + 1, vec.q.1, dn.1, thrshd.1, NULL, vce, kernel, kernel_type,C[as.logical(dat.centered$d == TRUE)])
      
      bn.0 <-  ((2 + 2 * (p+1)) * bn.const.0$V  / ( (2 * (p + 1) + 2 - 2 * (p+1)) * (bn.const.0$B^2 + scaleregul * bn.const.0$Reg.1) ) )^(1/(2 * p + 6))
      bn.1 <-  ((2 + 2 * (p+1)) * bn.const.1$V  / ( (2 * (p + 1) + 2 - 2 * (p+1)) * (bn.const.1$B^2 + scaleregul * bn.const.1$Reg.1) ) )^(1/(2 * p + 6))
      
      if (verbose) print(paste("bn.0 = ", bn.0, ", bn.1 = ", bn.1, sep = ""))
      if (verbose) {print("Constants for bn.0 and bn.1:"); print(paste("B.0 = ", bn.const.0$B, ", V.0 = ", bn.const.0$V, ", Reg.0 = ", bn.const.0$Reg.1));
        print(paste("B.1 = ", bn.const.1$B, ", V.1 = ", bn.const.1$V, ", Reg.1 = ", bn.const.1$Reg.1))}
      
      if (!is.null(bwcheck)){ # Bandwidth restrictions
        bn.0 <- max(bn.0, bw.min.0)
        bn.1 <- max(bn.1, bw.min.1)
        bn.0 <- min(bn.0, bw.max.0)
        bn.1 <- min(bn.1, bw.max.1)
      }
    }
    
    # Bandwidth for estimating (derivatives of) treatment effect using p-th degree model.
    
    if (bwselect == "mserd" | bwselect == "imserd"){
      
      hn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p, vec, dn.0, bn.0, thrshd.0, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      hn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p, vec, dn.1, bn.1, thrshd.1, vce, kernel, kernel_type, C[as.logical(dat.centered$d == TRUE)])
      
      hn <- ( (2 + 2 * deriv.sum) * (hn.const.0$V   + hn.const.1$V) /
                ( (2 * p + 2 - 2 * deriv.sum) * ( (hn.const.0$B + scalebiascrct * hn.const.0$Reg.2 - hn.const.1$B - scalebiascrct * hn.const.1$Reg.2)^2 +
                                                    scaleregul * hn.const.0$Reg.1 + scaleregul * hn.const.1$Reg.1) ) )^(1/(2 * p + 4))
      
      if (!is.null(bwcheck)) { # Bandwidth restrictions
        hn     <- max(hn, bw.min.0, bw.min.1)
        hn     <- min(hn, max(bw.max.0, bw.max.1))
      }
      
      hn.0 <- hn.1 <- hn
    }
    
    if (bwselect == "msetwo" | bwselect == "imsetwo"){
      
      hn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p, vec, dn.0, bn.0, thrshd.0, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      hn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p, vec, dn.1, bn.1, thrshd.1, vce, kernel, kernel_type, C[as.logical(dat.centered$d == TRUE)])
      
      hn.0 <- ( (2 + 2 * deriv.sum) * hn.const.0$V /
                  ( (2 * p + 2) * ( (hn.const.0$B + scalebiascrct * hn.const.0$Reg.2)^2 + scaleregul * hn.const.0$Reg.1) ) )^(1/(2 * p + 4))
      hn.1 <- ( (2 + 2 * deriv.sum) * hn.const.1$V /
                  ( (2 * p + 2) * ( (hn.const.1$B + scalebiascrct * hn.const.1$Reg.2)^2 + scaleregul * hn.const.1$Reg.1) ) )^(1/(2 * p + 4))
      
      if (!is.null(bwcheck)) { # Bandwidth restrictions
        hn.0     <- max(hn.0, bw.min.0)
        hn.1     <- max(hn.1, bw.min.1)
        hn.0     <- min(hn.0, bw.max.0)
        hn.1     <- min(hn.1, bw.max.1)
      }
    }
    
    results[i,c(1:2)] <- c(ev$x.1, ev$x.2)
    results[i,c(3:16)] <- c(hn.0, hn.0, hn.1, hn.1, eN.0, eN.1, hn.const.0$B, hn.const.1$B, hn.const.0$V, hn.const.1$V, hn.const.0$Reg.2, hn.const.1$Reg.2,
                            hn.const.0$Reg.1, hn.const.1$Reg.1)
  }
  
  if (bwselect == "imserd"){
    V.V <- mean(results$var.0) + mean(results$var.1)
    B.B <- mean( (results$bias.0 + scalebiascrct * results$reg.bias.0 - results$bias.1 - scalebiascrct * results$reg.bias.1)^2 + scaleregul * results$reg.var.0 + scaleregul * results$reg.var.1 )
    hIMSE <- ((2 + 2 * deriv.sum) * V.V / ( (2 * p + 2) * B.B ) )^(1/(2 * p + 4))
    results$h01 <- rep(hIMSE, dim(results)[1])
    results$h02 <- rep(hIMSE, dim(results)[1])
    results$h11 <- rep(hIMSE, dim(results)[1])
    results$h12 <- rep(hIMSE, dim(results)[1])
  }
  
  if (bwselect == "imsetwo"){
    V.V.0 <- mean(results$var.0)
    V.V.1 <- mean(results$var.1)
    B.B.0 <- mean( (results$bias.0 + scalebiascrct * results$reg.bias.0)^2 + scaleregul * results$reg.var.0)
    B.B.1 <- mean( (results$bias.1 + scalebiascrct * results$reg.bias.1)^2 + scaleregul * results$reg.var.1)
    hIMSE.0 <- ((2 + 2 * deriv.sum) * V.V.0 / ( (2 * p + 2) * B.B.0 ) )^(1/(2 * p + 4))
    hIMSE.1 <- ((2 + 2 * deriv.sum) * V.V.1 / ( (2 * p + 2) * B.B.1 ) )^(1/(2 * p + 4))
    results$h01 <- rep(hIMSE.0, dim(results)[1])
    results$h02 <- rep(hIMSE.0, dim(results)[1])
    results$h11 <- rep(hIMSE.1, dim(results)[1])
    results$h12 <- rep(hIMSE.1, dim(results)[1])
  }
  
  # Standardization (sd.1 = sd.2 = 1 if stdvar == FALSE)
  
  results$h01 <- results$h01 * sd.1
  results$h02 <- results$h02 * sd.2
  results$h11 <- results$h11 * sd.1
  results$h12 <- results$h12 * sd.2
  
  # Outputs
  
  bws <- results[,c("h01", "h02", "h11", "h12")]
  bws <- cbind(eval, bws)
  colnames(bws) <- c("b1","b2","h01", "h02", "h11", "h12")
  
  clustered <- !is.null(C)
  
  out        <- list(bws = bws, mseconsts = results,
                     opt = list(N=N, N.0 = N.0, N.1 = N.1,M.0 = M.0,
                                M.1 = M.1, neval=neval, p=p, deriv=deriv, tangvec = tangvec,
                                kernel=kernel.type, kernel_type = kernel_type,
                                bwselect=bwselect, method = method, bwcheck = bwcheck,
                                stdvars = stdvars, C= C, clustered = clustered,
                                vce = vce, masspoints = masspoints,
                                scaleregul = scaleregul, scalebiascrct = scalebiascrct))
  out$call   <- match.call()
  class(out) <- "rdbw2d"
  return(out)
}

sqrtm <- function(x) {
  ## Generate Basic informations of matrix x
  ## FIXME : should work for "Matrix" too, hence _not_  S <- as.matrix(x)
  d <- dim(x)
  if(length(d) != 2 || d[1] != d[2])
    stop(gettextf("'%s' must be a square matrix", "x"), domain=NA)
  
  ##MM: No need to really check here; we get correct error msg later anyway
  ##	  and don't need to compute det() here, in the good cases !
  ##	  if (det(x) == 0) stop("'x' is singular")
  n <- d[1]
  
  ##------- STEP 0: Schur Decomposition ---------------------------------------
  
  Sch.x <- Schur(Matrix(x)) ## <- {FIXME [Matrix]}
  ev <- Sch.x@EValues
  if(getOption("verbose") && any(abs(Arg(ev) - pi) < 1e-7))
    ## Let's see what works: temporarily *NOT* stop()ping :
    message(gettextf("'x' has negative real eigenvalues; maybe ok for %s", "sqrtm()"),
            domain=NA)
  
  S <- as.matrix(Sch.x@T)
  Q <- as.matrix(Sch.x@Q)
  
  ##---------STEP 1: Analyse block structure-----------------------------------
  if(n > 1L) {
    ## Count 2x2 blocks (as Schur(x) is the real Schur Decompostion)
    J.has.2 <- S[cbind(2:n, 1:(n-1))] != 0
    k <- sum(J.has.2) ## := number of non-zero SUB-diagonals
  } else k <- 0L
  
  ## Generate Blockstructure and save it as R.index
  R.index <- vector("list",n-k)
  l <- 1L
  i <- 1L
  while(i < n) { ## i advances by 1 or 2, depending on 1- or 2- Jordan Block
    if (S[i+1L,i] == 0) {
      R.index[[l]] <- i
    }
    else {
      i1 <- i+1L
      R.index[[l]] <- c(i,i1) # = i:(i+1)
      i <- i1
    }
    i <- i+1L
    l <- l+1L
  }
  if (is.null(R.index[[n-k]])) { # needed; FIXME: should be able to "know"
    ##message(gettextf("R.index[n-k = %d]] is NULL, set to n=%d", n-k,n), domain=NA)
    R.index[[n-k]] <- n
  }
  
  ##---------STEP 2: Calculate diagonal elements/blocks------------------------
  ## Calculate the root of the diagonal blocks of the Schur Decompostion S
  I <- diag(2)
  X <- matrix(0,n,n)
  for (j in seq_len(n-k)) {
    ij <- R.index[[j]]
    if (length(ij) == 1L) {
      X[ij,ij] <- if((.s <- S[ij,ij]) < 0) sqrt(.s + 0i) else sqrt(.s)
    }
    else {
      ev1 <- ev[ij[1]]
      r1 <- Re(sqrt(ev1)) ## sqrt(<complex>) ...
      X[ij,ij] <- r1*I + 1/(2*r1)*(S[ij,ij] - Re(ev1)*I)
    }
  }
  ##---------STEP 3: Calculate superdiagonal elements/blocks-------------------
  
  ## Calculate the remaining, not-diagonal blocks
  if (n-k > 1L) for (j in 2L:(n-k)) {
    ij <- R.index[[j]]
    for (i in (j-1L):1L) {
      ii <- R.index[[i]]
      sumU <- 0
      
      ## Calculation for 1x1 Blocks
      if (length(ij) == 1L & length(ii) == 1L) {
        if (j-i > 1L) for (l in (i+1L):(j-1L)) {
          il <- R.index[[l]]
          sumU <- sumU + {
            if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
            else		      X[ii,il] * X[il,ij]
          }
        }
        X[ii,ij] <- solve(X[ii,ii]+X[ij,ij],S[ii,ij]-sumU)
      }
      
      ## Calculation for	1x2 Blocks
      else if (length(ij) == 2 & length(ii) == 1L ) {
        if (j-i > 1L) for (l in(i+1L):(j-1L)) {
          il <- R.index[[l]]
          sumU <- sumU + {
            if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
            else		      X[ii,il] * X[il,ij]
          }
        }
        X[ii,ij] <- solve(t(X[ii,ii]*I + X[ij,ij]),
                          as.vector(S[ii,ij] - sumU))
      }
      ## Calculation for	2x1 Blocks
      else if (length(ij) == 1L & length(ii) == 2 ) {
        if (j-i > 1L) for (l in(i+1L):(j-1L)) {
          il <- R.index[[l]]
          sumU <- sumU + {
            if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
            else		      X[ii,il] * X[il,ij]
          }
        }
        X[ii,ij] <- solve(X[ii,ii]+X[ij,ij]*I, S[ii,ij]-sumU)
      }
      ## Calculation for	2x2 Blocks with special equation for solver
      else if (length(ij) == 2 & length(ii) == 2 ) {
        if (j-i > 1L) for (l in(i+1L):(j-1L)) {
          il <- R.index[[l]]
          sumU <- sumU + {
            if (length(il) == 2 ) X[ii,il] %*%  X[il,ij]
            else		      X[ii,il] %*% t(X[il,ij])
            
          }
        }
        tUii <- matrix(0,4,4)
        tUii[1:2,1:2] <- X[ii,ii]
        tUii[3:4,3:4] <- X[ii,ii]
        tUjj <- matrix(0,4,4)
        tUjj[1:2,1:2] <- t(X[ij,ij])[1L,1L]*I
        tUjj[3:4,3:4] <- t(X[ij,ij])[2L,2L]*I
        tUjj[1:2,3:4] <- t(X[ij,ij])[1L,2L]*I
        tUjj[3:4,1:2] <- t(X[ij,ij])[2L,1L]*I
        X[ii,ij] <- solve(tUii+tUjj, as.vector(S[ii,ij]-sumU))
      }
    } ## for (i in (j-1):1) ..
  } ## for (j in 2:(n-k)) ...
  
  ##------- STEP 4: Reverse the Schur Decomposition --------------------------
  ## Reverse the Schur Decomposition
  Q %*% X %*% solve(Q)
}

ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  #
  # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
  #
  if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if(!is.matrix(X)) X <- as.matrix(X)
  Xsvd <- svd(X)
  if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if(!any(Positive)) array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
}

mvrnorm <-
  function(n = 1, mu, Sigma, tol=1e-6, empirical = FALSE, EISPACK = FALSE)
  {
    p <- length(mu)
    if(!all(dim(Sigma) == c(p,p))) stop("incompatible arguments")
    if(EISPACK) stop("'EISPACK' is no longer supported by R", domain = NA)
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if(!all(ev >= -tol*abs(ev[1L]))) stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if(empirical) {
      X <- scale(X, TRUE, FALSE) # remove means
      X <- X %*% svd(X, nu = 0)$v # rotate to PCs
      X <- scale(X, FALSE, TRUE) # rescale PCs to unit variance
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    nm <- names(mu)
    if(is.null(nm) && !is.null(dn <- dimnames(Sigma))) nm <- dn[[1L]]
    dimnames(X) <- list(nm, NULL)
    if(n == 1) drop(X) else t(X)
  }

rd2d <- function(Y, X, t, b, h = NULL, deriv = c(0,0), tangvec = NULL,
                 p = 1, q = 2, kernel = c("tri","triangular","epa","epanechnikov","uni","uniform","gau","gaussian"),
                 kernel_type = c("prod","rad"), vce = c("hc1","hc0","hc2","hc3"),
                 masspoints = c("check", "adjust", "off"),C = NULL,
                 level = 95, cbands = TRUE, side = c("two", "left", "right"), repp = 1000,
                 bwselect = c("mserd", "imserd", "msetwo", "imsetwo", "user provided"),
                 method = c("dpi", "rot"), bwcheck = 50 + p + 1,
                 scaleregul = 3, scalebiascrct = 1, stdvars = TRUE){
  
  ######################## Input error handling ################################
  
  kernel <- match.arg(kernel)
  kernel_type <- match.arg(kernel_type)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  side <- match.arg(side)
  bwselect <- match.arg(bwselect)
  method <- match.arg(method)
  
  d <- t # renaming the variable
  
  exit <- 0
  
  # Check Y, X, d lengths
  if (length(Y) != length(d) || length(Y) != nrow(X)) {
    print("Y, d, and rows of X must have the same length")
    exit <- 1
  }
  
  # X must have 2 columns
  if (ncol(X) != 2) {
    print("X must have exactly 2 columns")
    exit <- 1
  }
  
  # d must be logical or contain only 0 and 1
  if (!(is.logical(d) || all(d %in% c(0, 1)))) {
    print("d must be a logical vector or a numeric vector containing only 0 and 1")
    exit <- 1
  }
  
  # b must be a matrix with 2 columns
  if (!(is.matrix(b) || is.data.frame(b)) || ncol(b) != 2) {
    print("b must be a matrix with 2 columns")
    exit <- 1
  }
  
  # h must be either a positive scalar or a matrix/data.frame with same rows as b and 4 columns
  if (!is.null(h)) {
    if (length(h) == 1) {
      if (!is.numeric(h) || h <= 0) {
        print("If h is a scalar, it must be a positive numeric value")
        exit <- 1
      }
    } else if (!(is.matrix(h) || is.data.frame(h)) ||
               nrow(h) != nrow(b) || ncol(h) != 4) {
      print("If h is not a scalar, it must be a matrix or data frame with the same number of rows as b and 4 columns")
      exit <- 1
    }
  }
  
  # deriv must be a numeric vector of length 2, and deriv[1] + deriv[2] <= p
  if (!is.numeric(deriv) || length(deriv) != 2) {
    print("deriv must be a numeric vector of length 2")
    exit <- 1
  } else if (sum(deriv) > p) {
    print("Sum of deriv components must be less than or equal to polynomial order p")
    exit <- 1
  }
  
  # tangvec, if provided, must be matrix/data.frame with same nrow as b and 2 columns
  if (!is.null(tangvec)) {
    if (!(is.matrix(tangvec) || is.data.frame(tangvec)) ||
        nrow(tangvec) != nrow(b) || ncol(tangvec) != 2) {
      warning("tangvec must be a matrix or data frame with same number of rows as b and 2 columns")
      exit <- 1
    }
  }
  
  # level must be numeric in (0, 100)
  if (!is.numeric(level) || level <= 0 || level >= 100) {
    print("level must be a numeric value between 0 and 100")
    exit <- 1
  }
  
  # repp must be a positive integer
  if (!is.numeric(repp) || repp < 1 || repp != as.integer(repp)) {
    print("repp must be a positive integer")
    exit <- 1
  }
  
  if (q < p){
    print("Parameter q must be no smaller than p. Please provide valid inputs.")
    exit <- 1
  }
  
  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }
  
  if (exit>0) stop()
  
  ############################ Data preparation ################################
  
  dat <- cbind(X[,1], X[,2], Y, d)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("x.1", "x.2", "y", "d")
  eval <- as.data.frame(b)
  colnames(eval) <- c("x.1", "x.2")
  neval <- dim(eval)[1]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2) & complete.cases(dat$y) & complete.cases(dat$d)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  N.0 <- dim(dat[dat$d == 0,])[1]
  N.1 <- dim(dat[dat$d == 1,])[1]
  
  if (is.null(p))         p <- 1
  kernel   <- tolower(kernel)
  
  e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
  deriv.sum <- deriv[1] + deriv[2]
  if (deriv.sum >= 1){
    e_deriv[,(factorial(deriv.sum+1)/(factorial(deriv.sum-1)* 2)) + deriv[2] + 1] <- 1
  } else {
    e_deriv[,1] <- 1
  }
  
  if (!is.null(tangvec)){
    warning("Tangvec provided. Ignore option deriv.")
    e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
    e_deriv[,2] <- tangvec[,1]
    e_deriv[,3] <- tangvec[,2]
    deriv <- c(1,0) # standardization for latter codes
    deriv.sum <- 1
  }
  
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"
  
  # Check for mass points
  
  M <- NULL; M.0 <- NULL; M.1 <- NULL
  
  
  # Only check for mass points if inputting bivariate coordinates.
  if (masspoints == "check" | masspoints == "adjust"){
    unique.const <- rd2d_unique(dat)
    unique <- unique.const$unique
    M.0 <- dim(unique[unique$d == 0,])[1]
    M.1 <- dim(unique[unique$d == 1,])[1]
    M <- M.0 + M.1
    mass <- 1 - M / N
    if (mass >= 0.2){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & (masspoints == "check" | masspoints == "adjust")) bwcheck <- 50 + p + 1
    }
  }
  
  min_sample_size <- bwcheck
  if (is.null(bwcheck)) min_sample_size <- 50 + p + 1
  
  if (N < min_sample_size){
    warning("Not enough observations to perform RDD calculations. ")
    stop()
  }
  
  ################################ Bandwidth ###################################
  
  if (is.null(h)){
    bws <- rdbw2d(Y = Y, X = X, t = d, b = b, p = p, deriv = deriv, tangvec = tangvec,
                  kernel = kernel, kernel_type = kernel_type,
                  bwselect = bwselect, method = method, vce = vce,
                  bwcheck = bwcheck, masspoints = masspoints,
                  C = C, scaleregul = scaleregul, scalebiascrct = scalebiascrct,
                  stdvars = stdvars)
    bws <- bws$bws
    hgrid <- cbind(bws[,3],bws[,4])
    hgrid.1 <- cbind(bws[,5],bws[,6])
  } else {
    bwselect <- "user provided"
    # standardize bandwidth
    if (length(h) == 1){
      hgrid <- matrix(h, nrow = neval, ncol = 2)
      hgrid.1 <- matrix(h, nrow = neval, ncol = 2)
    } else {
      hgrid <- cbind(h[,1],h[,2])
      hgrid.1 <- cbind(h[,3],h[,4])
    }
  }
  
  ###################### Point estimation and inference ########################
  
  count.q <- factorial(q + 2)/(factorial(q) * 2)
  count.p <- factorial(p + 2)/(factorial(p) * 2)
  
  e_deriv.q <- matrix(0, nrow = neval, ncol = count.q)
  e_deriv.q[,c(1:count.p)] <- e_deriv
  
  rdfit.p <- rd2d_fit_v2(dat, eval, e_deriv, deriv, p, hgrid, hgrid.1, kernel, kernel_type, vce, masspoints, C, bwcheck, unique)
  tau.hat.p <- rdfit.p$mu.1 - rdfit.p$mu.0
  se.hat.p <- sqrt(rdfit.p$se.0^2 + rdfit.p$se.1^2)
  h.01.p <- rdfit.p$h.0.x
  h.02.p <- rdfit.p$h.0.y
  h.11.p <- rdfit.p$h.1.x
  h.12.p <- rdfit.p$h.1.y
  eN.0.p <- rdfit.p$eN.0
  eN.1.p <- rdfit.p$eN.1
  
  rdfit.q <- rd2d_fit_v2(dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel, kernel_type, vce, masspoints, C, bwcheck, unique)
  tau.hat.q <- rdfit.q$mu.1 - rdfit.q$mu.0
  se.hat.q <- sqrt(rdfit.q$se.0^2 + rdfit.q$se.1^2)
  h.01.q <- rdfit.q$h.0.x
  h.02.q <- rdfit.q$h.0.y
  h.11.q <- rdfit.q$h.1.x
  h.12.q <- rdfit.q$h.1.y
  eN.0.q <- rdfit.q$eN.0
  eN.1.q <- rdfit.q$eN.1
  
  zvalues <- tau.hat.q/se.hat.q
  pvalues <- 2 * pnorm(abs(zvalues),lower.tail = FALSE)
  
  if (side == "two"){
    zval <- qnorm((level + 100)/ 200)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- tau.hat.q + zval * se.hat.q
  }
  if (side == "left"){
    zval <- qnorm(level / 100)
    CI.upper <- tau.hat.q + zval * se.hat.q
    CI.lower <- rep(-Inf, length(CI.upper))
  }
  if (side == "right"){
    zval <- qnorm(level / 100)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- rep(Inf, length(CI.lower))
  }
  
  # Covariance
  
  cov.hat.q <- NA
  cb.hat.q <- list(CI.l = CI.lower, CI.r = CI.upper, CB.l = rep(NA, length(CI.lower)), CB.r = rep(NA, length(CI.lower)))
  CB.lower <- NA
  CB.upper <- NA
  if (cbands){
    cov.hat.q <- rdbw2d_cov(dat, eval, e_deriv.q, deriv, q, cbind(h.01.q, h.02.q), cbind(h.11.q, h.12.q), kernel, kernel_type, vce, C)
    cb.hat.q <- rd2d_cb(tau.hat.q, cov.hat.q, repp, side, level)
    CB.lower <- cb.hat.q$CB.l
    CB.upper <- cb.hat.q$CB.r
  }
  
  clustered <- !is.null(C)
  
  ################################## Output ####################################
  
  main <- cbind(b[,1], b[,2], tau.hat.p, se.hat.p, tau.hat.q, se.hat.q, zvalues, pvalues,
                CI.lower, CI.upper, CB.lower, CB.upper, hgrid[,1], hgrid[,2],
                hgrid.1[,1], hgrid.1[,2], eN.0.p, eN.1.p)
  main <- as.data.frame(main)
  colnames(main) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q", "z", "P>|z|",
                      "CI.lower","CI.upper","CB.lower", "CB.upper", "h01", "h02",
                      "h11", "h12", "Nh0", "Nh1")
  
  main.A0 <- cbind(b[,1], b[,2], rdfit.p$mu.0, rdfit.p$se.0, rdfit.q$mu.0,rdfit.q$se.0, hgrid[,1], hgrid[,2], eN.0.p)
  main.A0 <- as.data.frame(main.A0)
  colnames(main.A0) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h01", "h02","Nh0")
  
  main.A1 <- cbind(b[,1], b[,2], rdfit.p$mu.1, rdfit.p$se.1, rdfit.q$mu.1,rdfit.q$se.1, hgrid.1[,1], hgrid.1[,2], eN.1.p)
  main.A1 <- as.data.frame(main.A1)
  colnames(main.A1) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h11", "h12","Nh1")
  
  rdmodel <- "rd2d"
  
  out <- list(results = main, results.A0 = main.A0, results.A1 = main.A1,
              opt=list(b = b, deriv = deriv, tangvec = tangvec, p = p, q = q, kernel=kernel.type, kernel_type = kernel_type, N=N, N.0 = N.0,
                       N.1 = N.1, M = M, M.0 = M.0, M.1 = M.1, neval=neval, bwselect = bwselect, method = method,
                       vce = vce, bwcheck = bwcheck, masspoints = masspoints, C = C, clustered = clustered,
                       scaleregul = scaleregul, scalebiascrct = scalebiascrct, stdvars = stdvars,
                       level = level, repp = repp, side = side,cbands = cbands,
                       h01 = hgrid[,1], h02 = hgrid[,2], h11 = hgrid.1[,1], h12 = hgrid.1[,2],
                       Nh0 = eN.0.p, Nh1 = eN.1.p), cov.q=cov.hat.q, rdmodel = rdmodel)
  out$call   <- match.call()
  class(out) <- "rd2d"
  
  return(out)
}

rd2d_fit_v2 <- function(dat, eval, deriv = NULL,o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                        kernel = "epa", kernel_type = "prod", vce = "hc1",
                        masspoints = "adjust", C = NULL, bwcheck = 50 + p + 1, unique = NULL){
  
  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  
  sd.x1 <- sd(dat$x.1)
  sd.x2 <- sd(dat$x.2)
  
  neval <- dim(eval)[1]
  
  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.data.frame(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.data.frame(hgrid.1)
  }
  
  if (ncol(hgrid.0) == 1){
    results <- data.frame(matrix(NA, ncol = 10, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0', 'h.1', 'mu.0', 'mu.1', 'se.0',
                           'se.1', 'eN.0', 'eN.1')
  } else {
    results <- data.frame(matrix(NA, ncol = 12, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0.x', 'h.0.y', 'h.1.x', 'h.1.y',
                           'mu.0', 'mu.1', 'se.0', 'se.1', 'eN.0', 'eN.1')
  }
  
  # Check for compatibility
  
  for (i in 1:neval){
    
    ev <- eval[i,]
    vec <- deriv[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]
    
    # Center data
    
    dat.centered <- dat[,c("x.1", "x.2", "y", "d")]
    dat.centered$x.1 <- dat.centered$x.1 - ev$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev$x.2
    # for product kernel, standardize the covariates so that they have the same sd
    dat.centered$dist <- pmax(abs(dat.centered$x.1/sd.x1), abs(dat.centered$x.2/sd.x2)) # infinity norm
    if (kernel_type == "rad") dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2) # Euclidean norm
    
    # Bandwidth restriction
    
    if (masspoints == "adjust"){
      unique.centered <- unique
      unique.centered$x.1 <- unique.centered$x.1 - ev$x.1
      unique.centered$x.2 <- unique.centered$x.2 - ev$x.2
      unique.centered$dist <- pmax(abs(unique.centered$x.1/sd.x1), abs(unique.centered$x.2/sd.x2)) # infinity norm
      if (kernel_type == "rad") unique.centered$dist <- sqrt(unique.centered$x.1^2 + unique.centered$x.2^2) # Euclidean norm
    }
    
    if (!is.null(bwcheck)){
      
      if (masspoints == "adjust"){
        sorted.0 <- sort(unique.centered[unique.centered$d == FALSE,]$dist)
        sorted.1 <- sort(unique.centered[unique.centered$d == TRUE,]$dist)
      } else{
        sorted.0   <- sort(dat.centered[dat.centered$d == FALSE,]$dist)
        sorted.1   <- sort(dat.centered[dat.centered$d == TRUE,]$dist)
      }
      
      bw.min.0   <- sorted.0[bwcheck]
      bw.min.1   <- sorted.1[bwcheck]
      bw.max.0   <- sorted.0[length(sorted.0)]
      bw.max.1   <- sorted.1[length(sorted.1)]
      
      # convert to the original if using product kernel
      if (kernel_type == "prod"){
        multiplier <- c(sd.x1, sd.x2)
      } else {
        multiplier <- c(1,1)
      }
      
      bw.min.0 <- bw.min.0 * multiplier
      bw.min.1 <- bw.min.1 * multiplier
      bw.max.0 <- bw.max.0 * multiplier
      bw.max.1 <- bw.max.1 * multiplier
      
      if (!is.null(hgrid.1)){
        h.0     <- pmax(h.0, bw.min.0)
        h.1     <- pmax(h.1, bw.min.1)
        h.0     <- pmin(h.0, bw.max.0)
        h.1     <- pmin(h.1, bw.max.1)
      } else{
        h.0 <- pmax(h.0, bw.min.0,bw.min.1)
        h.0 <- pmin(h.0, pmax(bw.max.0, bw.max.1))
        h.1 <- h.0
      }
    }
    
    fit.0.p <- rd2d_lm(dat.centered[dat.centered$d == 0,], h.0, p, vce, kernel = kernel, C = C[dat.centered$d == 0],
                       varr = TRUE, kernel_type = kernel_type)
    fit.1.p <- rd2d_lm(dat.centered[dat.centered$d == 1,], h.1, p, vce, kernel = kernel, C = C[dat.centered$d == 1],
                       varr = TRUE, kernel_type = kernel_type)
    mu.0 <- (vec %*% fit.0.p$beta)[1,1]
    mu.1 <- (vec %*% fit.1.p$beta)[1,1]
    
    # standardize
    if (length(h.0) == 1){
      h.0.x <- as.numeric(h.0)
      h.0.y <- as.numeric(h.0)
    }  else {
      h.0.x <- as.numeric(h.0[1])
      h.0.y <- as.numeric(h.0[2])
    }
    if (length(h.1) == 1){
      h.1.x <- as.numeric(h.1)
      h.1.y <- as.numeric(h.1)
    }  else {
      h.1.x <- as.numeric(h.1[1])
      h.1.y <- as.numeric(h.1[2])
    }
    
    # standard deviation
    invH.0 <- get_invH(c(h.0.x, h.0.y),p)
    se.0 <- matrix(vec, nrow = 1) %*% invH.0 %*% fit.0.p$cov.const %*% invH.0 %*% matrix(vec, ncol = 1) / (h.0.x * h.0.y)
    # se.0 <- matrix(vec, nrow = 1) %*% fit.0.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.0^(2 + 2 * o))
    se.0 <- sqrt(se.0[1,1])
    invH.1 <- get_invH(c(h.1.x, h.1.y),p)
    se.1 <- matrix(vec, nrow = 1) %*% invH.1 %*% fit.1.p$cov.const %*% invH.1 %*% matrix(vec, ncol = 1) / (h.1.x * h.1.y)
    # se.1 <- matrix(vec, nrow = 1) %*% fit.1.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.1^(2 + 2 * o))
    se.1 <- sqrt(se.1[1,1])
    
    # effective sample size
    eN.0 <- fit.0.p$eN
    eN.1 <- fit.1.p$eN
    
    if (ncol(hgrid.0) == 1){
      results[i,] <- c(ev[1], ev[2], h.0, h.1, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    } else {
      results[i,] <- c(ev[1], ev[2], h.0.x, h.0.y, h.1.x, h.1.y, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    }
  }
  
  return(results)
}

rdbw2d_cov <- function(dat, eval, deriv = NULL, o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                       kernel = "epa", kernel_type = "prod", vce = "hc2", C = NULL){
  
  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  
  neval <- dim(eval)[1]
  covs <- matrix(NA, nrow = neval, ncol = neval)
  halves.0 <- list()
  halves.1 <- list()
  inds.0 <- list()
  inds.1 <- list()
  
  clusters <- unique(C)
  g <- length(clusters)
  
  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }
  
  for (i in 1:neval){
    
    ev.a <- eval[i,]
    h.a.0 <- hgrid.0[i,]
    h.a.1 <- h.a.0
    if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
    deriv.vec.a <- deriv[i,]
    
    dat.centered <- dat[, c("x.1", "x.2", "y", "d")]
    
    dat.centered$x.1 <- dat.centered$x.1 - ev.a$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev.a$x.2
    dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)
    
    C.0 <- C[as.logical(dat.centered$d == 0)]
    C.1 <- C[as.logical(dat.centered$d == 1)]
    
    cov.half.consts.0 <- get_cov_half_v2(dat.centered[dat.centered$d == 0,], h.a.0, p, vce, kernel, kernel_type, C.0, clusters)
    cov.half.consts.1 <- get_cov_half_v2(dat.centered[dat.centered$d == 1,], h.a.1, p, vce, kernel, kernel_type, C.1, clusters)
    half.const.0 <- cov.half.consts.0$cov.half.const
    half.const.1 <- cov.half.consts.1$cov.half.const
    ind.0 <- cov.half.consts.0$ind
    ind.1 <- cov.half.consts.1$ind
    
    inds.0[[i]] <- ind.0
    inds.1[[i]] <- ind.1
    halves.0[[i]] <- half.const.0
    halves.1[[i]] <- half.const.1
    
  }
  
  for (i in 1:neval){
    for (j in i:neval){
      
      h.a.0 <- hgrid.0[i,]
      h.a.1 <- h.a.0
      if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
      h.b.0 <- hgrid.0[j,]
      h.b.1 <- h.b.0
      if (!is.null(hgrid.1)) h.b.1 <- hgrid.1[j,]
      
      if (length(h.a.0) == 1){
        h.a.0.x <- h.a.0
        h.a.0.y <- h.a.0
        h.a.1.x <- h.a.1
        h.a.1.y <- h.a.1
        h.b.0.x <- h.b.0
        h.b.0.y <- h.b.0
        h.b.1.x <- h.b.1
        h.b.1.y <- h.b.1
      } else{
        h.a.0.x <- h.a.0[1]
        h.a.0.y <- h.a.0[2]
        h.a.1.x <- h.a.1[1]
        h.a.1.y <- h.a.1[2]
        h.b.0.x <- h.b.0[1]
        h.b.0.y <- h.b.0[2]
        h.b.1.x <- h.b.1[1]
        h.b.1.y <- h.b.1[2]
      }
      deriv.vec.a <- deriv[i,]
      deriv.vec.b <- deriv[j,]
      
      half.0.a <- halves.0[[i]]
      half.0.b <- halves.0[[j]]
      half.1.a <- halves.1[[i]]
      half.1.b <- halves.1[[j]]
      
      invH.a.0 <- get_invH(c(h.a.0.x,h.a.0.y),p)
      invH.a.1 <- get_invH(c(h.a.1.x,h.a.1.y),p)
      invH.b.0 <- get_invH(c(h.b.0.x,h.b.0.y),p)
      invH.b.1 <- get_invH(c(h.b.1.x,h.b.1.y),p)
      
      if (is.null(C)){
        
        ind.0.a <- inds.0[[i]]
        ind.0.b <- inds.0[[j]]
        ind.1.a <- inds.1[[i]]
        ind.1.b <- inds.1[[j]]
        
        ind.0 <- ind.0.a * ind.0.b
        ind.0.a <- as.logical(ind.0[ind.0.a])
        ind.0.b <- as.logical(ind.0[ind.0.b])
        ind.1 <- ind.1.a * ind.1.b
        ind.1.a <- as.logical(ind.1[ind.1.a])
        ind.1.b <- as.logical(ind.1[ind.1.b])
        
        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% t(half.0.a[ind.0.a,,drop = "FALSE"]) %*% half.0.b[ind.0.b,,drop = "FALSE"] %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] / (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        # cov.0 <- cov.0[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% t(half.1.a[ind.1.a,,drop = "FALSE"]) %*% half.1.b[ind.1.b,,drop = "FALSE"] %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] / (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
        # cov.1 <- cov.1[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
      }
      
      if (!is.null(C)){
        k <- dim(deriv)[2]
        M.0 <- matrix(0, k, k)
        M.1 <- matrix(0, k, k)
        for (l in 1:g){
          M.0 <- M.0 + matrix(half.0.a[l,], ncol = 1) %*% matrix(half.0.b[l,], nrow = 1)
          M.1 <- M.1 + matrix(half.1.a[l,], ncol = 1) %*% matrix(half.1.b[l,], nrow = 1)
        }
        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% M.0 %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] /  (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% M.1 %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] /  (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
      }
      
      covs[i,j] <- cov.0 + cov.1
      covs[j,i] <- cov.0 + cov.1
    }
  }
  return(covs)
}

infl <- function(x, invG){
  result <- t(matrix(x, ncol = 1)) %*% invG %*% matrix(x, ncol = 1)
  return(result[1,1])
}

rd2d_lm <- function(dat, h, p, vce = "hc1", kernel = "epa", kernel_type = "prod",
                    C = NULL, varr = FALSE){
  
  dat <- dat[,c("x.1", "x.2", "y", "d", "dist")]
  
  # Variance and coefficients for a linear combination of (p+1)-th derivatives.
  
  h <- as.vector(as.matrix(h)) # if h is data frame, convert it to a vector
  
  # checks
  if (kernel_type == "prod"){ # product kernel
    if (length(h) == 1){
      h <- c(h,h)
    }
  }
  else{ # radius kernel
    if (length(h) == 2){
      h <- sqrt(h[1]^2 + h[2]^2)
    }
  }
  
  # weights
  if (kernel_type == "prod"){
    w <- W.fun(dat$x.1/c(h[1]), kernel) * W.fun(dat$x.2/c(h[2]), kernel) / c(h[1] * h[2])
  }
  else{
    w <- W.fun(dat$dist/c(h), kernel)/c(h^2)
  }
  
  if (length(h) == 1){
    h.x <- h; h.y <- h
  } else {
    h.x <- h[1]; h.y <- h[2]
  }
  
  ind <- as.logical(w > 0)
  
  eN <- sum(ind)
  
  ew <- w[ind]
  eY <- dat$y[ind]
  eC <- C[ind] # if C == NULL, eC == NULL.
  
  eu <- dat[ind, c("x.1", "x.2")]
  eu$x.1 <- eu$x.1/h.x
  eu$x.2 <- eu$x.2/h.y
  
  eR <- as.matrix(get_basis(eu,p))
  
  sqrtw_R <- sqrt(ew) * eR
  sqrtw_Y <- sqrt(ew) * eY
  
  w_R <- ew * eR
  
  invG <- qrXXinv(sqrtw_R)
  
  invH.p <- get_invH(c(h.x,h.y),p)
  
  H.p <- get_H(c(h.x,h.y),p)
  
  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)
  
  cov.const <- NA
  
  if (varr){
    
    resd <- (eY - (eR %*% H.p) %*% beta)[,1]
    
    lambda <- function(x){infl(x, invG)}
    
    if (vce=="hc0") {
      w.vce = 1
    } else if (vce=="hc1") {
      w.vce = sqrt(eN/(eN-factorial(p+2)/(factorial(p) * 2)))
    } else if (vce=="hc2") {
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = 1/(1-hii)
    }
    
    resd <- resd * w.vce
    
    sigma <- rd2d_vce(w_R, resd, eC, c(h.x, h.y))
    
    # sigma <-  t(resd * as.matrix(sqrtw_R)) %*% (ew * resd * as.matrix(sqrtw_R)) * h^2
    
    cov.const <- t(invG) %*% sigma %*% invG
  }
  
  return(list("beta" = beta, "cov.const" = cov.const, "eN" = eN))
}

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  mat <- crossprod(x)
  
  invMatrix <- tryCatch({
    chol2inv(chol(mat))
  },
  error = function(e) {
    if (grepl("leading minor of order", e$message)) {
      # If error is due to non-invertibility, issue warning and use generalized inverse
      warning("Calucating inverse of (t(X)%*%X), matrix is not positive-definite. Using generalized inverse.")
      return(ginv(mat))
    } else {
      # If it's another error, just stop and propagate the error
      stop(e)
    }
  })
  
  return(invMatrix)
}

W.fun = function(u,kernel){
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"
  
  if (kernel.type=="Epanechnikov") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel.type=="Uniform")      w =          0.5*(abs(u)<=1)
  if (kernel.type=="Triangular")   w =   (1-abs(u))*(abs(u)<=1)
  if (kernel.type=="Gaussian")     w =   dnorm(u)
  return(w)
}

rdbw2d_rot <- function(x,kernel.type, M){
  
  mu2K.squared <- NA
  l2K.squared <- NA
  
  if (kernel.type == "Epanechnikov"){mu2K.squared <- 1/6; l2K.squared <- 4/(3 * pi)}
  if (kernel.type == "Triangular"){mu2K.squared <- 3/20; l2K.squared <- 3/(2 * pi)}
  if (kernel.type == "Uniform"){mu2K.squared <- 1/4; l2K.squared <- 1/(pi)}
  if (kernel.type == "Gaussian"){mu2K.squared <- 1; l2K.squared <- 1/(4 * pi)}
  
  # Data cleaning
  
  x <- x[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
  x <- x[na.ok,]
  N <- dim(x)[1]
  
  # Estimate sample variance.
  
  cov.matrix <- cov(x[,c("x.1", "x.2")])
  
  D <- 2
  
  trace.const <- 1/( 2^(D+2) * pi^(D/2) * det(sqrtm(cov.matrix)) ) * ( 2*sum(diag(ginv(cov.matrix, 1e-20) %*% ginv(cov.matrix, 1e-20)))
                                                                       + (sum(diag(ginv(cov.matrix, 1e-20))))^2 )
  
  if (is.null(M)){
    hROT <-( (D * l2K.squared) / (N * mu2K.squared * trace.const) )^(1/(4+D))
  } else{
    hROT <-( (D * l2K.squared) / (M * mu2K.squared * trace.const) )^(1/(4+D))   # Adjust for mass points.
  }
  
  return(hROT)
}

rdbw2d_bw_v2 <- function(dat.centered, p, vec, dn, bn.1, bn.2 = NULL, vce, kernel, kernel_type, C){
  
  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "dist")]
  
  # Variance and coefficients for a linear combination of (p+1)-th derivatives.
  if (kernel_type == "prod"){
    w.v <- W.fun(dat.centered$x.1/c(dn), kernel) * W.fun(dat.centered$x.2/c(dn), kernel) / c(dn^2)
  } else{
    w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)
  }
  
  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)
  
  ew.v <- w.v[ind.v]
  eY.v <- dat.centered$y[ind.v]
  
  eC.v <- C[ind.v]
  
  eu.v <- dat.centered[ind.v, c("x.1", "x.2")]
  eu.v$x.1 <- eu.v$x.1/dn
  eu.v$x.2 <- eu.v$x.2/dn
  
  if (is.null(bn.2)){
    eR.v.aug <- as.matrix(get_basis(eu.v,p+1))
    eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
    eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]
  } else {
    eR.v.aug <- as.matrix(get_basis(eu.v,p+2))
    eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
    eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]
    eT.v <- eR.v.aug[, (factorial(p+1+2)/(factorial(p+1)*2) + 1) : (factorial(p+2+2)/(factorial(p+2)*2))]
  }
  
  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v
  
  w_R.v <- ew.v * eR.v
  
  invG.v <- qrXXinv(sqrtw_R.v)
  
  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)
  
  if (!is.null(bn.2)){
    sqrtw_eT.v <- sqrt(ew.v) * eT.v
    vec.t <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eT.v
    vec.t <- vec.t[1,]
    vec.t <- c(rep(0, factorial(p + 1 + 2)/(factorial(p + 1)*2)), vec.t)
  }
  
  invH.p <- get_invH(dn,p)
  
  H.p <- get_H(dn,p)
  
  beta.v <- invH.p %*% invG.v %*% t(sqrtw_R.v) %*% matrix(sqrtw_Y.v, ncol = 1)
  
  resd.v <- (eY.v - eR.v %*% (H.p %*% beta.v))[,1]
  
  lambda <- function(x){infl(x, invG.v)}
  
  if (vce=="hc0") {
    w.vce = 1
  } else if (vce=="hc1") {
    w.vce = sqrt(eN.v/(eN.v-factorial(p+2)/(factorial(p) * 2)))
  } else if (vce=="hc2") {
    hii <- apply(sqrtw_R.v, 1, lambda)
    w.vce = sqrt(1/(1-hii))
  } else if (vce=="hc3"){
    hii <- apply(sqrtw_R.v, 1, lambda)
    w.vce = 1/(1-hii)
  }
  
  resd.v <- resd.v * w.vce
  
  # sigma.v <-  t(resd.v * sqrtw_R.v) %*% (ew.v * resd.v * sqrtw_R.v) * dn^2
  
  sigma.v <- rd2d_vce(w_R.v, resd.v, eC.v, dn)
  
  V.V <- t(as.matrix(vec)) %*% t(invG.v) %*% sigma.v %*% invG.v %*% as.matrix(vec)
  V.V <- V.V[1,1]
  
  # Bias
  
  fit.ppls1 <- rd2d_lm(dat.centered, bn.1, p + 1, vce, kernel = kernel,
                       kernel_type = kernel_type, C = C, varr = TRUE)
  
  deriv.ppls1 <- fit.ppls1$beta
  B.B <- vec.q %*% deriv.ppls1
  B.B <- B.B[1,1]
  V.B <- matrix(vec.q, nrow = 1) %*% fit.ppls1$cov.const %*% matrix(vec.q, ncol = 1) / (bn.1^(2 + 2 * (p+1)))
  V.B <- V.B[1,1]
  
  Reg.v <- V.B
  
  Reg.b <- NA
  
  if (!is.null(bn.2)){
    fit.ppls2 <- rd2d_lm(dat.centered, bn.2, p + 2, vce, kernel = kernel,
                         kernel_type = kernel_type, C = C, varr = FALSE)
    deriv.ppls2 <- fit.ppls2$beta
    Reg.b <- dn * vec.t %*% deriv.ppls2
  }
  
  return(list("B" = B.B, "V" = V.V, "Reg.2" = Reg.b, "Reg.1" = Reg.v))
}

rd2d_fit_v2 <- function(dat, eval, deriv = NULL,o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                        kernel = "epa", kernel_type = "prod", vce = "hc1",
                        masspoints = "adjust", C = NULL, bwcheck = 50 + p + 1, unique = NULL){
  
  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  
  sd.x1 <- sd(dat$x.1)
  sd.x2 <- sd(dat$x.2)
  
  neval <- dim(eval)[1]
  
  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.data.frame(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.data.frame(hgrid.1)
  }
  
  if (ncol(hgrid.0) == 1){
    results <- data.frame(matrix(NA, ncol = 10, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0', 'h.1', 'mu.0', 'mu.1', 'se.0',
                           'se.1', 'eN.0', 'eN.1')
  } else {
    results <- data.frame(matrix(NA, ncol = 12, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0.x', 'h.0.y', 'h.1.x', 'h.1.y',
                           'mu.0', 'mu.1', 'se.0', 'se.1', 'eN.0', 'eN.1')
  }
  
  # Check for compatibility
  
  for (i in 1:neval){
    
    ev <- eval[i,]
    vec <- deriv[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]
    
    # Center data
    
    dat.centered <- dat[,c("x.1", "x.2", "y", "d")]
    dat.centered$x.1 <- dat.centered$x.1 - ev$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev$x.2
    # for product kernel, standardize the covariates so that they have the same sd
    dat.centered$dist <- pmax(abs(dat.centered$x.1/sd.x1), abs(dat.centered$x.2/sd.x2)) # infinity norm
    if (kernel_type == "rad") dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2) # Euclidean norm
    
    # Bandwidth restriction
    
    if (masspoints == "adjust"){
      unique.centered <- unique
      unique.centered$x.1 <- unique.centered$x.1 - ev$x.1
      unique.centered$x.2 <- unique.centered$x.2 - ev$x.2
      unique.centered$dist <- pmax(abs(unique.centered$x.1/sd.x1), abs(unique.centered$x.2/sd.x2)) # infinity norm
      if (kernel_type == "rad") unique.centered$dist <- sqrt(unique.centered$x.1^2 + unique.centered$x.2^2) # Euclidean norm
    }
    
    if (!is.null(bwcheck)){
      
      if (masspoints == "adjust"){
        sorted.0 <- sort(unique.centered[unique.centered$d == FALSE,]$dist)
        sorted.1 <- sort(unique.centered[unique.centered$d == TRUE,]$dist)
      } else{
        sorted.0   <- sort(dat.centered[dat.centered$d == FALSE,]$dist)
        sorted.1   <- sort(dat.centered[dat.centered$d == TRUE,]$dist)
      }
      
      bw.min.0   <- sorted.0[bwcheck]
      bw.min.1   <- sorted.1[bwcheck]
      bw.max.0   <- sorted.0[length(sorted.0)]
      bw.max.1   <- sorted.1[length(sorted.1)]
      
      # convert to the original if using product kernel
      if (kernel_type == "prod"){
        multiplier <- c(sd.x1, sd.x2)
      } else {
        multiplier <- c(1,1)
      }
      
      bw.min.0 <- bw.min.0 * multiplier
      bw.min.1 <- bw.min.1 * multiplier
      bw.max.0 <- bw.max.0 * multiplier
      bw.max.1 <- bw.max.1 * multiplier
      
      if (!is.null(hgrid.1)){
        h.0     <- pmax(h.0, bw.min.0)
        h.1     <- pmax(h.1, bw.min.1)
        h.0     <- pmin(h.0, bw.max.0)
        h.1     <- pmin(h.1, bw.max.1)
      } else{
        h.0 <- pmax(h.0, bw.min.0,bw.min.1)
        h.0 <- pmin(h.0, pmax(bw.max.0, bw.max.1))
        h.1 <- h.0
      }
    }
    
    fit.0.p <- rd2d_lm(dat.centered[dat.centered$d == 0,], h.0, p, vce, kernel = kernel, C = C[dat.centered$d == 0],
                       varr = TRUE, kernel_type = kernel_type)
    fit.1.p <- rd2d_lm(dat.centered[dat.centered$d == 1,], h.1, p, vce, kernel = kernel, C = C[dat.centered$d == 1],
                       varr = TRUE, kernel_type = kernel_type)
    mu.0 <- (vec %*% fit.0.p$beta)[1,1]
    mu.1 <- (vec %*% fit.1.p$beta)[1,1]
    
    # standardize
    if (length(h.0) == 1){
      h.0.x <- as.numeric(h.0)
      h.0.y <- as.numeric(h.0)
    }  else {
      h.0.x <- as.numeric(h.0[1])
      h.0.y <- as.numeric(h.0[2])
    }
    if (length(h.1) == 1){
      h.1.x <- as.numeric(h.1)
      h.1.y <- as.numeric(h.1)
    }  else {
      h.1.x <- as.numeric(h.1[1])
      h.1.y <- as.numeric(h.1[2])
    }
    
    # standard deviation
    invH.0 <- get_invH(c(h.0.x, h.0.y),p)
    se.0 <- matrix(vec, nrow = 1) %*% invH.0 %*% fit.0.p$cov.const %*% invH.0 %*% matrix(vec, ncol = 1) / (h.0.x * h.0.y)
    # se.0 <- matrix(vec, nrow = 1) %*% fit.0.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.0^(2 + 2 * o))
    se.0 <- sqrt(se.0[1,1])
    invH.1 <- get_invH(c(h.1.x, h.1.y),p)
    se.1 <- matrix(vec, nrow = 1) %*% invH.1 %*% fit.1.p$cov.const %*% invH.1 %*% matrix(vec, ncol = 1) / (h.1.x * h.1.y)
    # se.1 <- matrix(vec, nrow = 1) %*% fit.1.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.1^(2 + 2 * o))
    se.1 <- sqrt(se.1[1,1])
    
    # effective sample size
    eN.0 <- fit.0.p$eN
    eN.1 <- fit.1.p$eN
    
    if (ncol(hgrid.0) == 1){
      results[i,] <- c(ev[1], ev[2], h.0, h.1, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    } else {
      results[i,] <- c(ev[1], ev[2], h.0.x, h.0.y, h.1.x, h.1.y, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    }
  }
  
  return(results)
}

rdbw2d_cov <- function(dat, eval, deriv = NULL, o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                       kernel = "epa", kernel_type = "prod", vce = "hc2", C = NULL){
  
  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  
  neval <- dim(eval)[1]
  covs <- matrix(NA, nrow = neval, ncol = neval)
  halves.0 <- list()
  halves.1 <- list()
  inds.0 <- list()
  inds.1 <- list()
  
  clusters <- unique(C)
  g <- length(clusters)
  
  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }
  
  for (i in 1:neval){
    
    ev.a <- eval[i,]
    h.a.0 <- hgrid.0[i,]
    h.a.1 <- h.a.0
    if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
    deriv.vec.a <- deriv[i,]
    
    dat.centered <- dat[, c("x.1", "x.2", "y", "d")]
    
    dat.centered$x.1 <- dat.centered$x.1 - ev.a$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev.a$x.2
    dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)
    
    C.0 <- C[as.logical(dat.centered$d == 0)]
    C.1 <- C[as.logical(dat.centered$d == 1)]
    
    cov.half.consts.0 <- get_cov_half_v2(dat.centered[dat.centered$d == 0,], h.a.0, p, vce, kernel, kernel_type, C.0, clusters)
    cov.half.consts.1 <- get_cov_half_v2(dat.centered[dat.centered$d == 1,], h.a.1, p, vce, kernel, kernel_type, C.1, clusters)
    half.const.0 <- cov.half.consts.0$cov.half.const
    half.const.1 <- cov.half.consts.1$cov.half.const
    ind.0 <- cov.half.consts.0$ind
    ind.1 <- cov.half.consts.1$ind
    
    inds.0[[i]] <- ind.0
    inds.1[[i]] <- ind.1
    halves.0[[i]] <- half.const.0
    halves.1[[i]] <- half.const.1
    
  }
  
  for (i in 1:neval){
    for (j in i:neval){
      
      h.a.0 <- hgrid.0[i,]
      h.a.1 <- h.a.0
      if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
      h.b.0 <- hgrid.0[j,]
      h.b.1 <- h.b.0
      if (!is.null(hgrid.1)) h.b.1 <- hgrid.1[j,]
      
      if (length(h.a.0) == 1){
        h.a.0.x <- h.a.0
        h.a.0.y <- h.a.0
        h.a.1.x <- h.a.1
        h.a.1.y <- h.a.1
        h.b.0.x <- h.b.0
        h.b.0.y <- h.b.0
        h.b.1.x <- h.b.1
        h.b.1.y <- h.b.1
      } else{
        h.a.0.x <- h.a.0[1]
        h.a.0.y <- h.a.0[2]
        h.a.1.x <- h.a.1[1]
        h.a.1.y <- h.a.1[2]
        h.b.0.x <- h.b.0[1]
        h.b.0.y <- h.b.0[2]
        h.b.1.x <- h.b.1[1]
        h.b.1.y <- h.b.1[2]
      }
      deriv.vec.a <- deriv[i,]
      deriv.vec.b <- deriv[j,]
      
      half.0.a <- halves.0[[i]]
      half.0.b <- halves.0[[j]]
      half.1.a <- halves.1[[i]]
      half.1.b <- halves.1[[j]]
      
      invH.a.0 <- get_invH(c(h.a.0.x,h.a.0.y),p)
      invH.a.1 <- get_invH(c(h.a.1.x,h.a.1.y),p)
      invH.b.0 <- get_invH(c(h.b.0.x,h.b.0.y),p)
      invH.b.1 <- get_invH(c(h.b.1.x,h.b.1.y),p)
      
      if (is.null(C)){
        
        ind.0.a <- inds.0[[i]]
        ind.0.b <- inds.0[[j]]
        ind.1.a <- inds.1[[i]]
        ind.1.b <- inds.1[[j]]
        
        ind.0 <- ind.0.a * ind.0.b
        ind.0.a <- as.logical(ind.0[ind.0.a])
        ind.0.b <- as.logical(ind.0[ind.0.b])
        ind.1 <- ind.1.a * ind.1.b
        ind.1.a <- as.logical(ind.1[ind.1.a])
        ind.1.b <- as.logical(ind.1[ind.1.b])
        
        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% t(half.0.a[ind.0.a,,drop = "FALSE"]) %*% half.0.b[ind.0.b,,drop = "FALSE"] %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] / (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        # cov.0 <- cov.0[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% t(half.1.a[ind.1.a,,drop = "FALSE"]) %*% half.1.b[ind.1.b,,drop = "FALSE"] %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] / (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
        # cov.1 <- cov.1[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
      }
      
      if (!is.null(C)){
        k <- dim(deriv)[2]
        M.0 <- matrix(0, k, k)
        M.1 <- matrix(0, k, k)
        for (l in 1:g){
          M.0 <- M.0 + matrix(half.0.a[l,], ncol = 1) %*% matrix(half.0.b[l,], nrow = 1)
          M.1 <- M.1 + matrix(half.1.a[l,], ncol = 1) %*% matrix(half.1.b[l,], nrow = 1)
        }
        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% M.0 %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] /  (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% M.1 %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] /  (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
      }
      
      covs[i,j] <- cov.0 + cov.1
      covs[j,i] <- cov.0 + cov.1
    }
  }
  return(covs)
}

infl <- function(x, invG){
  result <- t(matrix(x, ncol = 1)) %*% invG %*% matrix(x, ncol = 1)
  return(result[1,1])
}

get_cov_half_v2 <- function(dat, h, p, vce = "hc1", kernel = "epa", kernel_type = "prod",
                            C = NULL, clusters = NULL){
  
  dat <- dat[,c("x.1", "x.2", "y", "d", "dist")]
  
  h <- as.vector(as.matrix(h)) # if h is data frame, convert it to a vector
  
  # checks
  if (kernel_type == "prod"){ # product kernel
    if (length(h) == 1){
      h <- c(h,h)
    }
  }
  else{ # radius kernel
    if (length(h) == 2){
      h <- sqrt(h[1]^2 + h[2]^2)
    }
  }
  
  # weights
  if (kernel_type == "prod"){
    w <- W.fun(dat$x.1/c(h[1]), kernel) * W.fun(dat$x.2/c(h[2]), kernel) / c(h[1] * h[2])
  }
  else{
    w <- W.fun(dat$dist/c(h), kernel)/c(h^2)
  }
  
  if (length(h) == 1){
    h.x <- h
    h.y <- h
  } else {
    h.x <- h[1]
    h.y <- h[2]
  }
  
  # Variance and coefficients for a linear combination of (p+1)-th derivatives
  
  ind <- as.logical(w > 0)
  
  eN <- sum(ind)
  
  ew <- w[ind]
  eY <- dat$y[ind]
  eC <- C[ind]
  
  eu <- dat[ind, c("x.1", "x.2")]
  eu$x.1 <- eu$x.1/h.x
  eu$x.2 <- eu$x.2/h.y
  
  eR <- as.matrix(get_basis(eu,p))
  
  sqrtw_R <- sqrt(ew) * eR
  sqrtw_Y <- sqrt(ew) * eY
  
  w_R <- ew * eR
  
  invG <- qrXXinv(sqrtw_R)
  
  invH.p <- get_invH(c(h.x,h.y),p)
  
  H.p <- get_H(c(h.x,h.y),p)
  
  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)
  
  resd <- abs(eY - eR %*% (H.p %*% beta))[,1]
  
  lambda <- function(x){infl(x, invG)}
  
  if (vce=="hc0") {
    w.vce <- 1
  } else if (vce=="hc1") { # TODO:Set to default
    w.vce <- sqrt(eN/(eN-factorial(p+2)/(factorial(p) * 2)))
  } else if (vce=="hc2") {
    hii <- apply(sqrtw_R, 1, lambda)
    w.vce <- sqrt(1/(1-hii))
  } else if (vce == "hc3"){
    hii <- apply(sqrtw_R, 1, lambda)
    w.vce <- 1/(1-hii)
  }
  
  resd <- resd * w.vce
  
  if (is.null(C)){
    
    sigma.half.const <-  (sqrt(ew) * resd * as.matrix(sqrtw_R)) * sqrt(h.x * h.y)
    
    cov.half.const <- sigma.half.const %*% invG
  }
  
  if (!is.null(C)){
    
    n <- length(eC)
    k <- dim(w_R)[2]
    
    g     <- length(clusters)
    w.w   <- ((n-1)/(n-k))*(g/(g-1))
    
    cov.half.const <- matrix(0,nrow = g, ncol = k)
    
    for (i in 1:g) {
      
      ind.vce <- as.logical(eC==clusters[i])
      w_R_i <- w_R[ind.vce,,drop=FALSE]
      resd_i <- resd[ind.vce]
      resd_i <- matrix(resd_i, ncol = 1)
      w_R_resd_i <- t(crossprod(w_R_i,resd_i)) * sqrt(h.x * h.y)
      
      cov.half.const[i,] <- as.vector(w_R_resd_i)
    }
    
    cov.half.const <- cov.half.const %*% invG
  }
  
  return(list("ind" = ind, "cov.half.const" = cov.half.const))
}

rd2d_cval <- function(cov, rep, side="two", alpha, lp=Inf) {
  tvec <- c()
  
  cval <- NA
  m <- dim(cov)[1]
  cov.t <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      cov.t[i,j] <- cov[i,j]/sqrt(cov[i,i] * cov[j,j])
    }
  }
  
  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)
  
  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- apply(sim, c(1), function(x){max(abs(x))})}
      else                 {tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})}
    } else if (side == "left") {
      tvec <- apply(sim, c(1), max)
    } else if (side == "right") {
      tvec <- apply(sim, c(1), min)
    }
  }
  
  if (!is.null(side)) {
    cval <- quantile(tvec, alpha/100, na.rm=T, names = F, type=2)
  }
  return(cval)
}

rd2d_pval <- function(tstat, cov, rep, side="two", lp=Inf) {
  tvec <- c()
  
  pval <- NA
  m <- dim(cov)[1]
  cov.t <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      cov.t[i,j] <- cov[i,j]/sqrt(cov[i,i] * cov[j,j])
    }
  }
  
  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)
  
  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- apply(sim, c(1), function(x){max(abs(x))})}
      else                 {tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})}
    } else if (side == "left") {
      tvec <- apply(sim, c(1), max)
    } else if (side == "right") {
      tvec <- apply(sim, c(1), min)
    }
  }
  
  if (!is.null(side)) {
    pval <- mean(tvec >= abs(tstat))
  }
  return(pval)
}

rd2d_cb <- function(mu.hat, cov.us, rep, side, alpha){
  
  # mu.hat: estimated (derivatives) of treatment effect
  # cov.us: estimated covariance matrix for treatment effects at all evaluation points
  # rep: number of repetitions for Gaussian simulation
  # side: "pos", "neg", or "two"
  # alpha: confidence level
  # If side == "two", returns upper and lower bounds of CI and CB
  # If side == "pos", returns upper bounds of CI and CB
  # If side == "neg", returns lower bounds of CI and CB
  
  se.hat <- sqrt(diag(cov.us))
  cval <- rd2d_cval(cov.us, rep = rep, side=side, alpha = alpha, lp=Inf)
  
  if (side == "two"){
    zval <- qnorm((alpha + 100)/ 200)
    CI.l <- mu.hat - zval * se.hat; CI.r <- mu.hat + zval * se.hat
    CB.l <- mu.hat - cval * se.hat; CB.r <- mu.hat + cval * se.hat
  }
  
  if (side == "left"){
    zval <- qnorm(alpha / 100)
    CI.r <- mu.hat + zval * se.hat; CI.l <- rep(-Inf, length(CI.r))
    CB.r <- mu.hat + cval * se.hat; CB.l <- rep(-Inf, length(CB.r))
  }
  if (side == "right"){
    zval <- qnorm(alpha / 100)
    CI.l <- mu.hat - zval * se.hat; CI.r <- rep(Inf, length(CI.l))
    CB.l <- mu.hat - cval * se.hat; CB.r <- rep(Inf, length(CB.l))
  }
  
  return(list(CI.l = CI.l, CI.r = CI.r, CB.l = CB.l, CB.r = CB.r))
}

get_basis <- function(u,p){
  u.x.1 <- u[,1]
  u.x.2 <- u[,2]
  result <- matrix(NA, nrow = dim(u)[1], ncol = factorial(p+2)/(factorial(p) * 2))
  result[,1] <- rep(1, dim(u)[1])
  count <- 2
  if (p >= 1){
    for (j in 1:p){
      for (k in 0:j){
        result[,count] <- u.x.1^(j-k) * u.x.2^k
        count <- count + 1
      }
    }
  }
  return(result)
}

get_H <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1
  
  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- h.x.1^(j-k) * h.x.2^k
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

get_invH <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1
  
  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- 1/(h.x.1^(j-k) * h.x.2^k)
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

rd2d_lm <- function(dat, h, p, vce = "hc1", kernel = "epa", kernel_type = "prod",
                    C = NULL, varr = FALSE){
  
  dat <- dat[,c("x.1", "x.2", "y", "d", "dist")]
  
  # Variance and coefficients for a linear combination of (p+1)-th derivatives.
  
  h <- as.vector(as.matrix(h)) # if h is data frame, convert it to a vector
  
  # checks
  if (kernel_type == "prod"){ # product kernel
    if (length(h) == 1){
      h <- c(h,h)
    }
  }
  else{ # radius kernel
    if (length(h) == 2){
      h <- sqrt(h[1]^2 + h[2]^2)
    }
  }
  
  # weights
  if (kernel_type == "prod"){
    w <- W.fun(dat$x.1/c(h[1]), kernel) * W.fun(dat$x.2/c(h[2]), kernel) / c(h[1] * h[2])
  }
  else{
    w <- W.fun(dat$dist/c(h), kernel)/c(h^2)
  }
  
  if (length(h) == 1){
    h.x <- h; h.y <- h
  } else {
    h.x <- h[1]; h.y <- h[2]
  }
  
  ind <- as.logical(w > 0)
  
  eN <- sum(ind)
  
  ew <- w[ind]
  eY <- dat$y[ind]
  eC <- C[ind] # if C == NULL, eC == NULL.
  
  eu <- dat[ind, c("x.1", "x.2")]
  eu$x.1 <- eu$x.1/h.x
  eu$x.2 <- eu$x.2/h.y
  
  eR <- as.matrix(get_basis(eu,p))
  
  sqrtw_R <- sqrt(ew) * eR
  sqrtw_Y <- sqrt(ew) * eY
  
  w_R <- ew * eR
  
  invG <- qrXXinv(sqrtw_R)
  
  invH.p <- get_invH(c(h.x,h.y),p)
  
  H.p <- get_H(c(h.x,h.y),p)
  
  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)
  
  cov.const <- NA
  
  if (varr){
    
    resd <- (eY - (eR %*% H.p) %*% beta)[,1]
    
    lambda <- function(x){infl(x, invG)}
    
    if (vce=="hc0") {
      w.vce = 1
    } else if (vce=="hc1") {
      w.vce = sqrt(eN/(eN-factorial(p+2)/(factorial(p) * 2)))
    } else if (vce=="hc2") {
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = 1/(1-hii)
    }
    
    resd <- resd * w.vce
    
    sigma <- rd2d_vce(w_R, resd, eC, c(h.x, h.y))
    
    # sigma <-  t(resd * as.matrix(sqrtw_R)) %*% (ew * resd * as.matrix(sqrtw_R)) * h^2
    
    cov.const <- t(invG) %*% sigma %*% invG
  }
  
  return(list("beta" = beta, "cov.const" = cov.const, "eN" = eN))
}

get_coeff <- function(dat.centered,vec, p,dn, kernel, kernel_type){
  
  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "dist")]
  if (kernel_type == "prod"){
    w.v <- W.fun(dat.centered$x.1/c(dn), kernel) * W.fun(dat.centered$x.2/c(dn), kernel) / c(dn * dn)
  }
  else{
    w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)
  }
  
  # w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)
  
  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)
  
  ew.v <- w.v[ind.v]
  eY.v <- dat.centered$y[ind.v]
  
  eu.v <- dat.centered[ind.v, c("x.1", "x.2")]
  eu.v$x.1 <- eu.v$x.1/dn
  eu.v$x.2 <- eu.v$x.2/dn
  
  eR.v.aug <- as.matrix(get_basis(eu.v,p+1))
  eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
  eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]
  
  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v
  
  invG.v <- qrXXinv(sqrtw_R.v)
  
  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)
  
  return(vec.q)
}


# new version with one h for each coordinate
rd2d_vce <- function(w_R, resd, eC, h){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  
  n <- length(eC)
  k <- dim(w_R)[2]
  M <- matrix(0, nrow = k, ncol = k)
  if (is.null(eC)){
    w.w <- 1
    M <-  crossprod(resd * as.matrix(w_R)) * h.x.1 * h.x.2
  }
  else{
    clusters = unique(eC)
    g     = length(clusters)
    w.w =((n-1)/(n-k))*(g/(g-1))
    for (i in 1:g) {
      ind=eC==clusters[i]
      # w_R_i = w_R[ind,,drop=FALSE]
      # Attempt to subset w_R with ind
      w_R_i <- tryCatch(
        {
          w_R[ind, , drop = FALSE]
        },
        error = function(e) {
          cat("Error: ", conditionMessage(e), "\n")
          cat("Dimensions of w_R: ", paste(dim(w_R), collapse = " x "), "\n")
          cat("Length of ind: ", length(ind), "\n")
          stop("Exiting due to the above error.")
        }
      )
      resd_i = resd[ind]
      resd_i <- matrix(resd_i, ncol = 1)
      w_R_resd_i = t(crossprod(w_R_i,resd_i))
      M = M + crossprod(w_R_resd_i,w_R_resd_i) * h.x.1 * h.x.2
    }
  }
  return(M * w.w)
}

rd2d_unique <- function(dat){
  
  dat <- dat[,c("x.1", "x.2", "y", "d")]
  
  ord <- order(dat[,1], dat[,2])
  
  dat <- dat[ord,]
  
  N <- dim(dat)[1]
  
  # if x has one or no element
  if (N == 0) return(list(unique = NULL, freq = c(), index = c()))
  if (N == 1) return(list(unique = dat, freq = 1, index = 1))
  
  # else
  uniqueIndex <- c(c(dat[2:N, 1] != dat[1:(N-1),1] | dat[2:N, 2] != dat[1:(N-1),2]), TRUE)
  unique <- dat[uniqueIndex,]
  nUnique <- dim(unique)[1]
  
  # all are distinct
  if (nUnique == N) return(list(unique=unique, freq=rep(1,N), index=1:N))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=N, index=N))
  
  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1
  
  return(list(unique=unique, freq=freq, index=(1:N)[uniqueIndex]))
}

# ######################### Estimation #############################
# 

panel.fcov.df <- panel.y.df %>% select(x,y,year,treat,fcover)

years <- c(2000:2001)
est <- list()
i <- 1
for (yr in years) {
panel.fcov.y.df <- panel.fcov.df[panel.fcov.df$year==yr,]
Y = panel.fcov.y.df$fcover
X = panel.fcov.y.df[,c("x","y")]
t= panel.fcov.y.df$treat
b = bpoints.df[1:10,2:3]
est[[i]] <- rd2d(Y = Y,
     X = X,
     t = t,
     b = b)


i <- i + 1

}
