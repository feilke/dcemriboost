blfitspat <- function(u,nu,m,time,allkep,factor,whichnb,lambda,kep_penal){  
  
  a1 <- 3.99
  a2 <- 4.78
  m1 <- 0.144
  m2 <- 0.0111
  D <- 0.1
  
  expmod <- function(beta) {
    c.t <- a1 * (exp(-beta[1L] * time) - exp(-m1 * time)) / (m1 - beta[1L]) 
    c.t <- c.t + a2 * (exp(-beta[1L] * time) - exp(-m2 * time)) / (m2 - beta[1L]) 
    c.t * D * beta[2L]
  }
    
  lower <- c(0.05, 0.01)
  upper <- c(20, 20)
  ps <- makeNumericParamSet("beta", len = 2L, lower = lower, upper = upper)
  
  starts <- generateDesign(n = 10, par.set = ps) 
    
  kk_esti <- NA
  
  howmanynb <- length(whichnb)
    
  penalfunc <- function(par){
    penal <- 0
    for(i in 1:howmanynb){
      vec <- c()
      for(j in 1:length(kep_penal[[i]])){
        vec <- c(vec,(par[1L]-kep_penal[[i]][j])^2)
      }
      minvec <- min(vec,na.rm=TRUE)
      penal <- penal + minvec
    }
    penal
  }
  
    objL2 <- function(par) {
      penal <- penalfunc(par)
      r <- (u - (expmod(par) + lambda * penal)) 
      crossprod(r) / 2
    }

    objRes <- function(par) {
      penal <- penalfunc(par)
      u - (expmod(par) + lambda * penal) 
    }


    zlist <- apply(X=starts,MARGIN=1,FUN=function(x) nls.lm(fn = objRes, par = c(x[1],x[2]), lower = lower, upper = upper, jac = NULL, nls.control(maxiter = 100)))

  loss <- lapply(zlist,FUN=function(x) objL2(x$par)) 
  beta <- lapply(zlist,FUN=function(x) x$par)
  
  w <- which.min(loss) 
  mod <- zlist[[w]] 

  kk_esti <- c(mod$par[1],nu*mod$par[2])  
  names(kk_esti) <- c("kep","Ktrans")

  kep <- mod$par[1]
  Ktrans <- mod$par[2]
  
  b<-((a1/(m1-kep))*(exp(-(time*kep))-exp(-(time*m1)))+(a2/(m2-kep))*(exp(-(time*kep))-exp(-(time*m2))))
    
  penal <- penalfunc(mod$par)

  fbl <- Ktrans*(D*(b)) + lambda * penal

  res <- (u-nu*fbl) 
  n <- length(res) 

  sse <- (1/n)*sum(res^2)

  if(m>1){
    
  dist_ok <- 0

  if(length(unique(allkep))==1){
    
    upperkep <- unique(allkep)*factor
    lowerkep <- unique(allkep)/factor
    
    if(kk_esti[1] >= upperkep || kk_esti[1] <= lowerkep){
      dist_ok <- 1 
    }
    
    if(dist_ok==0){
      solKtrans <- getKtrans(allkep, a1, m1, a2, m2, D, u, time, nu, m)
      sse <- solKtrans$sse
      fbl <- solKtrans$fbl
      kk_esti <- solKtrans$kk_esti
      n <- solKtrans$n
    }
  }
  
  if(length(unique(allkep))>1){
    upperkep_vec <- unique(allkep[!is.na(allkep)])*factor 
    lowerkep_vec <- unique(allkep[!is.na(allkep)])/factor
    
    dist_ok1 <- vector("numeric",length=length(upperkep_vec)) 
    
    for(l in 1:length(upperkep_vec)){
      if(kk_esti[1] >= upperkep_vec[l] || kk_esti[1] <= lowerkep_vec[l]){
        dist_ok1[l] <- 1       
      }
    }  
    if(any(dist_ok1==0)){
      solKtrans <- getKtrans(allkep, a1, m1, a2, m2, D, u, time, nu, m)
      sse <- solKtrans$sse
      fbl <- solKtrans$fbl
      kk_esti <- solKtrans$kk_esti
      n <- solKtrans$n
    }
  }
  }
  return(list("sse"=sse,"fbl"=fbl,"kk_esti"=kk_esti,"n"=n)) 
}


basefitspat <- function(u,nu,coeff,m,mstop,niter,fit,y,time,SSE,DF,IC,keps,factor,SSEbefore,ALLKEP,whichnb,lambda,kep_penal){ 
  
  sse <- 0
  kkcoef <- 0   
  ss <- NULL 
  df <- 0
  ic <- 0 
  allkep <- NULL
  
  if(m>1){
    allkep <- ALLKEP
  }
      
  blfit_p <- blfitspat(u=u,nu=nu,m=m,time=time,allkep=allkep,factor=factor,whichnb=whichnb,lambda=lambda,kep_penal=kep_penal)  
    
  sse <- blfit_p$sse 
  ss  <- blfit_p$fbl 
    
  allkep <- c(allkep,blfit_p$kk_esti[1])
  
  df <- 2*length(unique(allkep))  
        
  n <- blfit_p$n
  sigma_hat_sq <-  sse
  
  ic <- -n*log(n)+n*log(n*sigma_hat_sq)+ log(n)*df 

  kkcoef <- blfit_p$kk_esti 
  
  if(m>1){
    if (1-sse/SSEbefore < 1e-08){       
      mstop <- 1 
      return(list("f"=NA,"minsse"=NA,"coef"=NA,"df"=NA,"ic"=NA,"keps"=NA))
    }
  }
  
  if(mstop==0){ 
    minsse <- sse 
    ic_sel <- ic
    df_sel <- df
    coef <- kkcoef
    keps <- kkcoef[1]
    f <- ss 
    
    return(list("f"=f,"minsse"=minsse,"coef"=coef,"df"=df_sel,"ic"=ic_sel,"keps"=keps)) 
  }
}

boostspat <- function(mstop,niter,y,fit,nu,time,factor,M,SSEbefore,ALLKEP,whichnb,lambda,kep_penal){
  
  coeff <- vector(mode = "list", length = niter)
  keps <- vector("numeric", length=niter)
  
  DF <- rep(0,times=niter)
  IC <- rep(0,times=niter)
  
  SSE <- rep(NA,times=niter)

  m <- M 
    
    u <- ngradient(resp=y, pred=fit) 

    basess <- basefitspat(u,nu,coeff,m,mstop,niter,fit,y,time,SSE,DF,IC,keps,factor,SSEbefore,ALLKEP,whichnb,lambda,kep_penal) 

    if(!is.na(basess$minsse)){        
      
      fit <- fit + nu * basess$f 
      
      if (any(!is.finite(u[!is.na(u)]))) 
        stop("Infinite residuals: please decrease step-size nu in ",sQuote("boost_control"))

      coeff <- basess$coef 
      SSE <- basess$minsse 
      IC <- basess$ic 
      DF <- basess$df 
    }

  return(list("mstop"=mstop, "coeff"=coeff, "fit"=fit, "u"=u, "IC"=IC, "DF"=DF, "SSE"=SSE, "m"=m)) 
}

#' Spatial boosting
#' 
#' @param response numeric vector containing the measurements at time points specified in argument time
#' @param offset offset
#' @param niter number of iterations
#' @param nu step length factor
#' @param time numeric vector containing the time points at which measurements are available
#' @param factor factor by which the kep-values should differ at least
#' @param M current iteration
#' @param SSEbefore residual sum of squares in the last iteration
#' @param ALLKEP kep-estimates accepted in previous iterations
#' @param whichnb index for the neighbors of the voxel
#' @param lambda penalization parameter for spatial regularization
#' @param kep_penal kep-estimates of the neighbors of the voxel
#' 
#' @return response: numeric vector containing the measurements at time points specified in argument time
#' @return coeff: kep- and Ktrans-estimates
#' @return fit: numeric vector containing the current fit
#' @return u: numeric vector containing the current residuals
#' @return DF: degrees of freedom
#' @return IC: value of the information criterion (BIC)
#' @return SSE: residual sum of squares
#' @return m: current iteration
#' 
#' @export spatboost
#' @import minpack.lm
#' @import ParamHelpers
#' 
spatboost <- function(response,offset,niter,nu,time,factor,M,SSEbefore,ALLKEP,whichnb,lambda,kep_penal){ 
  
  mstop <- 0
  
  nu <- nu
  
  niter <- niter
  
  time <- time
  
  SSEbefore <- SSEbefore
  
  y <- response     
  
  if (is.null(offset)){
    fit <-  rep(0,times=length(y)) 
  }
  if (!is.null(offset)){
    fit <- offset 
  }
  
  boost_exec <- boostspat(mstop,niter,y,fit,nu,time,factor,M,SSEbefore,ALLKEP,whichnb,lambda,kep_penal)
  
  return(list("response"=response,"coeff"=boost_exec$coeff,"fit"=boost_exec$fit,"u"=boost_exec$u,"DF"=boost_exec$DF,"IC"=boost_exec$IC,"SSE"=boost_exec$SSE,"m"=boost_exec$m))  
}