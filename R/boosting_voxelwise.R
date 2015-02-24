blfit <- function(u,nu,m,time,allkep,factor){ 
  
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
  
  objL2 <- function(beta) {
    r <- (u - expmod(beta))
    crossprod(r) / 2
  }
  
  objRes <- function(beta) {
    u - expmod(beta)
  }
  
  jac <- function(beta){
    jac.beta1 <- -(beta[2L] * (D * (a1/(m1 - beta[1L])^2 * (exp(-(time * beta[1L])) - exp(-(time * m1))) - (a1/(m1 - beta[1L])) * (exp(-(time * beta[1L])) * time) + (a2/(m2 - beta[1L])^2 * (exp(-(time * beta[1L])) - exp(-(time * m2))) - (a2/(m2 - beta[1L])) * (exp(-(time * beta[1L])) * time)))))
    jac.beta2 <- -(D * ((a1/(m1 - beta[1L])) * (exp(-(time * beta[1L])) - exp(-(time * m1))) + (a2/(m2 - beta[1L])) * (exp(-(time * beta[1L])) - exp(-(time * m2)))))
    c(jac.beta1,jac.beta2)
  }
  
  lower <- c(0.05, 0.01)
  upper <- c(20, 20)
  ps <- makeNumericParamSet("beta", len = 2L, lower = lower, upper = upper)
  
  starts <- generateDesign(n = 10, par.set = ps) 
  
  kk_esti <- NA

  zlist <- apply(X=starts,MARGIN=1,FUN=function(x) nls.lm(fn = objRes, par = c(x[1],x[2]), lower = lower, upper = upper, jac = jac, nls.control(maxiter = 100)))
  
  loss <- lapply(zlist,FUN=function(x) objL2(x$par)) 
  beta <- lapply(zlist,FUN=function(x) x$par)

  w <- which.min(loss) 
  mod <- zlist[[w]] 

  kk_esti <- c(mod$par[1],nu*mod$par[2])
  names(kk_esti) <- c("kep","Ktrans")
  
  kep <- mod$par[1]
  Ktrans <- mod$par[2]
  
  fbl <- Ktrans*(D*((a1/(m1-kep))*(exp(-(time*kep))-exp(-(time*m1)))+(a2/(m2-kep))*(exp(-(time*kep))-exp(-(time*m2)))))

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
    upperkep_vec <- unique(allkep)*factor 
    lowerkep_vec <- unique(allkep)/factor
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

basefit <- function(u,nu,coeff,m,mstop,niter,fit,y,time,SSE,DF,IC,keps,factor){ 
  
  sse <- 0
  kkcoef <- 0   
  ss <- NULL 
  
  df <- 0
  ic <- 0
  
  allkep <- NULL
  
  if(m>1){

  un <- unlist(coeff[1:(m-1)])

  allkep <- un[names(un)=="kep"]
  }
      
  blfit_p <- blfit(u=u,nu=nu,m=m,time=time,allkep=allkep,factor=factor)  
    
  sse <- blfit_p$sse 
  ss  <- blfit_p$fbl 
  
  allkep <- c(allkep,blfit_p$kk_esti[1])
  dfkep <- 2*length(unique(allkep))  

  df <- dfkep
        
  n <- blfit_p$n
  sigma_hat_sq <-  sse
  
  ic <- -n*log(n)+n*log(n*sigma_hat_sq)+ log(n)*df
    
  kkcoef <- blfit_p$kk_esti
  
  if(m>1){ 
    if (1-sse/SSE[m-1] < 1e-08){ 
      mstop <- mstop + niter
      return(list("f"=NA,"minsse"=NA,"coef"=NA,"df"=NA,"ic"=NA,"keps"=NA))
    }
  }

  if(m > mstop){  
    minsse <- sse 
    ic_sel <- ic
    df_sel <- df
    coef <- kkcoef
    keps <- kkcoef[1]
    f <- ss 
    
    return(list("f"=f,"minsse"=minsse,"coef"=coef,"df"=df_sel,"ic"=ic_sel,"keps"=keps)) 
  }
}

boost <- function(mstop,niter,y,fit,nu,time,factor){
  
  coeff <- vector(mode = "list", length = niter)
  keps <- vector("numeric", length=niter)
  
  DF <- rep(0,times=niter)
  IC <- rep(0,times=niter)
  
  SSE <- rep(NA,times=niter)
  
  m <- (mstop + 1) 
  while(m <=(mstop + niter)){
    
    u <- ngradient(resp=y, pred=fit) 

    basess <- basefit(u,nu,coeff,m,mstop,niter,fit,y,time,SSE,DF,IC,keps,factor) 
    
    if(is.na(basess$minsse)){ 
      m <- mstop + niter + 1
    }
    
    if(!is.na(basess$minsse)){        
      
      fit <- fit + nu * basess$f 
      
      if (any(!is.finite(u[!is.na(u)]))) 
        stop("Infinite residuals: please decrease step-size nu in ",sQuote("boost_control"))
            
      coeff[[m]] <- basess$coef 
      SSE[m] <- basess$minsse
      IC[m] <- basess$ic
      DF[m] <- basess$df
      
      m <- m + 1
    }
  }
  mstop <- mstop + niter
  return(list("mstop"=mstop, "coeff"=coeff, "fit"=fit, "u"=u, "IC"=IC, "DF"=DF, "SSE"=SSE)) 
}

#' Voxelwise boosting
#' 
#' @param response numeric vector containing the measurements at time points specified in argument time
#' @param offset offset
#' @param niter number of iterations
#' @param nu step length factor
#' @param time numeric vector containing the time points at which measurements are available
#' @param factor factor by which the kep-values should differ at least
#' 
#' @return nkeps: number of selected kep-estimates
#' @return keps_un: unique selected kep-estimates
#' @return Ktrans: selected Ktrans-estimates
#' @return keps: kep-estimates

#' @export boostvw
#' @import minpack.lm
#' @import ParamHelpers
#' 
boostvw <- function(response,offset,niter,nu,time,factor){ 
  
  mstop <- 0
  
  nu <- nu
  
  niter <- niter
  
  time <- time
  
  y <- response     
  
  if (is.null(offset)){
    fit <-  rep(0,times=length(y)) 
  }
  if (!is.null(offset)){
    fit <- offset 
  }
  
  boost_exec <- boost(mstop,niter,y,fit,nu,time,factor)

  minic <- which.min(boost_exec$IC) 
  rescoeff <- boost_exec$coeff[1:minic] 
  
  keps <- NULL
  ktrans <- NULL
  for(l in 1:length(rescoeff)){
    keps <- c(keps,rescoeff[[l]][[1]])
    ktrans <- c(ktrans,rescoeff[[l]][[2]])
  }
  nkeps <- length(unique(keps))
  keps_un <- unique(keps)
  
  return(list("nkeps"=nkeps,"keps_un"=keps_un,"ktrans"=ktrans,"keps"=keps))
}