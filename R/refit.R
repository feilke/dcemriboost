#' Refit
#' 
#' @param response numeric vector containing the measurements at time points specified in argument time
#' @param response_true numeric vector containing the measurements without error at time points specified in vector time (for simulated measurements only)
#' @param time numeric vector containing the time points at which measurements are available
#' @param uniquekeps unique selected kep-estimates
#' 
#' @return len: number of unique selected kep-estimates
#' @return sse: residual sum of squares
#' @return sse_true: residual sum of squares for measurements without error (for simulated measurements only)
#' @return kk_esti: kep- and Ktrans-estimates
#' @export refit
#' @import quadprog
#'  
refit <- function(response,response_true=NULL,time,uniquekeps){
  
  y <- response 
  
  if(!is.null(response_true)){
    y_true <- response_true
  }

  keps <- uniquekeps
  len <- length(keps)
  
  D <- 0.1
  a1 <- 3.99
  a2 <- 4.78
  m1 <- 0.144
  m2 <- 0.0111
  
  model.weinmann <- function(time, kep, D, a1, a2, m1, m2)
  {
    erg <- D  * ((a1/(m1 - kep)) * (exp(-(time * 
                                          kep)) - exp(-(time * m1))) + (a2/(m2 - kep)) *
                   (exp(-(time * kep)) - exp(-(time * m2))))
    erg[time <= 0] <- 0
    return(erg)
  }
    
  x <- c()
  for (i in 1:len)
  {
    ft <- model.weinmann(time, keps[i], D, a1, a2, m1, m2)
    x <- cbind(x,ft)
  }
  
  y <- matrix(y,length(y),1) 
  
  dvec <- t(x)%*%y
  Dmat <- t(x)%*%x 
  
  Amat <- matrix(0,nrow=len,ncol=2*len)
  j <- 1
  i <- 1
  while(i<=nrow(Amat)){
    Amat[i,j] <- 1
    Amat[i,j+1] <- -1
    i <- i + 1
    j <- j + 2
  } 
  
  bvec <- matrix(rep(c(0.01,-20),times=len),2*len,1)
  
  z <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
  
  Ktrans <- z$solution 
  fbl <- vector("numeric",length=length(y))

  kk_esti <- c()
  
  for(i in 1:length(Ktrans)){
    kk_esti <- c(kk_esti,keps[i],Ktrans[i])
    fbl <- fbl + Ktrans[i]*(D*((a1/(m1-keps[i]))*(exp(-(time*keps[i]))-exp(-(time*m1)))+(a2/(m2-keps[i]))*(exp(-(time*keps[i]))-exp(-(time*m2)))))
  }

  res <- (y-fbl) 
  sse <- sum(res^2) 

  if(!is.null(response_true)){
    res_true <- (y_true-fbl) 
    sse_true <- sum(res_true^2) 
    return(list(len=len,sse=sse,sse_true=sse_true,kk_esti=kk_esti))
  }

  if(is.null(response_true)){
    return(list(len=len,sse=sse,kk_esti=kk_esti))
  }

}