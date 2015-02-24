ngradient <- function(resp, pred){ 
  return(resp - pred)
}

getKtrans <- function(allkep, a1, m1, a2, m2, D, u, time, nu, m){
  kep <- allkep[m-1]
  
  b <- ((a1/(m1-kep))*(exp(-(time*kep))-exp(-(time*m1)))+(a2/(m2-kep))*(exp(-(time*kep))-exp(-(time*m2))))
  anasol <- u%*%b/(D*(b%*%b))
  
  lower <- 0.01
  upper <- 20
  
  if(anasol<lower){
    anasol <- lower
  }
  
  if(anasol>upper){
    anasol <- upper
  }
  
  Ktrans <- anasol  
  
  kk_esti <- c(kep,nu*Ktrans)
  names(kk_esti) <- c("kep","Ktrans")
  
  fbl <- Ktrans*(D*(b))
  
  res <- (u-nu*fbl) 
  n <- length(res) 
  
  sse <- (1/n)*sum(res^2) 
  
  return(list("sse"=sse,"fbl"=fbl,"kk_esti"=kk_esti,"n"=n))
} 