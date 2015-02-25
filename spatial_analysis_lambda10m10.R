require(BatchExperiments)
require(parallel)
options(cores=40) 
multicore=TRUE
require(dcemriboost)
data(dcemrisim)

I <- dim(CONC_SIM_NOISE)[1] 
J <- dim(CONC_SIM_NOISE)[2] 
N <- I*J 
K <- dim(CONC_SIM_NOISE)[4] 

timep <- length(time)

xypix <- cbind(sort(rep(1:I,J)),rep(1:J,I))
colnames(xypix) <- c("x","y")
xysample <- xypix

xyblack <- xysample[which((xysample[,1]+xysample[,2])%%2==0, arr.ind=TRUE),1:2]
xywhite <- xysample[which((xysample[,1]+xysample[,2])%%2==1, arr.ind=TRUE),1:2]

xyblack_list <- c()
xywhite_list <- c()
for(i in 1:length(xywhite[,1]))
{
  xyblack_list[[i]] <- xyblack[i,]
  xywhite_list[[i]] <- xywhite[i,]
}

xy_list <- c()
for(i in 1:length(xysample[,1]))
{
  xy_list[[i]] <- xysample[i,]
}

Response <- array(NA, c(I*J,T=timep,K))
Response_true <- array(NA, c(I*J,T=timep,K))

for(k in 1:K)
{
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      Response[(i-1)*J+j,,k] <- CONC_SIM_NOISE[i,j,,k]
      Response_true[(i-1)*J+j,,k] <- CONC_SIM[i,j,]
    }
  }
}

NumberComp_true <- NUMBER_COMP[1:I,1:J,1]

summed_curves_allK <- list()

for(k in 1:K)
{
  summed_curve <- array(0, timep)
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      summed_curve <- summed_curve +  CONC_SIM_NOISE[i,j,,k]
    }
  }
  summed_curves_allK[[k]] <- summed_curve
}

MeanResponse <- array(NA, c(I*J,timep,K))
for(k in 1:K)
{
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      MeanResponse[(i-1)*J+j,,k] <- summed_curves_allK[[k]]/N
    }
  }
}

keps_pen_list <- list()
set.seed(182)
seeds <- rnorm(K,0,1)
for(k in 1:K){
  set.seed(seeds[k])
  model <- boostvw(response=MeanResponse[1,,k],offset=NULL,niter=100,nu=0.1,time=time,factor=5) 
  keps <- NULL
  ktrans <- NULL
  for(l in 1:length(model$ktrans)){
    keps <- c(keps,model$keps[[l]])
    ktrans <- c(ktrans,model$ktrans[[l]])
  }
  nkeps <- length(unique(keps)) 
  keps_un <- unique(keps)
  
  keps_penal <- vector("list",length=I*J)
  keps_pen_list[[k]] <- lapply(keps_penal,FUN=function(x) x <- keps_un)
}

keps_pen <- keps_pen_list 

WHICHNB <- mclapply(xy_list, WHICH_NB, I, J, Response[,,1]) 

xywhite_pix <- as.numeric(lapply(xywhite_list, FUN= function(x) indx(x,I,J)))
xyblack_pix <- as.numeric(lapply(xyblack_list, FUN= function(x) indx(x,I,J)))

lambda <- 10^(-10)

howoften <- 100 
set.seed(42)
se <- sample(1:100000,howoften*K)
seeds <- array(se,c(howoften,K))

response_allvoxels <- array(NA, c(I*J,timep,K)) 
Ktrans_allvoxels <- array(NA, c(I*J,howoften,K))
kep_allvoxels <- array(NA, c(I*J,howoften,K))
fit_allvoxels <- array(NA, c(I*J,timep,K)) 
u_allvoxels <- array(NA, c(I*J,timep,K)) 
DF_allvoxels <- array(NA, c(I*J,howoften,K))
IC_allvoxels <- array(NA, c(I*J,howoften,K))
SSE_allvoxels <- array(NA, c(I*J,howoften,K))

parallelfit <- function(index,xysample,Response,k,I,J,offset,M,SSE_allvoxels,fit_allvoxels,kep_allvoxels,WHICHNB,lambda,keps_pen,done,seeds){
  
  y <- Response[index,,k]
  SSEbefore <- NULL
  done <- 0
  
  if(M>1){
    SSEbefore <- SSE_allvoxels[index,M-1,k]
    offset <- fit_allvoxels[index,,k]
    ALLKEP <- kep_allvoxels[index,1:(M-1),k]
    
    keps_penalis <- kep_allvoxels[,1:(M-1),k]
    if(M==2){
      keps_penalis <- as.matrix(kep_allvoxels[,1:(M-1),k])
    }
    
    keps_pen <- as.list(apply(keps_penalis,MARGIN=1,FUN=function(x) unique(x))) 
    
    if(is.na(SSE_allvoxels[index,M-1,k])){
      done <- 1
    }
  }
  
  if(M==1){ 
    keps_pen <- keps_pen[[k]]
  }
  
  if(done==1){
    return(list(response=NA,coeff=NA,fit=NA,u=NA,DF=NA,IC=NA,SSE=NA,m=NA,done=done,k=k,index=index))
  }
  
  if(done==0){
    set.seed(seeds[M,k])
      
    model <- spatboost(response=y,offset=offset,niter=1,nu=0.1,time=time,factor=5,
                           M=M,SSEbefore=SSEbefore,ALLKEP=ALLKEP,whichnb=WHICHNB[[index]],lambda=lambda,
                           kep_penal=keps_pen[WHICHNB[[index]]]) 
    
    return(list(response=model$response,coeff=model$coeff,fit=model$fit,u=model$u,
                DF=model$DF,IC=model$IC,SSE=model$SSE,m=model$m,done=done,k=k,index=index))
  }
}

for(k in 1:K){  
    
    offset <- NULL
      
    for(M in 1:howoften){
      
      print("M")
      print(M)
       
      if(M<=howoften){
        
        res_factor5_white <- mclapply(xywhite_pix,FUN=parallelfit,xysample=xysample,Response=Response,k,I=I,J=J,offset=offset,M=M,SSE_allvoxels=SSE_allvoxels,fit_allvoxels=fit_allvoxels,kep_allvoxels=kep_allvoxels,WHICHNB=WHICHNB,lambda=lambda,keps_pen=keps_pen,done=done,seeds=seeds) 
      
        for(i in 1:length(res_factor5_white)){
          done <-  res_factor5_white[[i]]$done
          if(!is.null(res_factor5_white[[i]]$coeff[[1]])){
            Ktrans_allvoxels[xywhite_pix[i],M,k] <- res_factor5_white[[i]]$coeff[2]
            kep_allvoxels[xywhite_pix[i],M,k] <- res_factor5_white[[i]]$coeff[1] 
            response_allvoxels[xywhite_pix[i],,k] <- res_factor5_white[[i]]$response
            fit_allvoxels[xywhite_pix[i],,k] <- res_factor5_white[[i]]$fit
            u_allvoxels[xywhite_pix[i],,k] <- res_factor5_white[[i]]$u
            IC_allvoxels[xywhite_pix[i],M,k] <- res_factor5_white[[i]]$IC
            DF_allvoxels[xywhite_pix[i],M,k] <- res_factor5_white[[i]]$DF
            SSE_allvoxels[xywhite_pix[i],M,k] <- res_factor5_white[[i]]$SSE
          }
          if(is.null(res_factor5_white[[i]]$coeff[[1]]) || (done==1)){
            Ktrans_allvoxels[xywhite_pix[i],M,k] <- NA
            kep_allvoxels[xywhite_pix[i],M,k] <- NA
            response_allvoxels[xywhite_pix[i],,k] <- response_allvoxels[xywhite_pix[i],,k]
            fit_allvoxels[xywhite_pix[i],,k] <- fit_allvoxels[xywhite_pix[i],,k]
            u_allvoxels[xywhite_pix[i],,k] <- u_allvoxels[xywhite_pix[i],,k]
            IC_allvoxels[xywhite_pix[i],M,k] <- NA
            DF_allvoxels[xywhite_pix[i],M,k] <- NA
            SSE_allvoxels[xywhite_pix[i],M,k] <- NA
          }          
        }
        
        res_factor5_black <- mclapply(xyblack_pix,FUN=parallelfit,xysample=xysample,Response=Response,k,I=I,J=J,offset=offset,M=M,SSE_allvoxels=SSE_allvoxels,fit_allvoxels=fit_allvoxels,kep_allvoxels=kep_allvoxels,WHICHNB=WHICHNB,lambda=lambda,keps_pen=keps_pen,done=done,seeds=seeds)
        
        for(i in 1:length(res_factor5_black)){
          
          done <-  res_factor5_black[[i]]$done
          
          if(!is.null(res_factor5_black[[i]]$coeff[[1]])){
            Ktrans_allvoxels[xyblack_pix[i],M,k] <- res_factor5_black[[i]]$coeff[2]
            kep_allvoxels[xyblack_pix[i],M,k] <- res_factor5_black[[i]]$coeff[1] 
            response_allvoxels[xyblack_pix[i],,k] <- res_factor5_black[[i]]$response
            fit_allvoxels[xyblack_pix[i],,k] <- res_factor5_black[[i]]$fit
            u_allvoxels[xyblack_pix[i],,k] <- res_factor5_black[[i]]$u
            IC_allvoxels[xyblack_pix[i],M,k] <- res_factor5_black[[i]]$IC
            DF_allvoxels[xyblack_pix[i],M,k] <- res_factor5_black[[i]]$DF
            SSE_allvoxels[xyblack_pix[i],M,k] <- res_factor5_black[[i]]$SSE
          }
          if(is.null(res_factor5_black[[i]]$coeff[[1]])  || (done==1)){
            Ktrans_allvoxels[xyblack_pix[i],M,k] <- NA
            kep_allvoxels[xyblack_pix[i],M,k] <- NA
            response_allvoxels[xyblack_pix[i],,k] <- response_allvoxels[xyblack_pix[i],,k]
            fit_allvoxels[xyblack_pix[i],,k] <- fit_allvoxels[xyblack_pix[i],,k]
            u_allvoxels[xyblack_pix[i],,k] <- u_allvoxels[xyblack_pix[i],,k]
            IC_allvoxels[xyblack_pix[i],M,k] <- NA
            DF_allvoxels[xyblack_pix[i],M,k] <- NA
            SSE_allvoxels[xyblack_pix[i],M,k] <- NA
          }
        }
       }
    }
  }

Result_small_spat_lambda <- list(response_allvoxels=response_allvoxels,Ktrans_allvoxels=Ktrans_allvoxels,
                                 kep_allvoxels=kep_allvoxels,fit_allvoxels=fit_allvoxels,u_allvoxels=u_allvoxels,
                                 IC_allvoxels=IC_allvoxels,DF_allvoxels=DF_allvoxels,SSE_allvoxels=SSE_allvoxels)

M <- 10
numbercomp <- array(NA, c(I,J,K))
uniquekeps <- array(NA, c(I*J,M,K))
stoppingiter <- array(NA,c(I*J,K))

for(k in 1:K){ 
  
  kep_vals <- Result_small_spat_lambda$kep_allvoxels[,,k]
  ICs <- Result_small_spat_lambda$IC_allvoxels[,,k]
  
  for(i in 1:dim(kep_vals)[1]){
    x <- xysample[i,][[1]]
    y <- xysample[i,][[2]]
    
    w <- which.min(ICs[i,])
    keps <- kep_vals[i,1:w] 
    
    nkeps <- length(unique(keps))
    keps_un <- unique(keps)
    
    numbercomp[x,y,k] <- nkeps
    uniquekeps[i,1:nkeps,k] <- keps_un
  }
}

kepvals <- uniquekeps

unlink("refitspat-files", recursive = TRUE)

reg_refit = makeExperimentRegistry("refitspat", 
                                   packages = c(
                                     "ParamHelpers",
                                     "quadprog",
                                     "dcemriboost"
                                   )
)

addProblem(reg_refit, "refit", 
           static = list(xysample = xysample, Response = Response, Response_true = Response_true, I=I, J=J, time=time, kepvals=kepvals), 
           dynamic = function(static, pixel, image) {
             single_index <- static$xysample[pixel,] 
             list(
               y = static$Response[indx(single_index,static$I,static$J),,image],
               y_true = static$Response_true[indx(single_index,static$I,static$J),,image],
               keps = static$kepvals[indx(single_index,static$I,static$J),,image],
               image = image,
               pixel = pixel)
           }
)

addAlgorithm(reg_refit, id = "opt", fun = function(static, dynamic, ...) {
  
  ref <- refit(response=dynamic$y,response_true=dynamic$y_true,time=static$time,uniquekeps=dynamic$keps[which(!is.na(dynamic$keps))])
  
  return(list(len=ref$len,sse=ref$sse,sse_true=ref$sse_true,kk_esti=ref$kk_esti,image=dynamic$image,pixel=dynamic$pixel))  
})

pd_refit = makeDesign(id = "refit", exhaustive = list(image = 1:10, pixel = 1:2000))

addExperiments(reg_refit, prob.des = pd_refit)

IDs <- getJobIds(reg_refit)
ch <- chunk(IDs,n.chunks=10,shuffle=TRUE) 

submitJobs(reg_refit, ids=ch)
waitForJobs(reg_refit)

results_refit_list_allK <- list()
results_refit_list_allK[[1]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 1))
results_refit_list_allK[[2]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 2))
results_refit_list_allK[[3]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 3))
results_refit_list_allK[[4]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 4))
results_refit_list_allK[[5]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 5))
results_refit_list_allK[[6]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 6))
results_refit_list_allK[[7]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 7))
results_refit_list_allK[[8]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 8))
results_refit_list_allK[[9]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 9))
results_refit_list_allK[[10]] <- reduceResultsList(reg_refit,filterResults(reg_refit, fun=function(job, res) res$image == 10))

# average MSE, estimated average number of tissue compartments, 
# percentage of voxels for which the number of tissue compartments is correctly estimated

MSE_spat <- array(NA, c(I,J,K))

for(k in 1:length(results_refit_list_allK)){ 
  for(i in 1:length(results_refit_list_allK[[k]])){ 
    x <- xysample[i,][[1]]
    y <- xysample[i,][[2]]
    MSE_spat[x,y,k] <- (1/timep)*results_refit_list_allK[[k]][[i]]$sse_true 
  }

MSE_spat_vec <- vector("numeric",length=K)
for(k in 1:K){
  MSE_spat_vec[k] <- 1/(I*J) * sum(MSE_spat[,,k]) 
}

(mean(MSE_spat_vec))

av <- vector("numeric",length=K)
for(b in 1:K){
  tab <- table(numbercomp[,,b])
  avb <- 0
  for(n in 1:length(tab)){
    avb <- avb + tab[n]*n
  }
  av[b] <- avb/(I*J)
}

(av_number_comp_spat <- sum(av)/K)

sum_correct_boostspat <- 0

for(i in 1:K){
  w <- length(which(numbercomp[,,i] == NumberComp_true))
  sum_correct_boostspat <- sum_correct_boostspat + w
}
(percentage_correct_boostspat <- sum_correct_boostspat/(I*J*K))