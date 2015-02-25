require(BatchExperiments)
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

unlink("simvw-files", recursive = TRUE)

reg = makeExperimentRegistry("simvw", 
                             packages = c(
                               "ParamHelpers",
                               "minpack.lm",
                               "dcemriboost"
                             )
)

addProblem(reg, "simul", 
           static = list(xysample = xysample, Response = Response, I=I, J=J, time=time), 
           
           dynamic = function(static, pixel, image) {
             single_index <- static$xysample[pixel,]
             list(y = static$Response[indx(single_index,static$I,static$J),,image],
                  image = image,
                  pixel = pixel)
           }
)

addAlgorithm(reg, id = "opt", fun = function(static, dynamic, ...) {
  
  model <- boostvw(response=dynamic$y,offset=NULL,nu=0.1,time=static$time,niter=100,factor=5) 
  
  list(nkeps=model$nkeps,keps_un=model$keps_un,ktrans=model$ktrans,keps=model$keps,image=dynamic$image,pixel=dynamic$pixel)
  
})

pd = makeDesign(id = "simul", exhaustive = list(image = 1:10, pixel = 1:2000))

addExperiments(reg, prob.des = pd)

IDs <- getJobIds(reg)
ch <- chunk(IDs,n.chunks=20,shuffle=TRUE) 

submitJobs(reg, ids=ch)
waitForJobs(reg)

results_list_allK <- list()
results_list_allK[[1]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 1))
results_list_allK[[2]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 2))
results_list_allK[[3]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 3))
results_list_allK[[4]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 4))
results_list_allK[[5]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 5))
results_list_allK[[6]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 6))
results_list_allK[[7]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 7))
results_list_allK[[8]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 8))
results_list_allK[[9]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 9))
results_list_allK[[10]] <- reduceResultsList(reg,filterResults(reg, fun=function(job, res) res$image == 10))

# save(results_list_allK,file="results_list_allK.RData")

numbercomp <- array(NA, c(I,J,K))
M <- 10 
uniquekeps <- array(NA, c(I*J,M,K))

res_factor5_allK <- list()

for(k in 1:K){
  for(i in 1:length(results_list_allK[[k]])){ 
    x <- xysample[i,][[1]]
    y <- xysample[i,][[2]]
    numbercomp[x,y,k] <- results_list_allK[[k]][[i]]$nkeps
    m <- length(results_list_allK[[k]][[i]]$keps_un)
    uniquekeps[i,1:m,k] <- results_list_allK[[k]][[i]]$keps_un
  }
}

# save(uniquekeps,file="uniquekeps.RData")
  
kepvals <- uniquekeps

unlink("refitvw-files", recursive = TRUE)

reg_refit = makeExperimentRegistry("refitvw", 
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
ch <- chunk(IDs,n.chunks=5,shuffle=TRUE) 

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

# save(results_refit_list_allK,file="results_refit_list_allK.RData")

# average MSE, estimated average number of tissue compartments, 
# percentage of voxels for which the number of tissue compartments is correctly estimated

MSE_voxelwise <- array(NA, c(I,J,K)) 

for(k in 1:length(results_refit_list_allK)){ 
  for(i in 1:length(results_refit_list_allK[[k]])){
    x <- xysample[i,][[1]]
    y <- xysample[i,][[2]]
    MSE_voxelwise[x,y,k] <- (1/timep)*results_refit_list_allK[[k]][[i]]$sse_true 
  }
}

MSE_vw <- vector("numeric",length=K)
for(k in 1:K){
  MSE_vw[k] <- 1/(I*J) * sum(MSE_voxelwise[,,k]) 
}

(mean(MSE_vw))

av <- vector("numeric",length=K)
for(b in 1:K){
  tab <- table(numbercomp[,,b])
  avb <- 0
  for(n in 1:length(tab)){
    avb <- avb + tab[n]*as.numeric(names(tab)[n])
  }
  av[b] <- avb/(I*J)
}

(average_number_comp <- sum(av)/K)

sum_correct <- 0

for(i in 1:K){
  w <- length(which(numbercomp[,,i] == NumberComp_true))
  sum_correct <- sum_correct + w
}
(percentage_correct <- sum_correct/(I*J*K))