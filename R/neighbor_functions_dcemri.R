#' Get scalar index from two-dimensional index
#' 
#' @param index two-dimensional index
#' @param I first dimension of array containing simulated measurements
#' @param J second dimension of array containing simulated measurements
#' 
#' @return indexsc scalar index

#' @export indx

indx <- function(index,I,J)
{ 
  indexsc <-  (index[1]-1)*J + index[2]
  return(indexsc)
}

getNeighbors <- function(p, P1, P2, Response)
{
  if(length(p)!=2){print("wrong dimension of pixel coordinates")
                   print(p)}
  
  xp <- p[1]
  yp <- p[2]
  
  dimx <- P1
  dimy <- P2
  
  neighbors <- list()
  count <- 1
  if(xp-1 > 0 && !is.na(Response[indx(c(xp-1,yp),dimx,dimy),1]))
  {
    neighbors[[count]] <- c(xp-1,yp)
    count <- count+1
  }
  if(xp+1 <= P1 && !is.na(Response[indx(c(xp+1,yp),dimx,dimy),1]))
  {
    neighbors[[count]] <- c(xp+1,yp)
    count <- count+1
  }
  if(yp-1 >0 && !is.na(Response[indx(c(xp,yp-1),dimx,dimy),1]))
  {
    neighbors[[count]] <- c(xp,yp-1)
    count <- count+1
  }
  if(yp+1 <= P2 && !is.na(Response[indx(c(xp,yp+1),dimx,dimy),1]))
  {
    neighbors[[count]] <- c(xp,yp+1)
    count <- count+1
  }
  
  return(neighbors)    
}


ISNA <- function(index, dimx, dimy, Response){
  
  ISNAResponse <- 0 
  
  if(is.na(Response[indx(index,dimx,dimy),1])){
    ISNAResponse <- 1
  }
    
  return(ISNAResponse)
}

NNB <- function(index, dimx, dimy, Response){
   
  NumberNB <- NA
    
  if(!is.na(Response[indx(index,dimx,dimy),1])){
    neighborList <- getNeighbors(index,dimx,dimy,Response)
    NumberNB <- length(neighborList)
  }
  
  return(NumberNB)
}

#' Get neighboring voxels
#' 
#' @param index two-dimensional index
#' @param dimx first dimension of array containing simulated measurements
#' @param dimy second dimension of array containing simulated measurements
#' @param Response array containing the measurements for all voxels
#' 
#' @return nb_thispix giving the indices of the neighboring voxels

#' @export WHICH_NB
WHICH_NB <- function(index, dimx, dimy, Response){
  
    neighborList <- getNeighbors(index,dimx,dimy,Response)
    nb_thispix <- NULL
    
    for(i in 1:length(neighborList)){
      nb_thispix <- c(nb_thispix,indx(neighborList[[i]],dimx,dimy))
    }
    
  return(nb_thispix)
}
