## © Asger Hobolth, 2017

##------------------------------------------------------------------
## Time discretization:
## Building the Markov chain based on the joint description
##------------------------------------------------------------------
# setwd("C:/Users/asger/Projects/ParticleFilter/Programs/")
## Need expm package: install.packages("expm")
library("expm")
##------------------------------------------------------------------
##----------- Simonsen-Churchill rate matrix  ----------------------
##------------------------------------------------------------------
## Name: ARGRateM
## Input: rho: recombination rate
## Output: 8x8 Simonsen-Churchill rate matrix
## Comments: This is the rate matrix in Fig 1a in
##           Hobolth and Jensen (2014, TPB)
##------------------------------------------------------------------
ARGRateM <- function(rho){
  RateM <-
    matrix(c(0,rho,0,0,0,0,0,1,
             1,0,rho/2,1,1,0,0,0,
             0,4,0,0,0,1,1,0,
             0,0,0,0,0,rho/2,0,1,
             0,0,0,0,0,0,rho/2,1,
             0,0,0,2,0,0,0,1,
             0,0,0,0,2,0,0,1,
             0,0,0,0,0,0,0,0),nrow=8,ncol=8,byrow=TRUE)
  ## Get diagonals right (must sum to 0)
  for (rw in 1:8){
    RateM[rw,rw] <- -sum(RateM[rw,])
  }
  return(RateM)
}
##-------------------------------------------------------
## Build transition matrix
##------------------------------------------------------
## Name: TransMat
## Input: tm: Discrete times
##        rho: recombination rate
## Output: Transition matrix
TransMat <- function(tm,rho){
  ARGRateMat <- ARGRateM(rho)
  nInt <- length(tm) ## Number of intervals
  tm0 <- c(0,tm) ## tm0 is tm with time 0 added
  ##-------------------------------------
  JointMat <- matrix(0,nrow=nInt,ncol=nInt) ## Joint prb matrix
  #cat("Start calculating joint matrix...","\n")
  for (j in 1:(nInt-1)){  ## Left state
    #cat(j)
    for (k in j:nInt){  ## Right state
      if (j<k){
        JointMat[j,k] <-
          expm(tm0[j]*ARGRateMat)[1,1:3]%*%
          expm((tm[j]-tm0[j])*ARGRateMat)[1:3,c(5,7)]%*%
          expm((tm0[k]-tm[j])*ARGRateMat)[c(5,7),c(5,7)]%*%
          if (k==nInt) c(1,1) else
            expm((tm[k]-tm0[k])*ARGRateMat)[c(5,7),8]
        ## (expm can't handle infinity, and state 8 is absorbing)
        ## Symmetrize
        JointMat[k,j] <- JointMat[j,k]
      }
      if (j==k){
        JointMat[j,k] <-
          expm(tm0[j]*ARGRateMat)[1,1:3]%*%
          expm((tm[j]-tm0[j])*ARGRateMat)[1:3,8]
      }
    }
  }
  #cat("\n")
  ## Final entry
  JointMat[nInt,nInt] <- sum(expm(tm0[nInt]*ARGRateMat)[1,1:3])
  ## Again: expm can't handle infinity, and state 8 is absorbing
  #cat("Finished calculating joint rate matrix","\n")
  ## Transition matrix
  TransMat <- JointMat/rowSums(JointMat)
  return(TransMat)
}
##-------------------------------------------------------
## Fast version of transition matrix
##-------------------------------------------------------
## Name: TransMat
## Input: tm: Discrete times
##        rho: recombination rate
## Output: Transition matrix
FastTransMat <- function(tm,rho){
  #ARGRateMat <- ARGRateM(rho)
  LumpFromFour <- matrix(c(-(1+rho),rho,0,1,
                           1,-(3+rho/2),rho/2,2,
                           0,4,-6,2,
                           0,0,0,0),
                         nrow=4,ncol=4,byrow=TRUE)
  LumpFourToSeven <- matrix(c(-(1+rho),rho,0,0,1,
                              1,-(3+rho/2),rho/2,2,0,
                              0,4,-6,2,0,
                              0,0,0,-1,1,
                              0,0,0,0,0),
                            nrow=5,ncol=5,byrow=TRUE)
  nInt <- length(tm) ## Number of intervals
  tm0 <- c(0,tm) ## tm0 is tm with time 0 added
  ##-------------------------------------
  JointMat <- matrix(0,nrow=nInt,ncol=nInt) ## Joint prb matrix
  #cat("Start calculating joint matrix...","\n")
  for (j in 1:(nInt-1)){  ## Left state
    #cat(j)
    for (k in j:nInt){  ## Right state
      if (j<k){
        JointMat[j,k] <-
          #expm(tm0[j]*ARGRateMat)[1,1:3]%*%
          expm(tm0[j]*LumpFromFour)[1,1:3]%*%
          #expm((tm[j]-tm0[j])*ARGRateMat)[1:3,c(4,5,6,7)]%*%
          #rowSums(expm((tm[j]-tm0[j])*ARGRateMat)[1:3,c(4,5,6,7)])*
          expm((tm[j]-tm0[j])*LumpFourToSeven)[1:3,4]*
          #rep(exp(-(tm0[k]-tm[j])),4)*0.5*
          exp(-(tm0[k]-tm[j]))*0.5*
          #if (k==nInt) c(1,1) else
            #expm((tm[k]-tm0[k])*ARGRateMat)[c(5,7),8]
            (1-exp(-(tm[k]-tm0[k])))
        ## (expm can't handle infinity, and state 8 is absorbing)
        ## Symmetrize
        JointMat[k,j] <- JointMat[j,k]
      }
      if (j==k){
        JointMat[j,k] <-
          #expm(tm0[j]*ARGRateMat)[1,1:3]%*%
          expm(tm0[j]*LumpFromFour)[1,1:3]%*%
          #expm((tm[j]-tm0[j])*ARGRateMat)[1:3,8]
          expm((tm[j]-tm0[j])*LumpFourToSeven)[1:3,5]
      }
    }
  }
  #cat("\n")
  ## Final entry
  #JointMat[nInt,nInt] <- sum(expm(tm0[nInt]*ARGRateMat)[1,1:3])
  JointMat[nInt,nInt] <- sum(expm(tm0[nInt]*LumpFromFour)[1,1:3])
  ## Again: expm can't handle infinity, and state 8 is absorbing
  #cat("Finished calculating joint rate matrix","\n")
  ## Transition matrix
  TransMat <- JointMat/rowSums(JointMat)
  return(TransMat)
}
##----------------------------------------------------------------
## Debugging...
## Define time intervals, recombination rate and transition matrix
#nInt <- 10 ; tm <- -log(1-1:(nInt)/nInt)
#rho <- 0.3
#TransM <- TransMat(tm,rho)
#FTransM <- FastTransMat(tm,rho)
