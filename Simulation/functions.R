library(psych)
library(cvTools)
library(parallel)
library(tidyverse)
library(svd)
library(foreach)
kl <- function(beta,ui, test){
  calculated <- beta%*%t(beta) + diag(ui)
  calculated <- as.matrix(calculated)
  cov_f <- cor(test)
  KL <- log(det(calculated)) + 
    tr(solve(calculated)%*%cov_f)-
    log(det(cov_f)) - 
    ncol(test)
  KL <- KL*0.5
  return(KL)
}

qBest <- function(train, val, nfactor){
  kl <- tryCatch( {result <-  factanal(train, factors = nfactor, rotation = "varimax")
                  beta <- result$loadings
                  ui <- result$uniquenesses
                  kl(beta, ui, val)},
             error = function(e){
               print(e)
               print(nfactor)
               return(NA)
             }
  )
  #result <- fa(train, nfactors = nfactor, rotate = "varimax")
  return(kl)
}
#qBest_prime <- function(train,train_cltr, val,val_cltr, nfactor)

cv_evaluate <- function(nfolds, func, data, ncores = 1,param, ...){
  
  if (nfolds <= 1){
    stop("fold number must larger than 1")
  }
  #func <- match.fun(func)
  cl <- makeCluster(ncores)
  registerDoSEQ()
  folds <- cvFolds(NROW(data), K=nfolds)

  r <- foreach(i=1:nfolds,.combine = rbind) %dopar% {
    train <- data[folds$subsets[folds$which != i], ]
    validation <- data[folds$subsets[folds$which == i], ]
    result <- map_dbl(param, func, train = train, val = validation,...)
    return(result)
    stopCluster(cl)
  }
  stopCluster(cl)
  return(r)
}

nonZeroLoad <- function(load){
  load <- map(1:ncol(load), function(x) x[x != 0])
  return(load)
}


interSNP <- function(dat){
  SNPs <- attr(dat, "interact")
  return(SNPs)
}

tuneSfa <- function(train, val, nfactor){
  result <- fanc::fanc(train,
                       factors = nfactor, 
                       control = list(openmp = TRUE, num.threads = 8))
  lists <- c(result$loadings)
  
}


