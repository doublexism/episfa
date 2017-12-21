library(psych)
library(cvTools)
library(doParallel)
library(tidyverse)
library(svd)
library(foreach)
kullback <- function(beta,ui, test){
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
                  kullback(beta, ui, val)},
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
  cl <- makeCluster(ncores, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
  registerDoParallel(cl)
  folds <- cvFolds(NROW(data), K=nfolds)

  r <- foreach(i=1:nfolds,.combine = rbind) %dopar% {
    train <- data[folds$subsets[folds$which != i], ]
    validation <- data[folds$subsets[folds$which == i], ]
    result <- purrr::map_dbl(param, func, train = train, val = validation,...)
    return(result)
  }
  
  stopCluster(cl)
  return(r)
}

nonZeroLoad <- function(load, coef = TRUE){
  loadings <- plyr::alply(load, 2, function(x) x[x != 0])
  if (!coef){
  nonzeros <- loadings[lengths(loadings) !=0]
  SNPs <- map(nonzeros, names)
  return(SNPs)
  }
  return(loadings)
}


interSNP <- function(dat){
  SNPs <- attr(dat, "interact")
  return(SNPs)
}

tuneSfa <- function(train, val, nfactor){
  result <- fanc::fanc(train,
                       factors = nfactor)
  lists <- c(result$loadings)
}

cv.consistency <- function(episfa.obj, stat){
  snps <- map(episfa.obj, `[[`,stat)
  cvc <- map_
  
  
}

simPopLE_l2_sp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, m_eff = NULL){
          func <- simPopLE(n,
                     num_SNP = snp_num,
                     MAF = maf,  
                     main_effect = 1.5, 
                     interaction_effect =int_eff, 
                     margin_effect = m_eff,
                     cov_effect = 1.2, 
                     level = int_lev,
                     num_parents = 0, 
                     num_sib = 2, 
                     num_main = 0,
                     num_interact = int_num,
                     model = "logistic",
                     genetic_model = "additive",
                     p = p,
                     scale_weibull = 80,
                     shape_weibull = 4,
                     age_lower = 20,
                     age_higher = 80,
                     sex_effect = 1.2,
                     age_effect = 1.005,
                     age_mean = 50,
                     age_varb = 10,
                     age_varw = 5,
                     num_cov = 10,
                     cov_sigma = NULL,
                     population = FALSE)
          return(func)
}


simPopLE_l2_sp_wp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, m_eff = NULL){
  func <- simPopLE(n,
                   num_SNP = snp_num,
                   MAF = maf,  
                   main_effect = 1.5, 
                   interaction_effect =int_eff, 
                   margin_effect = m_eff,
                   cov_effect = 1.2, 
                   level = int_lev,
                   num_parents = 2, 
                   num_sib = 2, 
                   num_main = 0,
                   num_interact = int_num,
                   model = "logistic",
                   genetic_model = "additive",
                   p = p,
                   scale_weibull = 80,
                   shape_weibull = 4,
                   age_lower = 20,
                   age_higher = 80,
                   sex_effect = 1.2,
                   age_effect = 1.005,
                   age_mean = 50,
                   age_varb = 10,
                   age_varw = 5,
                   num_cov = 10,
                   cov_sigma = NULL,
                   population = FALSE)
  return(func)
}


interForm <- function(SNPs, data, model = "clogit"){
  require(survival,quietly = TRUE)
  main <- paste0(SNPs, collapse = "+")
  inter <- paste0(SNPs, collapse = ":")
  if (model == "clogit"){
  form <- paste0("Y ~",main, "+",inter,"+ strata(fid)") %>% 
    as.formula()
  res <- summary(clogit(form, data))
  } else if (model == "glm") {
    form <- paste0("Y~", main, "+", inter) %>% 
      as.formula()
    res <- summary(form, data, family = binomial(link = "logit"))
  } else if (model == "lm") {
    y <- paste0(SNPs[1],"~")
    if (length(SNPs) > 2){
    inter <-  paste0(SNPs[-1], collapse = ":")
    form <- paste0(y ,inter) %>% as.formula()
    } else {
    form <- paste0(y ,SNPs[2]) %>% as.formula()
    }
    res <- summary(lm(form, data))
  }
  return(res)
}

simVal <- function(dat, subset = NULL){
  inters <- interSNP(dat)
  if (is.null(subset)){
    subset <- 1:nrow(dat)
  }
  dat <- dat[subset,]
  models <- c("clogit", "lm")
  result_cc <- map(inters, interForm, data = dat, model = "clogit")
  result_co <- map(inters, interForm, data = dat, model = "lm")
  return(list(cc= result_cc, co = result_co))
}

episfa <- function(x, nfolds, type = "data", control = NULL){
  # function to extract best gamma and rho
  best_gr <- function(stat, nzl, n = 270){
    stat <- stat[1:n]
    nzl <- nzl[1:n]
    min_stat <- which(stat == min(stat))
    min_nzl <- min(nzl[min_stat])
    best <- min_stat[nzl[min_stat] == min_nzl][1]
    rho <- best %% 30
    if (rho == 0){
      gamma <- best %/% 30
      rho = 30
    } else {
      gamma <- best %/% 30 + 1
    }
    return(c(gamma, rho))
  }
  
  # function to extract best loadings
  best_loading <- function(x, gr){
      gr1 <- x$gamma[gr[1]]
      gr2 <- x$rho[gr[2],gr[1]]
      return(fanc::out(x,gamma = gr1, rho = gr2))
  }
  
  if(!is.matrix(x)){
    stop("x should be a matrix")
  }
  
  folds <- cvFolds(nrow(x), K=nfolds)
  # cross validation
  r <- foreach(i=1:nfolds) %do% {
    print(paste0("round",i))
    train <- x[folds$subsets[folds$which != i], ]
    validation <- x[folds$subsets[folds$which == i], ]
    if (type == "data"){
      result <- fanc(train, factors = round(0.05*ncol(x)))
    } else if(type == "covmat"){
        if (!is.null(control)){
          cov_mat <- cor(train) - cor(control) + diag(ncol(x))
        } else {
          cov_mat <- cor(train)
        }
      result <- fanc(factors = round(0.1*ncol(x)), covmat = cov_mat, n.obs = nrow(train))
    } else if(type == "poly"){
      result <- fanc(factors = round(0.1*ncol(x)), covmat = psych::polychoric(x)$rho, n.obs = nrow(train))
    } else {
      stop("type must be data or covmat or poly")
    }
    # extract statistics
    loadings <- reduce(result$loadings, c)
    ui <- plyr::alply(result$uniquenesses,c(1,3),`[`)
    kl <- map2_dbl(loadings, ui, kullback, test = validation)
    aic <- as.vector(result$AIC_dfnonzero)
    bic <- as.vector(result$BIC_dfnonzero)
    caic <- as.vector(result$CAIC_dfnonzero)
    nonzero <- as.vector(result$nonzero.loadings)
    # extract interactions
    best_kl <- best_gr(kl, nonzero)
    best_aic <- best_gr(aic, nonzero)
    best_caic <- best_gr(caic, nonzero)
    best_bic <- best_gr(bic, nonzero)
    best_kl_lasso <- best_gr(kl, nonzero, 30)
    best_caic_lasso <- best_gr(caic, nonzero, 30)
    best_aic_lasso <- best_gr(aic, nonzero, 30)
    best_bic_lasso <- best_gr(bic, nonzero, 30)
    # extract loadings
    loading_kl <- best_loading(result, best_kl)
    loading_aic <- best_loading(result, best_aic)
    loading_caic <- best_loading(result, best_caic)
    loading_bic <- best_loading(result, best_bic)
    loading_kl_lasso <- best_loading(result, best_kl_lasso)
    loading_caic_lasso <- best_loading(result, best_caic_lasso)
    loading_aic_lasso <- best_loading(result, best_aic_lasso)
    loading_bic_lasso <- best_loading(result, best_bic_lasso)
    # extract non zero loadings
    nz_kl <- nonZeroLoad(loading_kl$loadings, coef = FALSE)
    nz_aic <- nonZeroLoad(loading_aic$loadings, coef = FALSE)
    nz_caic <- nonZeroLoad(loading_caic$loadings, coef = FALSE)
    nz_bic <- nonZeroLoad(loading_bic$loadings, coef = FALSE)
    nz_kl_lasso <- nonZeroLoad(loading_kl_lasso$loadings, coef = FALSE)
    nz_caic_lasso <- nonZeroLoad(loading_kl_lasso$loadings, coef = FALSE)
    nz_aic_lasso <- nonZeroLoad(loading_aic_lasso$loadings, coef = FALSE)
    nz_bic_lasso <- nonZeroLoad(loading_bic_lasso$loadings, coef = FALSE)
    # return results
    return(list(kl  = kl,
                nonzero = nonzero,
                nz_kl = nz_kl, 
                nz_aic=nz_aic, 
                nz_caic = nz_caic,
                nz_bic = nz_bic, 
                nz_kl_lasso = nz_kl_lasso,
                nz_caic_lasso = nz_caic_lasso,
                nz_aic_lasso= nz_aic_lasso, 
                nz_bic_lasso = nz_bic_lasso,
                loading_kl = loading_kl,
                loading_aic = loading_aic,
                loading_caic = loading_caic,
                loading_bic = loading_bic,
                loading_kl_lasso = loading_kl_lasso,
                loading_aic_lasso = loading_aic_lasso,
                loading_bic_lasso = loading_bic_lasso))
  }
  return(r)
}
  
  episfa_sim <- function(co = FALSE , sim_control){}

svd_compress <- function(dat, depth = NULL){
  n <- nrow(dat)
  p <- ncol(dat)
  dat <- scale(dat)
  r <- svd(dat)
  if (is.null(depth)){
      var_p <- r$d**2/sum(r$d**2)*p
      depth <- sum(var_p > 1)
  }
  dat_comp <- r$u[,1:depth] %*% diag(r$d[1:depth]) %*% t(r$v[,1:depth])
  v_preserve <- diag(var(dat_comp))
  variance <- 1 - v_preserve
  dat_comp <- dat_comp + MASS::mvrnorm(n, rep(0,p), diag(variance))
  colnames(dat_comp) <- colnames(dat) 
  return(dat_comp)
}

cov_diff <- function(co, ctrl){
  co <- cor(co)
  ctrl <- cor(ctrl)
  ctrl[ctrl ==1 | abs(ctrl) < 0.05] <- 0
  co <- co - ctrl
  return(co)
}
