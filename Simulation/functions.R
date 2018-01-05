library(psych)
library(cvTools)
library(doParallel)
library(tidyverse)
library(svd)
library(foreach)
library(fanc)
tr <- psych::tr

log_trans <- function(Mat, offset){
  return(log(Mat + offset))
}

kullback <- function(beta,ui, test){
  calculated <- beta%*%t(beta) + diag(ui)
  calculated <- as.matrix(calculated)
  cov_f <- t(scale(test)) %*% scale(test)
  KL <- log(det(calculated)) + 
    tr(solve(calculated)%*%cov_f)
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

simPopLE_l2_sp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, m_eff = NULL, mode = NULL){
          func <- simPopLE(n,
                     num_SNP = snp_num,
                     MAF = maf,  
                     main_effect = 1.5, 
                     interaction_effect =int_eff, 
                     margin_effect = m_eff,
                     inter_mode = mode,
                     cov_effect = 1, 
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
                     sex_effect = 1,
                     age_effect = 1,
                     age_mean = 50,
                     age_varb = 10,
                     age_varw = 5,
                     num_cov = 10,
                     cov_sigma = NULL,
                     population = FALSE)
          return(func)
}


simPopLE_l2_sp_wp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, m_eff = NULL, mode = NULL){
  func <- simPopLE(n,
                   num_SNP = snp_num,
                   MAF = maf,  
                   main_effect = 1.5, 
                   interaction_effect =int_eff, 
                   margin_effect = m_eff,
                   inter_mode = mode,
                   cov_effect = 1.2, 
                   level = int_lev,
                   num_parents = 2, 
                   num_sib = 1, 
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
  result_co <- map(inters, interForm, data = dat[dat$Y == 1,], model = "lm")
  return(list(cc= result_cc, co = result_co))
}

cv.episfa <- function(x, nfolds,nfactor = NULL,sparsity = 0.05, type = "data", contrast = NULL, ...){
 
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
  if (is.null(nfactor)){
    nfactor <- sparsity*ncol(x)
  } 
  print(paste0("The size of matrix is ", dim(x)))
  print(paste0("The number of factors is ", nfactor))
  folds <- cvFolds(nrow(x), K=nfolds)
  # cross validation
  r <- foreach(i=1:nfolds) %do% {
    print(paste0("cross validate round",i))
    timestamp()
    train <- x[folds$subsets[folds$which != i], ]
    validation <- x[folds$subsets[folds$which == i], ]
    if (type == "data"){
      result <- fanc(train, factors = nfactor, ...)
    } else if(type == "covmat"){
        if (!is.null(contrast)){
          cov_mat <- cor(train) - cor(contrast) + diag(ncol(train))
        } else {
          cov_mat <- cor(train)
        }
      result <- fanc(factors = nfactor, covmat = cov_mat, n.obs = nrow(train),...)
    } else if(type == "poly"){
      result <- fanc(factors = nfactor, covmat = psych::polychoric(train)$rho, n.obs = nrow(train),...)
    } else {
      stop("type must be data or covmat or poly")
    }
    # extract statistics
    loadings <- reduce(result$loadings, c)
    ui <- plyr::alply(result$uniquenesses,c(1,3),`[`)
    print("begin to calc kl")
    kl <- map2_dbl(loadings, ui, kullback, test = validation)
    print("finished")
    aic <- as.vector(result$AIC_dfnonzero)
    bic <- as.vector(result$BIC_dfnonzero)
    caic <- as.vector(result$CAIC_dfnonzero)
    nonzero <- as.vector(result$nonzero.loadings)
    df <- as.vector(result$dfnonzero)
    # extended and high-dimentional bic
    ebic0.5 <- EBIC(bic, p = result$factors * ncol(train),df = df, gamma = 0.5)
    ebic0.75 <- EBIC(bic, p = result$factors *ncol(train), df = df, gamma = 0.75)
    ebic1 <- EBIC(bic, p = result$factors *ncol(train), df = df, gamma = 1)
    hbic0.5 <- HBIC(bic, p = result$factors * ncol(train),df = df, gamma = 0.5)
    hbic0.75 <- HBIC(bic, p = result$factors *ncol(train), df = df, gamma = 0.75)
    hbic1 <- HBIC(bic, p = result$factors *ncol(train), df = df, gamma = 1)
    # extract interactions
    best_kl <- best_gr(kl, nonzero)
    best_aic <- best_gr(aic, nonzero)
    best_caic <- best_gr(caic, nonzero)
    best_bic <- best_gr(bic, nonzero)
    best_ebic0.5 <- best_gr(ebic0.5, nonzero)
    best_ebic0.75<- best_gr(ebic0.75, nonzero)
    best_ebic1<- best_gr(ebic1, nonzero)
    best_hbic0.5<- best_gr(hbic0.5, nonzero)
    best_hbic0.75<- best_gr(hbic0.75, nonzero)
    best_hbic1<- best_gr(hbic1, nonzero)
    # extract loadings
    loading_kl <- best_loading(result, best_kl)
    loading_aic <- best_loading(result, best_aic)
    loading_caic <- best_loading(result, best_caic)
    loading_bic <- best_loading(result, best_bic)
    loading_ebic0.5 <- best_loading(result, best_ebic0.5)
    loading_ebic0.75 <- best_loading(result, best_ebic0.75)
    loading_ebic1 <- best_loading(result, best_ebic1)
    loading_hbic0.5 <- best_loading(result, best_hbic0.5)
    loading_hbic0.75 <- best_loading(result, best_hbic0.75)
    loading_hbic1 <- best_loading(result, best_hbic1)
    # extract non zero loadings
    nz_kl <- nonZeroLoad(loading_kl$loadings, coef = FALSE)
    nz_aic <- nonZeroLoad(loading_aic$loadings, coef = FALSE)
    nz_caic <- nonZeroLoad(loading_caic$loadings, coef = FALSE)
    nz_bic <- nonZeroLoad(loading_bic$loadings, coef = FALSE)
    nz_ebic0.5 <- nonZeroLoad(loading_ebic0.5$loadings, coef = FALSE)
    nz_ebic0.75 <- nonZeroLoad(loading_ebic0.75$loadings, coef = FALSE)
    nz_ebic1 <- nonZeroLoad(loading_ebic1$loadings, coef = FALSE)
    nz_hbic0.5 <- nonZeroLoad(loading_hbic0.5$loadings, coef = FALSE)
    nz_hbic0.75 <- nonZeroLoad(loading_hbic0.75$loadings, coef = FALSE)
    nz_hbic1 <- nonZeroLoad(loading_hbic1$loadings, coef = FALSE)
    # return results
    return(list(kl  = matrix(kl, nrow  = 30),
                bic = matrix(bic,nrow = 30),
                ebic0.5 = matrix(ebic0.5,nrow =30),
                ebic0.75 = matrix(ebic0.75, nrow = 30),
                ebic1 = matrix(ebic1, nrow = 30),
                hbic0.5 = matrix(hbic0.5, nrow = 30),
                hbic0.75 = matrix(hbic0.75, nrow = 30),
                hbic1 = matrix(hbic1, nrow = 30),
                nonzero = matrix(nonzero,nrow = 30),
                factors = result$factors,
                time = result$time,
                nz_kl = nz_kl, 
                nz_aic=nz_aic, 
                nz_caic = nz_caic,
                nz_bic = nz_bic, 
                nz_ebic0.5 = nz_ebic0.5,
                nz_ebic0.75 = nz_ebic0.75,
                nz_ebic1 = nz_ebic1,
                nz_hbic0.5 = nz_hbic0.5,
                nz_hbic0.75 = nz_hbic0.75,
                nz_hbic1 = nz_hbic1,
                loading_kl = loading_kl,
                loading_aic = loading_aic,
                loading_caic = loading_caic,
                loading_bic = loading_bic,
                loading_ebic0.5 = loading_ebic0.5,
                loading_ebic0.75 = loading_ebic0.75,
                loading_ebic1 = loading_ebic1,
                loading_hbic0.5 = loading_hbic0.5,
                loading_hbic0.75 = loading_hbic0.75,
                loading_hbic1 = loading_hbic1))
  }
  return(r)
}

isInter<- function(inter1, inter.list){
  consist <- map_dbl(inter.list, ~any(inter1 %in% .)) %>% sum
  return(consist)
}

cv.consistency <- function(episfa.obj, stat){
  inters.list <- map(episfa.obj, `[`, stat) %>% 
    unlist(recursive = FALSE) %>%
    unlist(recursive = FALSE)
  inters <- unique(inters.list) 
  inter_name <- map_chr(inters, paste0, collapse = "")
  consistency <- map_dbl(inters, isInter, inters.list) %>% 
    setNames(inter_name) %>%
    sort(decreasing = TRUE)
  return(consistency)
}

effRemove <- function(consistency, data){
  best <- names(consistency[[1]][1])
  snps <- str_extract_all(best,"SNP[0-9]*") %>% unlist() 
  cols <- setdiff(colnames(data), snps)
  dat <- data[,cols]
  return(dat)
}

episfa <- function(dat, nfolds, recursion = 5, criteria = "ebic",...){
  if (criteria  %in% c("ebic","hbic")){
    criteria <- paste0("nz_",criteria,c(1, 0.75, 0.5))
  } else {
    criteria <-  paste0("nz_",criteria)
  }
  inters <- character()
  ## recursion
  for(i in 1:recursion) {
    print(paste0("Searching for epistatic effects: number ",i,", ..."))
    result <- cv.episfa(dat, nfolds, nfactor = 1, type = "data",...)
    consistency <- map(1:length(criteria),~cv.consistency(result, criteria[.]))
    # consistency <- consist
    elements <- lengths(consistency)
    if (sum(elements) >= 1) {
      effect <- consistency[elements >= 1][[1]][1]
      name <- names(effect)
      print(paste0("Found interaction ",name))
      dat <- effRemove(consistency, dat)
       if (is.null(dat) || ncol(dat) < 2){
         break()
       }
      inters[name] <- effect
    } else {
      sprintf("no interaction found for criteria %s, algorithm halts", criteria)
      inters[i] <- NA
      break()
    }
  }
  return(inters)
}
    
episfa_sim <- function(n_rep = 100, recursion = 5, cvfolds = 10, ncores = NULL,sim_func = simPopLE_l2_sp, sim_control = list(), criteria = "ebic",...){
    if (is.null(ncores)){
      ncores <- detectCores(logical = FALSE)
    }
    sprintf("running on %i cores",ncores) %>% 
      print()
    ## make cluster
    cl <- makeCluster(ncores, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
    registerDoParallel(cl)
    ## hyperparameters
    num_interact <- sim_control[["int_num"]]
    global_funcs <- lsf.str(.GlobalEnv) %>% as.vector()
    ## looping
    benchmark <- foreach(i = 1:n_rep, 
                         .export = global_funcs, 
                         .packages = c("purrr","stringr","fanc","dplyr","sigmoid","truncnorm","data.table","cvTools", "foreach"),
                         .verbose = TRUE) %dopar% {
      #timing
      time_start <- Sys.time()
      # simulate data
      simdata <- do.call(sim_func,sim_control)
      inters <- interSNP(simdata) %>% map_chr(paste0, collapse = "") 
      snp_name <- colnames(simdata) %>% str_subset('^SNP[0-9]+$')
      df_co <- simdata[Y == 1,snp_name, with = FALSE] %>% as.matrix()
      # episfa runsss
      result_episfa <- episfa(df_co, cvfolds, recursion, criteria, ...)
      # inters <- c("SNP1SNP2","SNP3SNP4","SNP5SNP6")
      # result_episfa <- list(SNP1SNP2 = 3, SNP3SNP4 = 2, SNP5SNP7 = 1)
      ## false_positive 1: any unknown interaction effect is false dicovery
      false_positive <- 0
      false_positive_any <- 0
      false_positive_all <- 0
      ## true_positive 2: the proportion of false discovered effects among all discovered effects
      true_positive <- 0
      true_positive_any <- 0
      true_positive_all <- 0
      
      inter_names <-names(na.omit(result_episfa))
      if (!is.null(inter_names)){
        false_positive <- map_dbl(inter_names, interDiscover,inters) %>%
          sum()
        false_positive_any <- ifelse(false_positive > 0, 1, 0)
        false_positive_all <- ifelse(false_positive == length(inter_names), 1, 0)
        true_positive <- length(inter_names) - false_positive
        # any interaction were discovered
        true_positive_any <- ifelse(true_positive > 0 , 1, 0)
        # all discovered interaction are positive 
        true_positive_all <- ifelse(true_positive == num_interact, 1, 0)
      }
      time_end <- Sys.time()
      time_diff <- round(as.numeric(time_end - time_start),1)
      return(list(fp = false_positive,
             fpany = false_positive_any,
             fpall = false_positive_all,
             tp = true_positive,
             tpany = true_positive_any,
             tpall = true_positive_all,
             result = result_episfa,
             true_inter = inters,
             time = time_diff))
    }
    fp <- getListElement(benchmark,"fp") %>% sum()
    fpany <- getListElement(benchmark,"fpany") %>% sum()
    fpall <- getListElement(benchmark,"fpall") %>% sum()
    tp <- getListElement(benchmark,"tp") %>% sum()
    tpany <- getListElement(benchmark,"tpany") %>% sum()
    tpall <- getListElement(benchmark,"tpall") %>% sum()
    # alpha, power and ba
    alpha1 <- fpany/n_rep
    alpha2 <- fpall/n_rep
    alpha3 <- fp/(n_rep * recursion)
    power1 <- tpany/n_rep
    power2 <- tpall/n_rep
    if (num_interact != 0 ){
      power3 <- tp/(n_rep * num_interact)
      ba3 <- (1-alpha3 + power3)/2
    } else {
      ba3 <- NULL
    }
    ba1 <- (1-alpha1 + power1)/2
    ba2 <- (1-alpha2 + power2)/2
    stopCluster(cl)
    return(list(
      alpha1 = alpha1,
      alpha2 = alpha2,
      alpha3 = alpha3,
      power1 = power1,
      power2 = power2,
      power3 = power3,
      ba1 = ba1,
      ba2 = ba2,
      ba3 = ba3,
      benchmark = benchmark))
}

interDiscover <- function(candidate, inters){
  if (!is.null(inters)){
    in_inters <- str_detect(candidate, inters) %>%
      sum() 
    return(ifelse(in_inters == 0, 1 ,0))
  } else {
    return(1)
  }
}

getListElement <- function(.l, name, simplify = TRUE){
  element <- map(.l,`[[`,name)
  if (simplify){
    element <- unlist(element)
  }
  return(element)
}

simResults <- function(sim_control, sfa_control, n_rep = 100, recursion = 2, cvfolds = 5,ncores = NULL, sim_func = simPopLE_l2_sp, criteria = "ebic", save = TRUE){
  sim_param <- sim_control %>% 
    setNames(c("n", "snp_num", "maf","p","int_eff","int_lev","int_num")) %>%
    as.list
  scene <- episfa_sim(n_rep, 
                    recursion,
                    cvfolds, 
                    ncores,
                  sim_func,
                   sim_param, 
                   criteria,
                    sfa_control)
  name <- pmap_chr(sim_param, paste, sep = '-')
  if (save == TRUE){
    write_rds(scene,paste0("results/intern1p05", name,".rds"))
  }
  return(scene)
}

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

EBIC <- function(BIC, p, df, gamma){
  return(BIC + 2*gamma*log(p)*df)
}

HBIC <- function(BIC, p, df, gamma){
  return(BIC - log(p)*df + 2*gamma*log(p)*df)
}

partial <- function(y, X){
  X <- cbind(rep(1,nrow(X)),X)
  hat <- X %*% solve((t(X)%*%X)) %*% t(X) %*% y 
  res <- y - hat
  return(as.vector(res))
}


prodCol2 <- function(x,y = NULL, func = `*`){
  if (is.null(y)){
    y <- x
  }
  if (nrow(x) != nrow(y)){
    stop("x and y in 'prodCol' function ust have same number of rows")
  }
  if (is.data.frame(x) | is.data.frame(y)){
    x <- as.matrix(x)
    y <- as.matrix(y)
  }
  product <- plyr::alply(x, 2,function(col) col*y)
  product <- product%>% do.call(cbind,.)
  colnames(product) <- NULL
  return(product)
}

prodCol <- function(l, func = `*`){
  if (is.matrix(l)){
    return(prodCol2(l))
  }
  product <- reduce(l,prodCol2, func)
  return(product)
}

genoType <- function(vec){
  num <- length(vec)
  mgeno <- vec[seq(1,num,2)]
  pgeno <- vec[seq(2,num,2)]
  pargeno <- 3*mgeno+pgeno
  return(pargeno)
}

stanSurf <- function(sib, pat, type = "pat"){
  mean <- setDT(data.frame("sib" = sib, "pat" = pat))[,.(m = mean(sib)), pat]
  m <- mean$m
  pat_type <-mean$pat
  mean <- m[match(pat, pat_type)]
  return(mean)
}

corrS <- function(sib_geno, parent_geno){
  parental_type <- map(parent_geno, genoType) 
  sib_sub <- map2_dfc(sib_geno, parental_type, stanSurf)
  sib_geno <- sib_geno - sib_sub
  sib_geno <- as.matrix(sib_geno)
  X <- t(sib_geno) %*% sib_geno /nrow(sib_geno)
  return(X)
}

permuteMatrix <- function(mat, margin = 2){
  num_dim <- dim(mat)[margin]
  nums <- sample.int(num_dim,num_dim,replace = FALSE)
  mat <- mat[,nums]
  return(mat)
}


# sib_geno <- test1200[!is.na(mid),] %>% arrange(fid)
# sib_geno <- sib_geno[,1:200]
# pat_geno <- test1200[is.na(mid),] %>% arrange(fid)
# pat_geno <- pat_geno[,1:200]

# corrs <- corrS(sib_geno , pat_geno)
