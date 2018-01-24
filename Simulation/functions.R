library(psych)
library(cvTools)
library(doParallel)
library(tidyverse)
library(svd)
library(foreach)
library(fanc)
tr <- psych::tr
map <- purrr::map

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
sim_sib_stratify <- function(n, snp_num, int_eff, num_strata, mafs, ps, proportion){
  func <- sib_sim(num_strata,
                  mafs,
                  ps,
                  proportion,
                  simControl = list(N = n, 
                                    num_SNP = snp_num, 
                                    main_effect = 1.5,
                                    r = 0,
                                    interaction_effect = int_eff,
                                    cov_effect = 1, 
                                    level = int_lev,
                                    num_parents = 0, 
                                    num_sib = 2, 
                                    num_main = 0,
                                    num_interact = 1,
                                    model = "logistic",
                                    genetic_model = "additive",
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
                                    population = FALSE))
                                    
}

simPopLE_l2_sp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, m_eff = NULL, mode = NULL){
          func <- simPopLE(n,
                     num_SNP = snp_num,
                     MAF = maf,  
                     main_effect = 1.5, 
                     r = 0,
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

simPopLD_l2_sp <- function(n, snp_num, maf, p, int_eff, int_lev,int_num, ld, m_eff = NULL, mode = NULL){
  func <- simPopLE(n,
                   num_SNP = snp_num,
                   MAF = maf,  
                   main_effect = 1.5, 
                   r = ld,
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

cv.episfa <- function(x, nfolds,nfactor = NULL, type = "data", contrast = NULL, sfa_control = list(),...){
 
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
#  print(paste0("The number of factors is ", nfactor))
  folds <- cvFolds(nrow(x), K=nfolds)
  # cross validation
  r <- foreach(i=1:nfolds, .verbose = TRUE) %do% {
#    print(paste0("cross validate round",i))
    timestamp()
    train <- x[folds$subsets[folds$which != i], ]
    validation <- x[folds$subsets[folds$which == i], ]
    if (is.null(contrast) & type == "data"){
      result <- fanc(train, factors = nfactor, ...)
    } else if (!is.null(contrast) & type == "data"){
        n1 <- nrow(train)
        n2 <- nrow(contrast)
        n3 <- round((n1 - 3)*(n2 - 3)/(n1 + n2 - 6) + 3)
        cov_mat <- partial(train, contrast)
        result <- fanc(covmat = cov_mat, factors = nfactor, n.obs = n3,control = sfa_control, ...)
      } else if (type == "poly") {
      result <- fanc(factors = nfactor, covmat = psych::polychoric(train)$rho, n.obs = nrow(train),...)
    } else {
      stop("type must be data or covmat or poly")
    }
    # extract statistics
    loadings <- reduce(result$loadings, c)
    ui <- plyr::alply(result$uniquenesses,c(1,3),`[`)
#    print("begin to calc kl")
#    kl <- map2_dbl(loadings, ui, kullback, test = validation)
#    print("finished")
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
#    best_kl <- best_gr(kl, nonzero)
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
#    loading_kl <- best_loading(result, best_kl)
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
#    nz_kl <- nonZeroLoad(loading_kl$loadings, coef = FALSE)
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
    return(list(#kl  = matrix(kl, nrow  = 30),
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
#                nz_kl = nz_kl, 
                nz_aic=nz_aic, 
                nz_caic = nz_caic,
                nz_bic = nz_bic, 
                nz_ebic0.5 = nz_ebic0.5,
                nz_ebic0.75 = nz_ebic0.75,
                nz_ebic1 = nz_ebic1,
                nz_hbic0.5 = nz_hbic0.5,
                nz_hbic0.75 = nz_hbic0.75,
                nz_hbic1 = nz_hbic1,
#                loading_kl = loading_kl,
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
  best <- consistency
  snps <- str_extract_all(best,"SNP[0-9]*") %>% unlist() 
  cols <- setdiff(colnames(data), snps)
#  print(snps)
#  print(cols)
  dat <- data[,cols]
  return(dat)
}

episfa <- function(dat, nfolds, recursion = 1, criteria = "ebic",contrast = NULL, sfa_control = list(),...){
  criteria <- str_to_lower(criteria)
  if (criteria  %in% c("ebic","hbic")){
    criteria <- paste0("nz_",criteria,c(1, 0.75, 0.5))
  } else {
    criteria <-  paste0("nz_",criteria)
  }
  inters <- character()
  if (isSymmetric(dat)){
  type <- "covmat"
  } else {
    type <- "data"
  }
  ## recursion
  for(i in 1:recursion) {
#    print(paste0("Searching for epistatic effects: number ",i,", ..."))
    result <- cv.episfa(dat, nfolds, 1, type, contrast, sfa_control, ...)
    consistency <- map(1:length(criteria),~cv.consistency(result, criteria[.]))
#    print(consistency)
    # consistency <- consist
    elements <- lengths(consistency)
    if (sum(elements) >= 1) {
      effect <- consistency[elements >= 1][[1]][1]
      name <- names(effect)
      print(paste0("Found interaction ",name))
      dat <- effRemove(name, dat)
       if (is.null(dat) || ncol(dat) < 2){
         break()
       }
      inters[name] <- effect
    } else {
#      sprintf("no interaction found for criteria %s, algorithm halts", criteria)
      inters[i] <- NA
      break()
    }
  }
  return(inters)
}
    
episfa_sim <- function(n_rep = 100, recursion = 5, cvfolds = 10, ncores = NULL,sim_func = simPopLE_l2_sp, sim_control = list(), criteria = "ebic",ld = FALSE,sfa_control = list(),...){
    if (is.null(ncores)){
      ncores <- detectCores(logical = FALSE)
    }
#    sprintf("running on %i cores",ncores) %>% 
#      print()
    ## make cluster
    cl <- makeCluster(ncores,outfile = "info180124.txt",rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
    registerDoParallel(cl)
    ## hyperparameters
    num_interact <- sim_control[["int_num"]]
    global_funcs <- lsf.str(.GlobalEnv) %>% as.vector()
    ## looping
    benchmark <- foreach(i = 1:n_rep, 
                         .export = global_funcs, 
                         .packages = c("purrr","stringr","fanc","dplyr","sigmoid","truncnorm","data.table","cvTools", "foreach","psych"),
                         .verbose = TRUE) %dopar% {
      #timing
      time_start <- Sys.time()
      # simulate data
      simdata <- do.call(sim_func,sim_control)
      inters <- interSNP(simdata) %>% map_chr(paste0, collapse = "") 
      snp_name <- colnames(simdata) %>% str_subset('^SNP[0-9]+$')
      df_co <- simdata[Y == 1,snp_name, with = FALSE] %>% as.matrix()
      control <- simdata[Y == 0, snp_name, with = FALSE] %>% as.matrix()
      if (ld == TRUE){
        contrast <- control
      } else {
        contrast <- NULL
      }
      # episfa run
      result_episfa <- episfa(df_co, cvfolds, recursion, criteria, contrast, sfa_control,...)
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
    in_inters1 <- str_detect(candidate, paste0(inters,c("$|","[^0-9]"),collapse = "")) %>%
      sum() 
    in_inters2 <- str_detect(inters, paste0(candidate,c("$|","[^0-9]"), collapse = ""))  %>%
      sum() 
    in_inters <- in_inters1 + in_inters2
    #    print(in_inters1)
    #    print(in_inters2)
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

simResults <- function(sim_control, sfa_control, n_rep = 100, recursion = 2, cvfolds = 5,ncores = NULL, sim_func = simPopLE_l2_sp,ld=FALSE, criteria = "ebic", save = TRUE){
  #  print(length(sim_control))
  sim_param <- sim_control %>% as.list
  scene <- episfa_sim(n_rep, 
                      recursion,
                      cvfolds, 
                      ncores,
                      sim_func,
                      sim_param, 
                      criteria,
                      ld,
                      sfa_control)
  name <- pmap_chr(sim_param, paste, sep = '-')
  if (save == TRUE){
    write_rds(scene,paste0("results/intern1p05", name,".rds"))
  }
  return(scene)
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


partial <- function(case, control){
  cormat_case <- cor(case) %>% fisherz()
  cormat_control <- cor(control) %>% fisherz()
  cor_cc <- (cormat_case - cormat_control) %>% fisherz2r() 
  diag(cor_cc) <- 1
  return(cor_cc)
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

corrS <- function(dat, famid = "fid", bw = "w"){
  SNP_names <- str_subset(colnames(dat),"^SNP[0-9]*$")
  df_co <- setDT(dat)[Y == 1, SNP_names, with = FALSE] 
  df_co_fid <- dat[[famid]][dat$Y == 1]
  df_co_strata <- map_dfc(setDT(dat)[,SNP_names, with = FALSE],stratas)
#  dat_surf <- map2_dfc(dat[,SNP_names,with = FALSE], df_co_strata, stratScale, df_co_fid)
  if (bw == "w"){
    dat_surf <- map2_dfc(df_co, df_co_strata, stratScale, df_co_fid)
  } else {
    dat_surf <- map2_dfc(df_co, df_co_strata, stratMean, df_co_fid)
  }
  return(cor(dat_surf))
}

corr <- function(vec1,vec2){
  non_zero <- (vec1*vec2) != 0
  vec1 <- vec1[non_zero]
  vec2 <- vec2[non_zero]
  r <- vec1 %*% vec2 /(sqrt((vec1 %*% vec1)*(vec2 %*% vec2)))
  return(r)
}

stratas <- function(vec){
  num <- length(vec)
  s1geno <- vec[seq(1,num,2)]
  s2geno <- vec[seq(2,num,2)]
  screen <- s1geno != s2geno
  non_info <- which(!screen)
  strata <- numeric(num/2)
  strata[which(screen)] <- s1geno[screen] + s2geno[screen]
  strata[non_info] <- 9
  return(strata)
} 

stratScale <- function(vec, strata, fid){
  index <- match(fid, unique(fid))
  strata_co <- strata[index]
#  strata <- rep(strata,each = 2)
  vec[strata_co == -1] <- vec[strata_co == -1] - mean(vec[strata_co == -1])
  vec[strata_co == 0] <- vec[strata_co == 0] - mean(vec[strata_co == 0])
  vec[strata_co == 1] <- vec[strata_co == 1] - mean(vec[strata_co == 1])
  vec[strata_co == 9] <- 0
  return(vec)
}

stratMean <- function(vec, strata, fid){
  index <- match(fid, unique(fid))
  strata_co <- strata[index]
#  strata <- rep(strata,each = 2)
  vec[strata_co == -1] <- mean(vec[strata_co == -1])
  vec[strata_co == 0] <- mean(vec[strata_co == 0])
  vec[strata_co == 1] <- mean(vec[strata_co == 1])
  return(vec)
}
corRank <- function(cormat){
  order_b <- rank(cormat[upper.tri(cormat)])
  n <- length(order_b)
  rankb <- matrix(0, nrow(cormat), ncol(cormat))
  rankb[upper.tri(rankb)] <- order_b
  rankb[lower.tri(rankb)] <- t(rankb)[lower.tri(rankb)]
  diag(rankb) <- NA
  return(n - rankb + 1)
}

integrateCorr <- function(dat,nsnps,Y = "Y", famid = "fid", weights = NULL){
  if (is.null(weights)){
    weights <-c(1/sqrt(2),1/sqrt(2))
  } else if (length(weights) != 2){
    stop("length of weights should be 2")
  }
  corr_b <- corrS(dat, famid, bw = "b")
  corr_w <- corrS(dat, famid, bw = "w")
  zw <- fisherz(corr_w)*sqrt(sum(dat[[Y]]) - 3)
  zb <- ((corRank(corr_b) - 0.5)/choose(nsnps,2)) %>% qnorm()
  z_int <- (zw*weights[2] - zb*weights[1])*sqrt(2)/sqrt(sum(dat[[Y]]) - 3) 
  corr_int <- fisherz2r(z_int)
  diag(corr_int) <- 1
  return(corr_int)
}

permuteMatrix <- function(mat, margin = 2){
  num_dim <- dim(mat)[margin]
  nums <- sample.int(num_dim,num_dim,replace = FALSE)
  mat <- mat[,nums]
  return(mat)
}



## fam-mdr
as.SNPs <- function(dat, num_snps){
  SNPs <- dat[,1:num_snps] + 1
  SNPs <- as.data.frame(SNPs)
  return(SNPs)
}

as.pedigree <- function(dat){
  pedigree <- dat[,c("fid","sid","mid","faid","sex"), with = FALSE]
  pedigree$sid[pedigree$sid %in% c(1,2)] <- paste0(pedigree$fid[pedigree$sid %in% c(1,2)],pedigree$sid[pedigree$sid %in% c(1,2)])
  pedigree <- as.data.frame(pedigree)
  return(pedigree)
}

as.pheno <- function(dat){
  cov <- paste0("cov",1:10)
  pheno <- dat[,c("Y","sex","age",cov),with = FALSE] %>%
    as.data.frame()
  return(pheno)
}

kinship_sib <- function(dat){
  n <- nrow(dat)/2
  mat <- matrix(c(0.5,0.25,0.25,0.5),2,2)
  kin <- bdiag(map(1:n, function(x) mat)) %>% as.matrix()
  return(kin)
}

formBuild <- function(outcome, num_SNPs, covs = NULL, strata = NULL){
  SNPs <- paste0("SNP",1:num_SNPs) %>% c(covs)
  if (is.null(strata)){
    form <- paste0(outcome,"~", paste0(SNPs,collapse = "+")) %>%
      as.formula()
  } else {
    form <- paste0(outcome,"~", paste0(SNPs,collapse = "+"), "+strata(",strata,")") %>%
      as.formula() 
  }
  return(form)
}

addParents <- function(pedfile){
  pedfile$mid <- paste0(pedfile$mid,"p")
  pedfile$faid <- paste0(pedfile$faid, "p")
  pedfile$sid <- as.character(pedfile$sid)
  mo <- pedfile[,c("fid","mid"),with=FALSE] 
  fa <- pedfile[,c("fid","faid"),with=FALSE]
  colnames(mo) <-c("fid","sid")
  colnames(fa) <-c("fid","sid")
  parents <- bind_rows(mo,fa) %>% distinct()
  pedfile <- pedfile %>% union_all(parents)
  return(pedfile)
}

permuteY <- function(Y){
  Y <- sample(Y, length(Y), replace = FALSE)
  return(Y)
}

fammdr <- function(dat,null_p = NULL, P = 0.1){
  SNP_names <- colnames(dat) %>% str_subset("^SNP[0-9]*$")
  nsnp <- length(SNP_names)
  SNPS <- as.SNPs(dat, nsnp)
  pedigree <- as.pedigree(dat)
  phenotype <- as.pheno(dat)
  SNPs.factor <- map_dfc(SNPS, as.factor)
  genopheno <- bind_cols(SNPs.factor,phenotype)
  # kin <- kinship(pedigree[[2]][1:4],pedigree[[3]][1:4],pedigree[[4]][1:4])
  
  form <- formBuild("Y",nsnp,strata = "fid",rem = TRUE)
  Yfit <- lme4::glmer(form, data = genopheno,family = binomial(link = "logit"),nAGQ = 0,control = glmerControl(calc.derivs = FALSE)) %>%
    summary()
  yres <- Yfit$residuals
  
  for (j in 1:ncol(SNPS)) SNPS[,j] <- as.factor(SNPS[,j])
  EXPO <- list()
  EXPO[["0"]] <- (SNPS==0)
  EXPO[["1"]] <- (SNPS==1)
  EXPO[["2"]] <- (SNPS==2)
  mdr <- MBMDR(yres,EXPO,SNPS,ESTRAT =NULL,PVAL=P,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)
  n_col_mdr <- ncol(mdr)
  m <- which.min(as.numeric(mdr[,n_col_mdr]))
  SNPs <- mdr[m,(n_col_mdr - 4):1]
  inters <- paste0("SNP",SNPs, collapse = "")
  p_min <-  mdr[m,n_col_mdr] %>% as.numeric()
  sprintf("Sample: interaction %s : %f",  inters, p_min) %>% 
    print()
  if (!is.null(null_p)){
    if (p_min <= null_p){
      p_val <- c(p_min) %>% setNames(inters)
      return(p_val)
    } else {
      return(NULL)
    }
  } else {
    return(p_min)
  }
  #  null_p <- numeric(permutation)
  #  dist <- 1
  #  print(p_min)
  ## permutation
  #  for (perm_n in 1:permutation){
  #    yres <- permuteY(yres)
  #    mdr <- MBMDR(yres,EXPO,SNPS,ESTRAT =NULL,PVAL=P,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)
  #   null_p[perm_n] <- mdr[which.min(mdr[,6]),6] %>% as.numeric()
  # print(null_p)
  #    if (null_p[perm_n] < p_min){
  #      dist <- dist + 1
  #    } 
  #   if (dist >= 2){
  #      break
  #   }
  #  }
  #  p_permute <- 1 - sum(null_p >= p_min)/permutation
  #  if (p_permute <= 0.05){
  #    return(c(inters = p_min))
  #  } else {
  #    print("no interaction found")
  #   return(NULL)
  # }
}


fammdr_sim <- function(n_rep = 100,P = 0.5,null_p = NULL, ncores = NULL,sim_func = simPopLE_l2_sp, sim_control = list(), verbose = TRUE){
  if (is.null(ncores)){
    ncores <- detectCores(logical = FALSE)
  }
  #  sprintf("running on %i cores",ncores) %>% 
  #    print()
  ## make cluster
  cl <- makeCluster(ncores,outfile = "log0111.txt",rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
  registerDoParallel(cl)
  ## hyperparameters
  num_interact <- sim_control[["int_num"]]
  
  global_funcs <- lsf.str(.GlobalEnv) %>% as.vector()
  ## looping
  tryCatch({
    benchmark <- foreach(i = 1:n_rep, 
                         .export = global_funcs, 
                         .packages = c("purrr","stringr","GenABEL","dplyr","sigmoid","truncnorm","data.table","cvTools", "foreach","Matrix"),
                         .verbose = verbose) %dopar% {
                           #timing
                           time_start <- Sys.time()
                           # simulate data
                           simdata <- do.call(sim_func,sim_control)
                           inters <- interSNP(simdata) %>% map_chr(paste0, collapse = "") 
                           snp_name <- colnames(simdata) %>% str_subset('^SNP[0-9]+$')
                           # fammdr runsss
                           result_mdr <- fammdr(simdata, null_p, P)
                           if (is.null(null_p)){
                             return(result_mdr)
                           }
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
                           
                           inter_names <-names(na.omit(result_mdr))
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
                                       result = result_mdr,
                                       true_inter = inters,
                                       time = time_diff))
                         }
    if (is.null(null_p)){
      benchmark <- as.numeric(benchmark) %>% sort()
      p95 <- benchmark[round(n_rep * 0.95)]
      stopCluster(cl)
      return(p95)
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
    alpha3 <- fp/(n_rep)
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
  },
  error = function(e){
    print(e)
    stopCluster(cl)
  })
  
}


simResults_mdr <- function(sim_control,null_p = NULL, n_rep = 100,  P = 0.1,ncores = NULL, sim_func = simPopLE_l2_sp, save = TRUE, verbose = TRUE){
  sim_param <- sim_control %>% as.list
  
  scene <- fammdr_sim(n_rep, 
                      P, 
                      null_p,
                      ncores,
                      sim_func,
                      sim_param,
                      verbose
  )
  name <- pmap_chr(sim_param, paste, sep = '-')
  if (save == TRUE){
    write_rds(scene,paste0("results/intern1p05_mdr", name,".rds"))
  }
  return(scene)
}  


# sib_geno <- test1200[!is.na(mid),] %>% arrange(fid)
# sib_geno <- sib_geno[,1:200]
# pat_geno <- test1200[is.na(mid),] %>% arrange(fid)
# pat_geno <- pat_geno[,1:200]

# corrs <- corrS(sib_geno , pat_geno)
