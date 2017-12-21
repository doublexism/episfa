library(sigmoid)
library(truncnorm)
library(foreach)
library(tidyverse)
library(data.table)

# simulate a population under linkage equilibrium
simPopLE <- function(N, 
                     num_SNP, 
                     MAF,  
                     main_effect, 
                     interaction_effect, 
                     margin_effect = NULL,
                     cov_effect = 1.2, 
                     level = 2, 
                     num_parents = 0, 
                     num_sib = 2, 
                     num_main = 0,
                     num_interact = 1,
                     model = "logistic",
                     genetic_model = "additive",
                     p = 0.05,
                     scale_weibull = 80,
                     shape_weibull = 4,
                     age_lower = 20,
                     age_higher = 80,
                     sex_effect = 1.5,
                     age_effect = 1.05,
                     age_mean = 50,
                     age_varb = 10,
                     age_varw = 5,
                     num_cov = 10,
                     cov_sigma = NULL,
                     population = TRUE
                     ) {
  
  # N: number of affected familes
  # num_SNP: number of SNPs
  # MAF: minor allele frequency
  # margin_effect: marginal effect of the interacted SNPs, OR
  # interaction: the effect size of interaction, OR
  # level: level of interaction
  # num_parents: number of available parents in the dataset
  # num_sib: number of siblings in each family
  # num_main: number of marginal effects, default to levels
  # num_interact: number of interaction effects
  # model: underlying disease model, either logistic or the cox model
  # p: the population disease prevalence in logistic model
  # scale_weibull: scale parameter of the  weibull model
  # shape_weibull: shape parameter of the  weibull model
  # age_lower: the lower limit of age
  # gender_effect: effect size of male gender in both models
  # age_effect: effect size of every 1 year increase in age in logistic model
  # seed: random seed;
  
  if(length(level)>1 && length(interaction_effect) == 1){
    interaction_effect <- rep(interaction_effect, length(level))
  }
  if (length(level)>1 && length(num_interact) == 1){
    num_interact <- rep(num_interact, length(level))
  }
  ## if number of margin or interaction effect exceeded number of SNPs, stop
    if (num_main > num_SNP | sum(num_interact*level) > num_SNP | sum(level) > num_SNP){
      stop("more SNP needed to specify this model")
    }
    if ((length(level) %% length(interaction_effect)) != 0 ){
      stop("Incorrect number of interaction effects levels was specified.")
    }
  if(length(num_interact) > length(level)){
    stop("Incorrect number of interaction effects were specified")
  }
  if(is.null(cov_sigma)){
    cov_sigma <- diag(rep(1, num_cov))
  }
  ## covariate effects
  cov_beta <- c(log(cov_effect), -log(cov_effect)) %>% rep(each=num_cov/2)
  ## set seed
  # set.seed(seed)
  ## forward simulation probabilites
  geno_mat <- data.frame(cat = 1:9,
                     a = c(1,0.5,0,0.5,0.25,0,0,0,0),
                     b = c(0,0.5,1,0.5,0.5,0.5,1,0.5,0),
                     c = c(0,0,0,0,0.25,0.5,0,0.5,1))
  ## generate sibling genotypes
  childGenoProb <- function(mother, father, gmat, num_sib){
    cat = 3*mother + father + 5
    break1 <- gmat$a[cat] %>% rep(each = num_sib)
    break2 <- (gmat$b[cat]+gmat$a[cat]) %>% rep(each = num_sib)
    p <- runif(2*length(cat))
    p[p<= break1] <- -1
    p[p > break1 & p<= break2] <- 0
    p[p > break2] <- 1
    return(p)
  }
  
  ## generate interaction terms
  interTerm <- function(vars, dat){
    inter <- dat[vars] %>% reduce(`*`) %>% as.data.frame()
    return(inter)
  }
  
  ## get the list of interaction level for each SNPs
  interLevel <- function(snp, snp_list){
    levels <- lengths(snp_list)[map_lgl(snp_list, ~snp %in% .)]
    return(levels)
  }
 
  ## calculate offsets for marginal effects
  offsetEffect <- function(snp_levels, level, effect, maf = MAF){
    snp_effects <- effect[match(snp_levels, level)] %>% log()
    offset <- map2_dbl(snp_levels, snp_effects, function(x,y) y*(2*maf -1)**(x-1)) %>% 
      sum()
    return(offset)
  }
  ## calculate variance for interaction effects
  var_interact <- function(n,maf){
    pn1 <- (1-maf)**2
    pp1 <- maf**2
    ex <- 2*maf - 1
    var <- (pn1+pp1)**n - ex**(2*n)
    return(var)
  }
  ## outcome calculation under different assumptions
  # logistic
  logistic_func <- function(X, beta){
    P <- sigmoid(X %*% beta)
    print(mean(P))
    y <- (runif(nrow(X)) < P) %>% as.numeric()
    return(y)
  }
  # age-on-set
  cox_func <- function(F0, X, beta){
    P <- F0*exp(X %*% beta)
    P[P>1] <- 1
    y <- (runif(nrow(X)) < P) %>% as.numeric()
    return(y)
  }
  # sampling
  sampling <- function(dat, fids, N){
    fids <- fids%>% sample(N)
    sample <- subset(dat, fid %in% fids)
    return(sample)
  }
  # weibull cumulative risk function
  weibull_pf <- function(x) pweibull(x, shape_weibull, scale_weibull)
  
  ## calculate basic parameters
  ## sample size calculation
  if (model == "cox") {
    p <-  weibull_pf(age_mean-age_lower)
  }
  
  p_affected <- 1 - (1-p)**num_sib
  if (population == TRUE){
    FG_num = N / p_affected * 100 * 2 
  }
  FG_num <- (qnbinom(p = 0.999999, size = N, prob = p_affected) + N)*2
  ## genotype frequencies
  MA_freq <- c((1-MAF)**2, 2*MAF*(1-MAF), MAF**2)

  ## parent genotype
  genotype_parents <- map(rep(FG_num, num_SNP), ~sample(x=c(-1,0,1),size = ., prob = MA_freq, replace=TRUE))
  SNP_names <- paste0("SNP", 1:num_SNP)
  names(genotype_parents) <- SNP_names
  genotype_parents <- as.data.frame(genotype_parents)
 
  ## parent dataframe
  df_parent <- data.frame(id = 1:FG_num, sex = rep(0:1, each = FG_num/2)) %>% bind_cols(genotype_parents)
  df_f <- df_parent %>% filter(sex == 1) %>% sample_n(FG_num/2) %>% mutate(fid = 1:(FG_num/2))
  df_m <- df_parent %>% filter(sex == 0) %>% sample_n(FG_num/2) %>% mutate(fid = 1:(FG_num/2))
  
  ## generate sibling data
  df_sib <- map2_dfc(df_m[SNP_names], df_f[SNP_names], childGenoProb, gmat = geno_mat, num_sib = num_sib)
  df_sib$mid <- rep(df_m$id, each = num_sib)
  df_sib$faid <- rep(df_f$id, each = num_sib)
  df_sib$fid <- rep(df_m$fid, each = num_sib)
  df_sib$sid <- rep(1:num_sib, FG_num/2)
  df_sib$sex <- rbinom(nrow(df_sib), 1, 0.5)
  
  if(num_parents == 1){
    df_sib <- df_sib %>% bind_rows(df_m)
  } else if(num_parents == 2){
    df_sib <- df_sib %>% bind_rows(df_m, df_f)
  }
  
  ## add genetic models
  if (genetic_model == "recessive"){
    df_sib[SNP_names] <-  map(df_sib[SNP_names], function(x)  as.numeric(x > 0))
  } else if (genetic_model == "dorminent") {
    df_sib[SNP_names] <-  map(df_sib[SNP_names], function(x)  as.numeric(x > -1))
  }
  
  ## generate age through truncated normal distribution
  m_age <- rtruncnorm(FG_num/2, mean = age_mean,sd = sqrt(age_varb), a = age_lower, b = age_higher)
  df_sib$age <- map(m_age, 
                    ~rtruncnorm(num_sib, ., sd = sqrt(age_varw), a = age_lower, b = age_higher)) %>%
    unlist()
  df_sib$t_age <- df_sib$age - age_mean
  df_sib$intercept <-  1
  
  ## generate interaction effects
  if (sum(num_interact) != 0){
    # select interacting SNPs
    SNPs <- map(rep(level, num_interact),~sample(SNP_names, size = .))
    SNPs_margin <- unique(unlist(SNPs))
    SNP_idx <- match(SNPs_margin, SNP_names)
    # default effect size for the main effect 
    main <- log(main_effect)
    # generate marginal effects
    if(!is.null(margin_effect)){
      snp_levels <- map(SNPs_margin, interLevel, snp_list = SNPs)
      offset_effect <- map_dbl(snp_levels, offsetEffect, level = level,effect = interaction_effect)
      main <- log(margin_effect) - offset_effect
    }
    # add interaction term to data
    inter_names <- map_chr(SNPs, paste0, collapse = "")
    inter_df <-  map_dfc(SNPs, interTerm, df_sib) %>% setNames(inter_names) 
    df_sib <- inter_df %>% bind_cols(df_sib,.)
  } else {
    inter_names <- NULL
    SNP_idx <- NULL
    SNPs_margin <- NULL
    SNPs <- NULL
  }
  if (num_main > 0){
    SNP_main <- setdiff(1:num_SNP, SNP_idx) %>% `[`(1:num_main)
  } else {
    SNP_main <- NULL
  } 
  
  ## generate covariates
  cov_names <- paste0("cov",1:num_cov)
  covariates <- MASS::mvrnorm(n = nrow(df_sib), rep(0, num_cov), Sigma = cov_sigma) %>% 
    as.data.frame() %>%
    setNames(cov_names) 
  df_sib <- bind_cols(df_sib, covariates)
  
  ## generate outcomes
  if (model == "logistic"){
  X <- as.matrix(df_sib[c("intercept", SNP_names, "sex", "t_age", inter_names, cov_names)])
   beta <- c(logit(p), 
            rep(0, num_SNP), 
            log(sex_effect), 
            log(age_effect), 
            rep(log(interaction_effect), num_interact),
            cov_beta
            )

  if (!is.null(SNP_idx)){
    beta[SNP_idx + 1] <- main
  }
  beta[SNP_main + 1] <- log(main_effect)
  # calculate variance of linear combination
  # var_int <- map_dbl(rep(level, num_interact), var_interact, maf=MAF)
  # print(var_int)
  variance <- c(0, rep(0, num_SNP), 0, age_varb+age_varw ,rep(0, num_interact), diag(cov_sigma))
  #total variance and multiplier
  beta_sq <- beta**2
  total_variance <- sum(variance %*% beta_sq)
  lambda_sq <- pi/8
  multiplier <- sqrt(1+lambda_sq*total_variance)
  
  #calculate offsets
  if(is.null(margin_effect)){
    offset <- (num_main+length(SNPs_margin))*log(main_effect)*(2*MAF - 1) + 
      log(sex_effect)*0.5 +
      sum(num_interact*log(interaction_effect)*(2*MAF - 1)**level)
  } else {
    offset <- (num_main)*log(main_effect)*(2*MAF - 1) +  
      sum(main)*(2*MAF - 1) +
      log(sex_effect)*0.5 +
      sum(num_interact*log(interaction_effect)*(2*MAF - 1)**level)
  }
  # update beta
  beta[1] <- (logit(p))*multiplier - offset
  df_sib$Y <- logistic_func(X, beta)
  
  } else if (model == "cox"){
    F0 <- weibull_pf(df_sib$age - age_lower)
    X <- as.matrix(df_sib[c(SNP_names, "sex", inter_names, cov_names)])
    beta <- c(rep(log(main_effect), num_main),
              rep(0, num_SNP - num_main), 
              log(sex_effect), 
              rep(log(interaction_effect), num_interact),
              cov_beta)
    df_sib$Y <- cox_func(F0, X, beta)
  }
  ## subsampling
  fids <- setDT(df_sib)[!is.na(mid)][, .(n=sum(Y)), by = "fid"][n > 0]$fid
  # check if number of patients met, other wise simutate again the population again
  if (length(fids) < N){
    return(simPopLE(N,
                    num_SNP, 
                    MAF,  
                    main_effect, 
                    interaction_effect, 
                    margin_effect = NULL,
                    cov_effect,
                    level, 
                    num_parents,
                    num_sib, 
                    num_main,
                    num_interact,
                    model,
                    genetic_model,
                    p,
                    scale_weibull,
                    shape_weibull,
                    age_lower,
                    age_higher,
                    sex_effect,
                    age_effect,
                    age_mean,
                    age_varb,
                    age_varw,
                    num_cov,
                    cov_sigma,
                    population
                    ))
  }
  if (population == FALSE) {
    df_sib <- sampling(df_sib, fids, N)
  } else {
    attr(df_sib, "cases") <- fids
  }
  #add interaction snps to dataframe
  attr(df_sib,"interact") <- SNPs
  return(df_sib)
}

sib_sim <- function(n, name, ...){
  times(n) %do% 
    simPopLE(...) %>%
    write_rds(paste0("data/",name,".rds")) 
    return(paste0(name, " completed!"))
    gc()
}

