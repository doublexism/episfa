source("D:\\Fangshan\\thesis\\interaction\\Simulation\\simulate data LE.R")
nsnps <- 50
test <- simPopLE(1000, num_SNP = nsnps, 
                     MAF = 0.1,  
                     main_effect = 1.5, 
                     interaction_effect =2, 
                     margin_effect = NULL,
                     cov_effect = 2, 
                     level = 2,
                     num_parents = 0, 
                     num_sib = 2, 
                     num_main = 1,
                     num_interact = 5,
                     model = "logistic",
                     genetic_model = "additive",
                     p = 0.05,
                     scale_weibull = 80,
                     shape_weibull = 4,
                     age_lower = 20,
                     age_higher = 80,
                     sex_effect = 1.2,
                     age_effect = 1.05,
                     age_mean = 50,
                     age_varb = 10,
                     age_varw = 5,
                     num_cov = 10,
                     cov_sigma = NULL)

## validate simulated data
SNPs <- paste0("SNP",1:5, collapse = "+")%>%paste("Y ~",.) %>% as.formula()
results <- glm(Y~SNP60*SNP156, data=test, family = binomial(link = "logit")) 
summary(results)
library(survival)
result <- clogit(Y~SNP20*SNP42+strata(fid), data = test)
summary(result)
lm(SNP20 ~ SNP42, data = test) %>% summary()
lm(SNP20 ~ SNP42, data = test[Y == 1,]) %>% summary()
lm(SNP20 ~ SNP42, data = test[Y == 0,]) %>% summary()

result <- clogit(Y~SNP10*SNP82+strata(fid), data = test)
summary(result)
lm(SNP10 ~ SNP82, data = test) %>% summary()


result <- clogit(Y~SNP97+SNP168+SNP154+SNP97:SNP168:SNP154+strata(fid), data = test)
summary(result)
lm(SNP97 ~ SNP168:SNP154, data = test) %>% summary()

result <- clogit(Y~SNP157+SNP38+SNP62+SNP157:SNP38:SNP62+strata(fid), data = test)
summary(result)

df_co <-test %>% dplyr::filter(Y == 1)
df_ctrl <- test %>% dplyr::filter(Y == 0)
lm(SNP60 ~ SNP156, data = df_co) %>% summary()
lm(SNP60 ~ SNP156, data = test[Y==0,]) %>% summary()

library(fanc)
library(stringr)
SNP_names <- paste0("SNP", 1:nsnps)
Mat_x <- as.matrix(test[,..SNP_names])


dat <- test[,1:nsnps]
cor_dat <- cor(dat)
cor_dat_p <- abs(cor_dat)
fa.parallel(abs(cor_dat) , n.obs= 2000)
cor_case <- cor(df_co[,1:nsnps])
cor_control <- cor(df_ctrl[,1:nsnps])
cor_2 <- cor_dat - cor_control + diag(rep(1,nsnps))

result <- factanal(factors = 5, rotation = "varimax", covmat = cor_2, n.obs = 2000)
fa.parallel(cor_2, n.obs = 2000)

result1 <- fanc(factors = 10, covmat = cor_dat, n.obs = nrow(dat),control = list(openmp = TRUE, num.threads = 8))
result2 <- fanc(factors = 10, covmat = cor_2, n.obs = nrow(dat), control = list(openmp = TRUE, num.threads = 8))

nfacs <- cv_evaluate(10, qBest,dat,ncores = 3, 1:37)

cors <- foreach(i = 1:100) %do% {
  a <- sample_frac(test, 0.7)
  correlation <- cor(a[Y == 1, 1:nsnps]) - cor(a[Y == 0, 1:nsnps]) + diag(rep(1,50))
  return(correlation)
}
cor_bs <- reduce(cors, `+`)/100


