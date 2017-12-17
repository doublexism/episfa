source("D:\\Fangshan\\thesis\\interaction\\Simulation\\simulate data LE.R")


test <- simPopLE(1000, num_SNP = 100, 
                     MAF = 0.2,  
                     main_effect = 1.5, 
                     interaction_effect =2, 
                     margin_effect = NULL,
                     cov_effect = 1.5, 
                     level = 2, 
                     num_parents = 0, 
                     num_sib = 2, 
                     num_main = 5,
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

result <- clogit(Y~SNP60*SNP156+strata(fid), data = test)
summary(result)
lm(SNP60 ~ SNP156, data = test) %>% summary()

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
SNP_names <- paste0("SNP", 1:100)
Mat_x <- as.matrix(test[,..SNP_names])

qBest <- function(nfactor, data){
  result <- fa(data, nfactors = nfactor, rotate = "varimax")
  beta <- result$loadings
  ui <- result$uniquenesses
  kl <- kl(beta, ui, data)
  return(kl)
}
timestamp()
fa_num <-  cv_evaluate(10, qBest, test, ncores =4, 1:100)
timestamp()


rbenchmark::benchmark(
  result1 <- fanc(Mat_x,factors = 10, control = list(openmp = TRUE, num_threads = 4)),
  result2 <- fanc(Mat_x,factors = 10, control = list(openmp = TRUE, num_threads = 8)),
  result3 <- fanc(Mat_x,factors = 10),
  replications = 1
)
lm(SNP9~SNP49, data = test) %>% summary()
lm(SNP9~SNP49, data = df_co) %>% summary()
lm(SNP9~SNP49, data = df_ctrl) %>% summary()

dat <- test[,1:100]
cor_dat <- cor(dat)
cor_case <- cor(df_co[,1:100])
cor_control <- cor(df_ctrl[,1:100])
cor_2 <- cor_dat - cor_control + diag(rep(1,100))

result <- factanal(factors = 5, rotation = "varimax", covmat = cor_2, n.obs = 2000)
fa.parallel(cor_2, n.obs = 2000)

result1 <- fanc(factors = 10, covmat = cor_dat, n.obs = nrow(dat),control = list(openmp = TRUE, num.threads = 8))
result2 <- fanc(factors = 10, covmat = cor_2, n.obs = nrow(dat), control = list(openmp = TRUE, num.threads = 8))


nfacs <- cv_evaluate(10, qBest,dat,ncores = 3, 1:37)
