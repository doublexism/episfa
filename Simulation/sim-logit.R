library(sigmoid)
library(tidyverse)
library(foreach)
sim <- function(prevalence, n, n_var, beta){
  maf <-  0.1
  probs <- c((1-maf)**2, 2*maf*(1-maf), maf**2)
  x <- sample(c(-1,0,1), n_var*n, prob = probs, replace = TRUE) %>%
    matrix(ncol = n_var)
  x <- cbind(rep(1, n),x, x[,1]*x[,2], x[,3]*x[,4], x[,5]*x[,6])
  beta <- c(logit(prevalence), beta)
  offset <- sum(beta[2:(n_var+1)]*((2*maf - 1))) + beta[n_var+2]*(1-2*maf)**2 + beta[n_var+3]*(1-2*maf)**2+ beta[n_var+4]*(2*maf - 1)**3
  beta[1] <- beta[1]-offset
  p <- sigmoid(x %*% beta) 
  y <- runif(n) < p
  df <- cbind(x,y) %>% as.data.frame()
  colnames(df) <- c("int",paste0("x",1:n_var),"x12","x23","x123","y")
  return(df)
}

df <- sim(0.05, 50000, 100, c(rep(0,100),1,1,1))
df_co <- df[df$y == 1,]

df_co_mat <- as.matrix(df_co[2:101])

cor_df_co <- cor(df_co_mat)

df_co_mat_s <- plyr::aaply(df_co_mat,2, function(v) (v-mean(v))/sd(v)) %>% t()


r <- svd(df_co_mat_s)
sig <- sum(r$d[91:100]**2/nrow(df_co_mat_s))/100
guassian_noise <- rnorm(length(df_co_mat_s),0,sig) %>% matrix(ncol = ncol(df_co_mat_s))

depth <- 90
df_co_mat_sr <- r$u[,1:90] %*% diag(r$d[1:90]) %*% t(r$v[,1:90])
df_co_mat_sr <- df_co_mat_sr + guassian_noise

fanc <- cv_evaluate(nfolds = 10, func = qBest, data = df_co_mat_sr, ncores = 3, 1:20)
fa.parallel(df_co_mat_sr)

results <- factanal(df_co_mat_sr, factors = 20, rotation = "varimax")

cor_df_co1 <- cor(df_co_mat_s)
cor_df_co <- cor(df_co_mat_sr)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()



psych::fa.parallel(df_co_mat_s)
psych::fa.parallel(df_co_mat_sr)
  
result <- princomp(df_co_mat_sr)
res <- eigen(cor(df_co_mat_sr))


result$loadings %>% dim()
eigens <- result$scores %*% result$loadings

result$loadings[1,]


df <- df_co %>% bind_rows(sample_n(df[df$y == 0,], 2500))

lm(x3 ~ x4 + x5, data=df_co) %>% summary()
lm(x3 ~ x5, data=df_co) %>% summary()
glm(y ~ x1*x2 + x3*x4, data = df,family = "binomial")%>% summary()

result <- princomp(df_co[2:11])

lm(x2 ~ x3, data = df[which(df$y == 0 & df$x1 ==1),]) %>% summary()
lm(x2 ~ x3*y, data = df[which(df$x1 ==1),]) %>% summary()

glm(y~ x1+x2+x3+x4+x1:x2, data = df, family = "binomial") %>% summary()
glm(y~ x1:x2, data = df, family = "binomial") %>% summary()


lm(x1 ~ x2:x3, data = df)%>% summary()
lm(x1 ~ x2:x3, data = df[y ==1,]) %>% summary()
lm(x1 ~ x2:x3, data = df[y ==0,]) %>% summary()

lm(x1 ~ x2 , data=df)  %>% summary()
lm(x1 ~ x2 , data=df[df$y == 0,])  %>% summary()
lm(x1 ~ x2 , data=df[df$y == 1,])  %>% summary()

lm(x1 ~ x3 , data=df)  %>% summary()

lm(x1 ~ x2 + x3)

lm(x1 ~ x2:x3, data = df) %>% summary()

lm(x1 ~ x2:x3, data = df) %>% summary()

lm(x1 ~ x2*x3, data = df[y ==1,]) %>% summary()
lm(x1 ~ x2:x3, data = df[y ==0,]) %>% summary()

lm(x2 ~ x3, data = df[which(df$y == 1 & df$x1 ==-1),]) %>% summary()

## cross-validation
dfSplit <- function(data, n){
  ntest <- round(nrow(data)/n,0)
  ntimes <- rep(ntest, n-1) %>% c(nrow(data) - (n-1)*ntest)
  divide_n <- rep(1:n, times = ntimes) 
  matlist <- split(data, divide_n)
  return(matlist)
}


ntest <- round(0.2*nrow(df_co_mat), 0)

df_co_t1 <- df_co_mat[1:ntrain,]
df_co_v1 <- df_co_mat[(ntrain+1):nrow(df_co_mat),]

nfactors <- fa(df_co_v1, nfactors = 1)

result <- fanc(df_co_t1 ,factors = 2, control = list(length.gamma = 6,max.gamma = 10,max.rho=0.05, length.rho = 30))

## find best number of factors
qBest <- function(nfactor, train, val){
  result <- fa(train, nfactors = nfactor, rotate = "varimax")
  beta <- result$loadings
  ui <- result$uniquenesses
  kl <- kl(beta, ui, val)
  return(kl)
}



nfactors <- map_dbl(1:9, qBest, train= df_co_mat_v, val = df_co_mat_t)

kl <- function(beta,ui, test){
  calculated <- beta%*%t(beta) + diag(ui) + diag()
  calculated <- as.matrix(calculated)
  cov_f <- t(test) %*% test /nrow(test) 
  KL <- log(det(calculated)) + tr(solve(calculated)%*%cov_f)
  return(KL)
}

kl_cv <- function(result, test){
  loadings <- unlist(result$loadings)
  ui <- result$uniquenesses %>% plyr::alply(c(1,3))
  KL_mat <- map2_dbl(loadings, ui, kl, test = test) %>% matrix(ncol = 6)
  return(KL_mat)
}

KL <- kl_cv(result, df_co_v1)


cv <- pfa_cv(9, 20, df_co, factors = 2, control = list(length.rho = 20, max.gamma = 10))
cv2 <- pfa_cv(8, 20, df_co, factors = 2, control = list(length.rho = 20, max.gamma = 10))


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
  } 
  stopCluster(cl)
  return(r)
}

depth = 9
s <- svd(df_co_mat)
df_co_mat_compressed <- s$u[,1:depth] %*% diag(s$d)[1:depth,1:depth] %*% t(s$v[,1:depth])
df_co_mat_compressed <- df_co_mat_compressed + MASS::mvrnorm(8849, rep(0,10), diag(rep(0.1,10)))

kl_cv <- cv_evaluate(10, qBest, df_co_mat_compressed, ncores =2, 1:9)
m <- colMeans(kl_cv)

r <- pca(df_co_mat, nfactors = 9)
df_co_mat_compressed <- r$scores %*% t(r$loadings) 
df_co_mat_compressed <- df_co_mat_compressed + MASS::mvrnorm(8849, rep(0,10), diag(rep(1,10)))
