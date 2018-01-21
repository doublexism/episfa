test1 <- simPopLE_l2_sp(1000, 50, 0.5, 0.05, 2, 2, 1)
test2 <- simPopLD_l2_sp(1000, 50, 0.5, 0.05, 1, 2, 1,0.5)

test2 <- simPopLE_l2_sp(1000, 100, 0.2, 0.05, 2, 3, 1)
test3 <- simPopLE_l2_sp(1000, 100, 0.2, 0.05, 2, 4, 1)
test4 <- simPopLE_l2_sp(1000, 100, 0.2, 0.05, 2, 5, 1)
df_co1 <- as.matrix(test1[Y == 1, 1:50])
control1 <- as.matrix(test1[Y == 0, 1:50])

sidf_co2 <- as.matrix(test2[Y == 1, 1:100])
df_co2log <- log_trans(df_co2, offset = 2)


df_co3<- as.matrix(test3[Y == 1, 1:100])
df_co4 <- as.matrix(test4[Y == 1, 1:100])

result1 <- episfa(df_co1, 5, nfactor= 3, type = "data")

result2 <- episfa(df_co2, 10,nfactor= 3, type = "data")
result2log <- episfa(df_co2log, 10,nfactor= 3, type = "data")


result3 <- episfa(df_co3, 10,type = "data")
result4 <- episfa(df_co4, 10,type = "data")


result1.simple <- fanc(df_co1,factors=3)
#result1.prenet <- result3 <- fanc(df_co1,factors=2,  type = "prenet")
result2.simple <- fanc(as.matrix(test[Y==1, 1:100]),factors=5, control = list(Delta = 0.01))

SNP_names <- str_subset(colnames(test), "^SNP[0-9]*$")
val <- simVal(test1)
clogit(Y ~ SNP18*SNP67 + strata(fid), data = test1) %>% summary()


df_co <- as.matrix(test[Y == 1, 1:20])
df_ctrl <- as.matrix(test[Y == 0, 1:20])
co_df <- cor(df_co)
co_ctrl <- cor(df_ctrl)
co_dif <- co_df - co_ctrl + diag(rep(1, ncol(df_co)))

cov_df_cltr <- cor(test[Y==0, 1:20])
df_co_compress <- svd_compress(df_co,18)
cov_df_co <- cor(df_co_compress)

result1 <- episfa(df_co, 10,type = "covmat", control = df_ctrl)
result2 <- episfa(df_co, 10,type = "data")


fa.parallel(df_co)
fa.parallel(df_co_compress)
factanal(df_co_compress, factors = 3, rotation = "varimax")
factanal(df_co, factors = 3, rotation = "varimax")
factanal(df_co_s, factors = 3, rotation = "varimax")

d <- svd(df_co_s)
v <- as.vector(d$d[19:20]**2 %*% d$v[19:20,]**2 /1044)

noise <- MASS::mvrnorm(nrow(df_co_s),rep(0, 20), Sigma = diag(v))


lm(SNP17 ~ SNP64 , data = test1[Y==1]) %>% summary()
lm(SNP20 ~ SNP13 , data = test[Y==0]) %>% summary()

lm(SNP20 ~ SNP13 * Y, data = test) %>% summary()

map(result1,`[`, 21)

kl_2_2 <- matrix(result1[[2]]$kl,nrow = 30)

result <- fanc(factors=2, covmat = co_dat, n.obs = 2000)

test1 <- sample_frac(test, 1, replace = TRUE)

lm(SNP3~SNP4, data = test[Y==1]) %>% summary()
clogit(Y ~ SNP18 + SNP37 + SNP18:SNP37 +strata(fid), data = test1) %>% summary()
clogit(Y ~ SNP1 + SNP2 + SNP1:SNP2 +strata(fid), data = test1) %>% summary()

small_control <- control1[,1:4]

pcor_control <- cor2pcor(cor(small_control))

lm(SNP1 ~ SNP1 + SNP2 + SNP3 + SNP4, data = test1[Y == 0,1:4])
cor_case1 <- cor(df_co1)
cor_control1 <- cor(control1)
cor_cc <- diag(ncol(df_co1))+cor_case1 - cor_control1

result1 <- fanc(df_co1, factors = 1, control = (eta = 0.025))

result2 <- fanc(covmat = cor_cc, n.obs = 1000, factors = 1, control = (eta = 0.025))


cor_1 <- cor(SNP1[strata1 != 9 & strata2 != 9 & test1$Y == 1], SNP2[strata1 != 9 & strata2 != 9 & test1$Y == 1])

test1_test <- test1_sort[c("fid","SNP1","Y","SNP2")]

test1_sort <- test1 %>% arrange(fid,desc(Y))
n_y <- test1[,.(n  = sum(Y)),fid]$n
test1_co <- test1_sort[seq(1,nrow(test1_sort),2), 1:50]
test1_control <- test1_sort[seq(2,nrow(test1_sort),2),1:50]
test1_diff <- (test1_co + test1_control)/2
cor_diff <- cor(test1_diff[n_y == 1,])
constant <- rep(1, length(n_y[n_y ==1]))
lm(SNP2 ~ SNP1, data=test1_diff[n_y ==1,]) %>% summary()
glm(constant ~ 0 + SNP1*SNP2, data = test1_diff[n_y == 1,],family = binomial(link = "logit")) %>% summary()

clogit(Y ~ SNP1*SNP2 + strata(fid), data = test1[rep(n_y == 1, each = 2)]) %>% summary()


cor(SNP1[test1$Y == 1], SNP2[ test1$Y == 1])

order_b <- rank(corr_b[upper.tri(corr_b)])
rankb <- matrix(0, nrow(corr_b), ncol(corr_b))
rankb[upper.tri(rankb)] <- order_b
rankb[lower.tri(rankb)] <- t(rankb)[lower.tri(rankb)]

corr_int <- integrateCorr(test1, nsnps = 50, weight = c(1,0))
r <- fanc(covmat = corr_int, factors = 1, n.obs = 1028)


mat <- MASS::mvrnorm(n = 1000, c(0,0),Sigma = matrix(c(1,0.3,0.3,1),2))

#cutoff <- c(0.81,0.99) %>% qnorm()
cutoff <- c(0.25,0.75) %>% qnorm()

index1 <- mat > cutoff[2] 
index2 <- mat < cutoff[2] & mat > cutoff[1]
index3 <- mat < cutoff[1]
mat[index1] <- 1
mat[index2] <- 0
mat[index3] <- -1
cor(mat)

test3 <- sib_sim(num_strata  = 1,
                 MAF_strata = 0.1,
                 P_strata = 0.05,
                 proportion_strata = 1,
                 N = 1000,num_SNP = 50, MAF = 0.1,main_effect = 1.5, interaction_effect = 2, p = 0.05)


clogit(Y ~ SNP1*SNP2+strata(fid), data = test3) %>% summary()


test3 <- sib_sim(2,c(0.5, 0.25),c(0.05,0.05) ,c(0.5,0.5),
                 list(N = 1000,num_SNP = 50, MAF = 0.1,main_effect = 1.5, interaction_effect = 2, p = 0.05))

cor_case_test3 <- cor(test3[Y == 1, 1:50]) 
cor_control_test3 <- cor(test3[Y == 0, 1:50]) 

cor_test3<- partial(test3[Y == 1, 1:50], test3[Y == 0, 1:50])

r <- fanc(covmat = cor_test3, n.obs = 500,factors = 1)
