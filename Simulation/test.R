test1 <- simPopLE_l2_sp(1000, 100, 0.1, 0.05, 2, 2, 1)
test2 <- simPopLE_l2_sp(1000, 100, 0.1, 0.05, 2, 3, 1)
test3 <- simPopLE_l2_sp(1000, 100, 0.1, 0.05, 2, 4, 1)
test4 <- simPopLE_l2_sp(1000, 100, 0.1, 0.05, 2, 5, 1)
df_co1 <- as.matrix(test1[Y == 1, 1:100])
df_co2 <- as.matrix(test2[Y == 1, 1:100])
df_co3<- as.matrix(test3[Y == 1, 1:100])
df_co4 <- as.matrix(test4[Y == 1, 1:100])
result1 <- episfa(df_co1, 10,type = "data")
result2 <- episfa(df_co2, 10,type = "data")
result3 <- episfa(df_co3, 10,type = "data")
result4 <- episfa(df_co4, 10,type = "data")

SNP_names <- str_subset(colnames(test), "^SNP[0-9]*$")
val <- simVal(test)



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


lm(SNP20 ~ SNP13+ cov1+sex+age , data = test[Y==1]) %>% summary()
lm(SNP20 ~ SNP13 , data = test[Y==0]) %>% summary()

lm(SNP20 ~ SNP13 * Y, data = test) %>% summary()

map(result2, `[`, 6)

kl_2_2 <- matrix(result1[[2]]$kl,nrow = 30)

result <- fanc(factors=2, covmat = co_dat, n.obs = 2000)

test1 <- sample_frac(test, 1, replace = TRUE)

clogit( Y ~ SNP3 + SNP6  +SNP3:SNP6+ strata(fid), data = test) %>% summary()