test1 <- simPopLE_l2_sp(1000, 50, 0.1, 0.05,2.5, 2, 1)

test2 <- simPopLE_l2_sp(1000, 50, 0.2, 0.05, 2, 3, 1)
test3 <- simPopLE_l2_sp(1000, 100, 0.2, 0.05, 2, 4, 1)
test4 <- simPopLE_l2_sp(1000, 100, 0.2, 0.05, 2, 5, 1)
df_co1 <- as.matrix(test1[Y == 1, 1:50])

df_co2 <- as.matrix(test2[Y == 1, 1:100])
df_co2log <- log_trans(df_co2, offset = 2)


df_co3<- as.matrix(test3[Y == 1, 1:100])
df_co4 <- as.matrix(test4[Y == 1, 1:100])

result1 <- episfa(df_co1, 5,recursion = 2)

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
clogit(Y ~ SN + SNP70 + SNP85+SNP20:SNP70:SNP85 +strata(fid), data = test1) %>% summary()


##  
SNPS <- as.SNPs(test1, 50)
pedigree <- as.pedigree(test1)
phenotype <- as.pheno(test1)
SNPs.factor <- map_dfc(SNPS, as.factor)
genopheno <- bind_cols(SNPs.factor,phenotype)
# kin <- kinship(pedigree[[2]][1:4],pedigree[[3]][1:4],pedigree[[4]][1:4])
kin <- kinship_sib(test1)

rownames(kin) <- 1:nrow(genopheno)
colnames(kin) <- 1:nrow(genopheno)
form <- formBuild("Y",50,strata = NULL)

timestamp()
Yfit = GenABEL::polygenic(form,kin,genopheno,trait.type="binomial", quiet = TRUE)
timestamp()

forms <- formBuild("Y",50,strata = "fid", rem=TRUE)
Yfit_logit <- clogit(forms, data = genopheno)
timestamp()
Yfit_lmer <- lme4::glmer(forms, data = genopheno,family = binomial(link = "logit"),nAGQ = 0,control = glmerControl(calc.derivs = FALSE))
timestamp()
Yfit_lmer_summmary <- summary(Yfit_lmer)

yres <- Yfit$pgresidualY
yres_lmer <- Yfit_lmer_summmary$residuals

for (j in 1:ncol(SNPS)) SNPS[,j] <- as.factor(SNPS[,j])
EXPO <- list()
EXPO[["0"]] <- (SNPS==0)
EXPO[["1"]] <- (SNPS==1)
EXPO[["2"]] <- (SNPS==2)


#SNPS <- SNPS[1:10]
rbenchmark::benchmark(
  mdr <- MBMDR_o(yres,EXPO,SNPS,ESTRAT =NULL,PVAL=0.1,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)
  ,replications = 1)

mdr[] <- map(mdr, as.character)
mdr$`Pr(>|t|)` <- as.numeric(mdr$`Pr(>|t|)`)
index <- mdr[which.min(mdr$`Pr(>|t|)`),][c("V1","V2")] %>% 
  as.numeric() %>% sort()
selected.SNPs <- SNPS[index]




pair <- t(combn(1:50,2))

clogit(Y ~ SNP1 * SNP2 +strata(fid), data = test1) %>% summary()

lm(yres_lmer ~ SNP1 + SNP2, data = genopheno) %>% summary()
# permutation

mdr1 <- MBMDR(yres,EXPO,SNPS,ESTRAT =NULL,PVAL=0.1,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)
mdr2.2 <- MBMDR(yres_lmer,EXPO,SNPS,ESTRAT =NULL,PVAL=0.1,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)

ymin <- min(as.numeric(mdr2.2[,6]))
i <- 0

ps <- numeric()
ts <- numeric()
repeat{
i <- i+1
yres_p <- permuteY(yres_lmer)
mdr2 <- MBMDR(yres_p,EXPO,SNPS,ESTRAT =NULL,PVAL=0.01,dimen=2,first.model=NULL,AJUST=0,list.models=NULL,correction=F)
m6 <- which.min(as.numeric(mdr2[,6]))
m5 <- which.max(as.numeric(mdr2[,5])**2)
p <- mdr2[m6,6] %>% as.numeric()
t <- mdr2[m5,5] %>% as.numeric()
ps[i] <- p
ts[i] <- t
if (i >= 20) break
}

fid <- pedigree[selection, "fid"]

lm(form,data = genopheno[selection,]) %>% summary()

timestamp()
scene1 <- simResults_mdr(n_rep = 50,P = 0.5,permutation = 20,save = FALSE, sim_control = sim_params1[[2]],ncores = 50)
timestamp()


