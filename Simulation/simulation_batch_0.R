n <- c(1000, 2000, 4000)
snps <- c(50, 100, 200)
maf <- c(0.1, 0.25, 0.5)
p <- c(0.05)
or <- c(1.5, 2, 2.5)
int_lev <- c(2)
int_num <- c(1)

sim_params1 <- expand.grid(n, snps, maf, p, or, int_lev, int_num) %>% 
  setnames(c("n", "snp_num", "maf","p","int_eff","int_lev","int_num")) %>%
  arrange(maf, snp_num, n) %>%
  t() %>%
  as.data.frame()

scene12 <- map(sim_params1[1:2], simResults, sfa_control = list(eta = 0.025), n_rep = 5)

