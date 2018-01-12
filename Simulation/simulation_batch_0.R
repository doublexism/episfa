n <- c(1000, 2000, 4000)
snps <- c(50, 100, 200)
maf <- c(0.1, 0.25, 0.5)
p <- c(0.05)
or <- c(1.5, 2, 2.5)
int_lev <- c(2)
int_lev2 <- c(3,4)
int_num <- c(1)

sim_params1 <- expand.grid(n, snps, maf, p, or, int_lev, int_num) %>% 
  setnames(c("n", "snp_num", "maf","p","int_eff","int_lev","int_num")) %>%
  arrange(maf, snp_num, n) %>%
  t() %>%
  as.data.frame()

sim_params2 <- expand.grid(n, snps, maf, p, or, int_lev2, int_num) %>% 
  setnames(c("n", "snp_num", "maf","p","int_eff","int_lev","int_num")) %>%
  arrange(maf, snp_num, n) %>%
  t() %>%
  as.data.frame()


scene12 <- map(sim_params1, simResults, sfa_control = list(eta = 0.025))
scene13_4 <- map(sim_params2, simResults, sfa_control = list(eta = 0.005), recursion = 1)

write_rds(scene12,"scene12.rds")
alpha1 <- map_dbl(scene12, `[[`, "power3")

