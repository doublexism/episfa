n <- c(1000, 2000, 4000)
snps <- c(50, 100, 200)
snps2 <- c(50, 100)
maf <- c(0.1, 0.25, 0.5)
p <- c(0.05)
or <- c(1, 1.5, 2, 2.5)
or2 <- c(1, 2.5)
int_lev <- c(3)
int_lev2 <- c(4)
int_num <- c(1)

sim_params1 <- expand.grid(n, snps2, maf, p, or2, int_lev, int_num) %>% 
  setnames(c("n", "snp_num", "maf","p","int_eff","int_lev","int_num")) %>%
  arrange(int_lev, maf, snp_num, n, int_eff) %>%
  t() %>%
  as.data.frame()

scene13_mdr <- list()
for (batch in 1:ncol(sim_params1)){
  sim_param <- sim_params1[[batch]]
  
  if (sim_param[5] == 1){
    null_p <- simResults_mdr(sim_param,null_p =NULL, n_rep = 50, P = 0.1,ncores = 50, sim_func = simPopLE_l2_sp, save = TRUE, verbose = TRUE)
  } else {
    scene <- simResults_mdr(sim_param, null_p, n_rep = 50,  P = 0.1,ncores = 50, sim_func = simPopLE_l2_sp, save = TRUE, verbose = TRUE)
    scene13_mdr[[batch]] <- scene
  }
}

