sim_param1 <- list(n = 1000, snp_num = 50, maf = 0.25, p = 0.05, int_eff = 1.5, int_lev = 2, int_num = 1, m_eff = NULL, mode = NULL)
sim_param2 <- list(n = 1000, snp_num = 50, maf = 0.25, p = 0.05, int_eff = 2, int_lev = 2, int_num = 1, m_eff = NULL, mode = NULL)

sfa_control1 <- list(eta = 0.01)
sfa_control2 <- list(eta = 0.025)
sfa_control3 <- list(eta = 0.05)

scene10 <- episfa_sim(n_rep = 5, recursion = 2,cvfolds = 10, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "ebic",sfa_control1)
scene101 <- episfa_sim(n_rep = 10, recursion = 2,cvfolds = 10, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "kl",sfa_control1)
scene102 <- episfa_sim(n_rep = 10, recursion = 2,cvfolds = 10, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "hbic",sfa_control1)


scene11 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "ebic",sfa_control1)
scene12 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "ebic",sfa_control2)
scene13 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param1, criteria = "ebic",sfa_control3)
scene21 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param2, criteria = "ebic",sfa_control1)
scene22 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param2, criteria = "ebic",sfa_control2)
scene23 <- episfa_sim(n_rep = 100, recursion = 2,cvfolds = 5, sim_func = simPopLE_l2_sp, sim_control = sim_param2, criteria = "ebic",sfa_control3)


write_rds(list(scene11,scene12,scene13,scene21,scene22,scene23), "scene1_result.rds")

