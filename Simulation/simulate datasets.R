source("simulation\\simulate data LE.R")

dis_model <- c("cox")
model <- c("additive", "dominant", "recessive")
n_sib <- c(500, 1000)
MAF <- c(0.1,  0.4)
snpNum <- c(10, 50, 200)
margin <- c(1, 1.05, 1.20)
inter<- c(1.5)

combs <- expand.grid(dis_model, model, n_sib, MAF, snpNum, margin, inter,stringsAsFactors = FALSE)%>% 
  setNames(c("model","gmodel","nsib","maf","nsnp","margin","inter")) %>%
  arrange(model,inter, gmodel,desc(nsib),maf,nsnp, margin)
combs$name <- pmap_chr(combs, paste, sep = '-')

combs1 <- combs[1:36,]

## simulate first 36 scenorios
pmap(combs1[17:36,], ~sib_sim(n=1000, 
                     name = ..8,
                     model = ..1, 
                     genetic_model = ..2, 
                     N = ..3,
                     MAF = ..4,
                     num_SNP = ..5,
                     margin_effect = ..6,
                     interaction_effect = ..7
                     ))
