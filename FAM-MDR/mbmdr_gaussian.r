# MB-MDR code by Victor
# modified to handle Gaussian data: family=Gaussian, link=id
# solved generated warnings: exclude=NA, without ""
# added check for e being both 0 and 1 in first step
# possibly turn e missing into 0 in second step (commented out for now)

MBMDR <- function(Y,EXPO,ESTRAT,SNPS,PVAL=0.1,dimen,first.model=NULL,output="mbmdr.out",AJUST=0,list.models=NULL,correction=T)
{
  
  #Functions
  NextModel <- function(models=NULL){
    if(is.null(models)){
      switch(class(list.models),
             NULL = { ifelse(is.null(first.model), models <- list(c(ncol(SNPS):(ncol(SNPS)-dimen+1)),1,0,0), models <- list(first.model,1,0,0)) },
             integer = { models <- list(list.models,1,1,list.models) },
             numeric = { models <- list(list.models,1,1,list.models) },
             matrix = {	models <- list(list.models[1,],1,nrow(list.models),list.models) },
             character = {	aux <- as.matrix(read.table(list.models))
             models <- list(aux[1,],1,nrow(aux),aux)}
      )
    }
    else{
      if(models[[3]]>0){
        if(models[[2]]==models[[3]]) return(NULL)
        models[[2]] <- models[[2]] + 1
        models[[1]] <- models[[4]][models[[2]],]
      }
      else{
        n <- length(models[[1]])
        model <- models[[1]][n:1]
        for(i in 1:n){
          if (model[i]>i){
            k <- model[i]-1
            model[1:i] <- c((k-i+1):k)
            models[[1]] <- model[n:1]
            lastmodel <- FALSE
            break
          }
          lastmodel <- TRUE
        }
        if(lastmodel) models <- NULL
      }
    }
    return(models)
  }
  
  Regresion <- function(resp,exposed,datos=NULL,estrat=NULL,ajust,correction){
    # tab <- table(resp,exposed,exclude=NA)
    nab <- c(n_sub - sum(exposed), sum(exposed))
    if (min(nab) <= 1){
      return(c(0,1))
    } else if (min(nab) <= 1000){
      r <- lm(resp ~ exposed) %>% summary()
      return(r$coef["exposed",][3:4])
    }
    
    mu <- mean(resp[exposed == 1]) - mean(resp[exposed == 0])
    v <- var(resp[exposed == 1])/nab[2] + var(resp[exposed == 0])/nab[1]
    t <- mu/sqrt(v)
    df_t <- v**2/((var(resp[exposed])/nab[2])**2/(nab[2]-1) +  (var(resp[!exposed])/nab[1])**2/(nab[1]-1))
    p <- ifelse(t >0, (1 - pt(t,df_t))*2,pt(t,df_t)*2)
    return(c(t,p))
    
    
    if(correction ){
      if(ajust==0) auxgl <- logistf(resp~.,data=data.frame(resp,exposed),pl=F)
      if(ajust>0){
        if(ajust==1) auxgl <- logistf(resp~.,data=data.frame(resp,exposed,estrat),pl=F)
        if(ajust==2) auxgl <- logistf(resp~.,data=data.frame(resp,exposed,datos),pl=F)
        if(ajust==3) auxgl <- logistf(resp~.,data=data.frame(resp,exposed,datos,estrat),pl=F)
      }
      return(auxgl)
    }
    else{
      if(ajust==0) auxgl <- lm(resp~.,data=data.frame(resp,exposed))
      if(ajust>0){
        if(ajust==1) auxgl <- lm(resp~.,data=data.frame(resp,exposed,estrat))
        if(ajust==2) auxgl <- lm(resp~.,data=data.frame(resp,exposed,datos))
        if(ajust==3) auxgl <- lm(resp~.,data=data.frame(resp,exposed,datos,estrat))
      }
      return(auxgl)
    }
  }
  
  Exposed <- function(snps,genotip){
    k <- 1
    genotip <- as.character(genotip)
    e <- EXPO[[genotip[k]]][,snps[k]]
    while(k<length(snps)){
      k <- k+1
      e <- as.logical(e * EXPO[[genotip[k]]][,snps[k]])
    }
    return(e)
  }
  
  Exposition <- function(snps,expolist){
    switch(class(expolist),
           matrix = { es <- apply(expolist, 1, function(x){Exposed(snps,x)})
           ifelse(ncol(es)>1, e <- as.logical(rowSums(es)),  e <- es)},
           numeric = { e <- Exposed(snps,expolist)}
    )
    return(e)
  }
  
  FirstStep <- function(model){
    GenotipReg <- function(x){
      #      if(tab[x[1]]==0) return(c(0,0,0,1))
      if(tab[x[1]]==0) return(c(0,1))
      e <- 1*Exposition(model,x[-1])
      if (sum(e==0,na.rm=TRUE)==0 | sum(e==1,na.rm=TRUE)==0) return(c(0,1))
      r <- Regresion(resp=Y,exposed=e,datos=SNPS[,model],estrat=ESTRAT,ajust=AJUST,correction=correction)
      if (class(r)[1]=="logistf"){
        coefic <- r$coefficients[2]
        sd.coef <- diag(r$var)[2]^0.5
        return(c(coefic, sd.coef, coefic/sd.coef, r$prob[2]))
      }
      #     return(summary(r)$coef["exposed",])
      return(r)
    }   
    tab <- table(SNPS[,model],exclude=NA)
    #    print(tab)
    dimens <- dim(tab)
    n <- length(dimens)
    part <- matrix(,prod(dimens),(n + 3))
    part[,1] <- 1:prod(dimens)
    aux1 <- c(dimens,1)
    aux2 <- c(1,dimens)
    for (i in 1:n){
      #    part[,i] <- rep(rep(as.numeric(dimnames(tab)[[i]]),rep(prod(aux1[(i+1):(n+1)]),length(dimnames(tab)[[i]]))),prod(aux2[1:i]))
      #    part[,i] <- rep(rep(1:dimens[i],rep(prod(aux1[(i+1):(n+1)]),length(dimnames(tab)[[i]]))),prod(aux2[1:i]))
      part[,i+1] <- rep(rep(as.numeric(dimnames(tab)[[i]]),rep(prod(aux2[1:i]),length(dimnames(tab)[[i]]))),prod(aux1[(i+1):(n+1)]))
      #    part[,i] <- rep(rep(1:dimens[i],rep(prod(aux2[1:i]),length(dimnames(tab)[[i]]))),prod(aux1[(i+1):(n+1)]))
    }
    reg <- apply(part,1,function(x) GenotipReg(x[1:(n+1)]))
    part[,(n+2):(n+3)] <- t(reg)
#    print(part)
    return(part)
  }
  
  
  #Main program
  #  result <- data.frame(stringsAsFactors = FALSE)
  n_sub <- length(Y)
  models <- NextModel()
  model <- models[[1]]
  npair <- choose(ncol(SNPS),dimen)
  #times(npair) %do% {
  result <-   foreach(preach = 1:npair, .combine = "rbind", .verbose = FALSE) %do% {
    result.l <- NULL
    result.h <- NULL
    part <- FirstStep(model)
    #second step
    h.list <- part[(part[,(dimen+2)]>0 & part[,(dimen+3)]<=PVAL),(2:(1+dimen))]
    l.list <- part[(part[,(dimen+2)]<0 & part[,(dimen+3)]<=PVAL),(2:(1+dimen))]
    if(!is.na(l.list[1])){
      #      print("L")
      switch(class(l.list), matrix={N <- nrow(l.list)}, numeric={N <- 1})
      l.e <- 1*Exposition(model,l.list)
      #l.e=ifelse(is.na(l.e),0,l.e)
      l.r <- Regresion(resp=Y,exposed=l.e,datos=SNPS[,model],estrat=ESTRAT,ajust=AJUST,correction=correction)
      if (class(l.r)[1]=="logistf"){
        coefic <- l.r$coefficients[2]
        sd.coef <- diag(l.r$var)[2]^0.5
        l.regout <- c(coefic, sd.coef, coefic/sd.coef, l.r$prob[2])
      }
      else{
        #      l.regout <- summary(l.r)$coef["exposed",]
        l.regout <- l.r
      }
      #next model
      result.l <-  c(model,"L",N,l.regout)
    }
    
    if(!is.na(h.list[1])){
      #      print("H")
      switch(class(h.list), matrix={N <- nrow(h.list)}, numeric={N <- 1})
      h.e <- 1*Exposition(model,h.list)
      h.r <- Regresion(resp=Y,exposed=h.e,datos=SNPS[,model],estrat=ESTRAT,ajust=AJUST,correction=correction)
      if (class(h.r)[1]=="logistf"){
        coefic <- h.r$coefficients[2]
        sd.coef <- diag(h.r$var)[2]^0.5
        h.regout <- c(coefic, sd.coef, coefic/sd.coef, h.r$prob[2])
      }
      else{
        #        h.regout <- summary(h.r)$coef["exposed",]
        h.regout <- h.r
      }
      #next model
      result.h <- c(model,"H",N,h.regout)
    }
    #    print(model)
    
    models <- NextModel(models)
    model <- models[[1]]
    return(rbind(result.l, result.h))
  }
  
  return(result)
}