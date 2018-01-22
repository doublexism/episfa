# data files
pedfile="pedigree1.ped"
phenofil="pheno1.dat"
mapfile="simulation.map"

# solver settings
npermute = 1000                    # fixed number of permutations, 1000 is recommended
alpha = 0.05                       # type-I error
Pval = 0.1                         # significance level at first step of mbmdr
set.seed(1982)                     # for reproducibility

# loading libraries and code
library(kinship2)
library(GenABEL)
source("mbmdr_gaussian.r")

# removing all possible existing files first as a safety measure
file.remove("permoutput.txt")
file.remove("simulation.raw")
file.remove("permuted.txt")
file.remove("permchisq.txt")
file.remove("rawoutput.txt")

# loading data and bringing in GenABEL format
rawfile="simulation.raw"
convert.snp.ped(pedfile, mapfile, rawfile)
simulation.GenABEL = load.gwaa.data(phenofile = phenofil, genofile = rawfile, force=F,makemap=F,sort=F)
pedigree=read.table(pedfile)
pedsize=nrow(pedigree)
nsnps=(ncol(pedigree)-6)/2

# minor allele count and handling missing genotype data 
allelic = function(k){
geno=pedigree[,(5+2*k):(6+2*k)]
allelic=rowSums(geno==2)-(geno[,1]==0 & geno[,2]==0)            # -1 for missing, 0,1,2 gives count of variant allele
}

# preparing MB-MDR
SNPS = matrix(0,nrow=pedsize,ncol=nsnps)
for (k in 1:nsnps){
 SNPS[,k] = allelic(k)
}
SNPS = as.data.frame(SNPS)
SNPS[SNPS==-1] = NA

# GRAMMAR calculations
# here one has to include main effect and/or covariate adjustments in the polygenic model statement
pkin = kinship(pedigree[,2],pedigree[,3],pedigree[,4])
maineff1 = as.factor(SNPS[,1])  
maineff2 = as.factor(SNPS[,2])  
Yfit = polygenic(trait~maineff1+maineff2,pkin,simulation.GenABEL,trait.type="gaussian")
resi = Yfit$pgresidualY

# discarding individuals with missing new phenotypes
pedsize=sum(!is.na(resi))
SNPS=SNPS[!is.na(resi),]
resi=resi[!is.na(resi)]
                  
# calling modified version of Victor's MB-MDR code
Dim = 2
ajust = 0
rawoutputfile = paste("rawoutput.txt")
correction = F      # routine logistf will not be used
Y = resi
for (j in 1:ncol(SNPS)) SNPS[,j] <- as.factor(SNPS[,j])
EXPO <- list()
EXPO[["0"]] <- (SNPS==0)
EXPO[["1"]] <- (SNPS==1)
EXPO[["2"]] <- (SNPS==2)
MBMDR(Y,EXPO,ESTRAT,PVAL=Pval,dimen=Dim,first.model=NULL,output=rawoutputfile,AJUST=ajust,list.models=NULL,correction=correction)

# from row index to snp pair
npairs = nsnps*(nsnps-1)/2
pair = matrix(0,npairs,2)
index=1
for (i in nsnps:2){
 for (j in (i-1):1){
  pair[index,]=c(i,j)
  index=index+1
 }
}

# from snp pair to index
index = function(pair){
 n = pair[1]
 m = pair[2]
 count = nsnps*(nsnps-1)/2 - n*(n-1)/2 + n - m
 return(count)
} 

# maximum chisquare test statistic
results = read.table(rawoutputfile,header=F,sep=";")
Wald = results$V7
chisq = Wald^2
chisqmax = max(chisq)
pairfound = results[which.max(chisq),1:2]
indexfound = as.numeric(index(pairfound))

# fixed number of permutations
permute = matrix(nrow=npermute,ncol=pedsize)
for (k in 1:npermute){
 permute[k,] = sample(1:pedsize)
}
permoutfile = "permuted.txt"
permoutchisq = "permchisq.txt"
chisqmaxperm = numeric(npermute)

# permutations and their maximum chisquare test statistics
for (p in 1:npermute){
 Y = resi[permute[p,]]
 write.table(NULL,file=permoutfile,quote=F)
 MBMDR(Y,EXPO,ESTRAT,PVAL=Pval,dimen=Dim,first.model=NULL,output=permoutfile,AJUST=ajust,list.models=NULL,correction=correction)
 resultsperm = read.table(permoutfile,sep=";")
 Waldperm = resultsperm$V7
 chisqperm = Waldperm^2
 resultsperm=cbind(resultsperm,chisqperm)
 chisqpermbypair=numeric(npairs)
 for (k in 1:npairs){
  chisqpermbypair[k]=max(resultsperm[resultsperm[,1]==pair[k,1]&resultsperm[,2]==pair[k,2],9],0)
 }
 write.table(t(chisqpermbypair),file=permoutchisq,row.names=F,col.names=F,append=T)
 chisqmaxperm[p] = max(chisqpermbypair)
}

# permutation p-value of maximal test statistic
permpval = length(chisqmaxperm[chisqmaxperm>chisqmax]) / npermute
resultsperm = cbind(pairfound,permpval)
write.table(resultsperm,file="permoutput.txt",sep="\t",row.names=F,col.names=F)

# remove unnecessary files
file.remove(rawfile)
file.remove(permoutfile)