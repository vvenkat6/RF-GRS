# Script for simulating data and evaluating results for Polygenic Risk score
# Training and Validation Data is got from selecting chr21 from the TWB dataset
# Base Data is got from chr21 from TWB Dataset

# check if all necessary packages are installed
#list.of.packages <- c("randomForest","fmsb","pROC","glmnet","gtools","Rlab","Metrics","dplyr","ranger","bigsnpr","data.table","magrittr","lassosum","methods","parallel","measures","tidyr")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# Libraries Used
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(fmsb))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(Rlab))
suppressPackageStartupMessages(library(Metrics))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(bigsnpr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(lassosum))
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(measures))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(neuralnet))
suppressPackageStartupMessages(library(sigmoid))
suppressPackageStartupMessages(library(keras))
suppressPackageStartupMessages(library(tensorflow))
suppressPackageStartupMessages(library(tidyverse))

# Function for simulation of phenotype
SimPheno <- function(Data, model = 1, causalN = 100,strength= 1,causal,seed = 1){
	# Key for model: 1 - No interaction, 2 - Interaction , 3 - Non-linear interaction, 4 - Combination
	# Base data will be written to a file outside and Target data will be returned from this function
	# b0 needs to be picked in a way to have a balanced dataset
	print(seed)
        print("Simulating phenotypes")
	if (strength == 1){
	b1log = 1.5
	b0 = 0
	}
	if (strength == 2){
        b1log = 1.1
	b0 = 0
        }
	if (strength == 3){
        b1log = 1.002
	b0 = 0
        }

	Data <- scale(Data[])
	pi <- rep(NA,dim(Data)[1])
	set.seed(seed)
        b1 = rnorm(causalN,mean=log(b1log),sd=0.1)
        # Randomly make half of them -ve
        set.seed(seed)
        b1 = b1 * sample(c(-1,1),causalN,replace =TRUE)
	
	if(model == 1){
		# Non-interacting model
		set.seed(seed)
		sum <- rowSums(sweep(Data[,causal], MARGIN = 2, b1,"*"))
		pi <- (exp(b0 + sum) / (1+exp(b0 + sum)))
		set.seed(seed)
		Pheno <- rbern(dim(Data)[1],pi)
		return(Pheno)	
	}
	if(model == 2){
		# This model indicates only interaction model
		set.seed(seed)
		bdata.xaction <- matrix(, nrow = nrow(Data), ncol =((causalN)/2)*((causalN)/2))
		for (i in 1:nrow(Data)){
                	bdata.xaction[i,] <- Data[i,causal[1:50]] %*% t(Data[i,causal[51:100]])
        	}
		set.seed(seed)
		causal.xaction <- sort(sample(ncol(bdata.xaction), size = ncol(bdata.xaction)/2))
		set.seed(seed)
		b3 <- rnorm(length(causal.xaction),mean=log(b1log + 0.5),sd=0.2)
		set.seed(seed)
		b3 <- b3 * sample(c(-1,1),length(causal.xaction),replace =TRUE)	
		sum<- rowSums(sweep(bdata.xaction[,causal.xaction],MARGIN = 2, b3,"*"))
		pi <- (exp(b0 + sum) / (1+exp(b0 + sum)))
		set.seed(seed)
		Pheno <- rbern(dim(Data)[1],pi)
		return(Pheno)
	}
	if(model == 3){
                # This model indicates only interaction model
                set.seed(seed)
                bdata.xaction <- matrix(, nrow = nrow(Data), ncol =((causalN)/2)*((causalN)/2))
                for (i in 1:nrow(Data)){
                        bdata.xaction[i,] <- Data[i,causal[1:50]] %*% t(Data[i,causal[51:100]])
                }
                set.seed(seed)
                causal.xaction <- sort(sample(ncol(bdata.xaction), size = ncol(bdata.xaction)/2))
                set.seed(seed)
                b3 <- rnorm(length(causal.xaction),mean=log(b1log + 0.5),sd=0.2)
                set.seed(seed)
                b3 <- b3 * sample(c(-1,1),length(causal.xaction),replace =TRUE)
                sum<- rowSums(sweep(bdata.xaction[,causal.xaction],MARGIN = 2, b3,"*")) + rowSums(sweep(Data[,causal[1:50]],MARGIN = 2, b1[1:50],"*"))
                pi <- (exp(b0 + sum) / (1+exp(b0 + sum)))
                set.seed(seed)
                Pheno <- rbern(dim(Data)[1],pi)
                return(Pheno) 
	}
}

args = commandArgs(trailingOnly=TRUE)
seed = as.numeric(args[1])
modelN = as.numeric(args[2])
causalN = as.numeric(args[3])
strength = as.numeric(args[4])
#Sys.sleep(seed * 5 * modelN * strength)
#seed = seed + 100
#modelN = 3
causalN =100
system(paste("mkdir -p logs",paste(seed,modelN,sep="_"),sep=""))

# Read in data
tmpfile <- tempfile()
bedfile = "/home3/DrYao/FINALest/Simulation/SimFiles/TWB_21.bed"
snp_readBed(bedfile, backingfile = tmpfile)
obj.SNP <- snp_attach(paste0(tmpfile, ".rds"))
str(obj.SNP, max.level = 2, strict.width = "cut")
Data   <- obj.SNP$genotypes

setwd(paste("logs",paste(seed,modelN,sep="_"),sep=""))
print("changin directory")

set.seed(seed)
causal <- sort(sample(ncol(Data), size = causalN))
#causal <- sort(sample(ncol(Data), size = causalN/2))
#causal <- sort(c(causal,causal+1))

# Call function to simulate phenotype
Pheno <-SimPheno(Data,modelN,causalN,strength,causal,seed)
print(table(Pheno))
split = table(Pheno)[1]/(table(Pheno)[1]+table(Pheno)[2])
nCausal <- c(1:ncol(Data))[-causal]
causal <- sort(c(causal,nCausal))

# Separate data into Target and Base Data
a <- which(Pheno == 1)
set.seed(seed)
sel.a <- sample(a,500)
b <- which(Pheno == 0)
set.seed(seed)
sel.b <- sample(b,500)
sel <- sort(c(sel.a,sel.b))
tsel <- obj.SNP$fam[sel,][1:2]
tData <- Data[sel,causal]
tPheno <- Pheno[sel]
write.table(tsel,file="tsel",row.names=F,col.names=F,quote=F)
write.table(obj.SNP$map[causal,]$marker.ID,file="causalSNPs",row.names=F,col.names=F,quote=F)

system("plink --bfile /home3/DrYao/SIM2023/SimFiles/TWB_21 --allow-no-sex --keep tsel --out TargetData --make-bed --extract causalSNPs")
system("plink --bfile /home3/DrYao/SIM2023/SimFiles/TWB_21 --allow-no-sex --remove tsel --out BaseData --make-bed --extract causalSNPs")
# Write Pheno to files
dat <- read.table("BaseData.fam")
dat$V6 <- Pheno[-sel] + 1
write.table(dat,file="BaseData.fam",quote = F,row.names = F,col.names = F,sep=" ")
dat <- read.table("TargetData.fam")
dat$V6 <- tPheno + 1
write.table(dat,file="TargetData.fam",quote = F,row.names = F,col.names = F,sep=" ")

## Make summary statistic file
system("module load plink")
system("plink --bfile BaseData --logistic beta --ci 0.5 --allow-no-sex --out PLINK")

# Pre-processing
# Create and read pc files
system("plink --bfile TargetData --pca 5 --out TargetData")
phenotype <- read.table("TargetData.fam", header=F)
colnames(phenotype) <- c("FID","IID","a1","a2","a3","PHENO")
pcs <- read.table("TargetData.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:5))
pheno <- merge(phenotype, pcs, by=c("FID","IID"))

# Read in data
tmpfile <- tempfile()
bedfile = "TargetData.bed"
snp_readBed(bedfile, backingfile = tmpfile)
obj.SNP <- snp_attach(paste0(tmpfile, ".rds"))
G <- obj.SNP$genotypes
CHR <- obj.SNP$map$chromosome
POS <- obj.SNP$map$physical.pos
sumstats <- bigreadr::fread2("PLINK.assoc.logistic")
sumstats$lpval <- -log10(sumstats$P)

# Split data into train and validation set varying with seed
train.ratio <- 0.8
set.seed(seed)
ind.train <- sample(nrow(G), round(train.ratio*nrow(G)))
ind.test <- setdiff(rows_along(G), ind.train)

NCORES = 4
#################################################################################################################################################
# C+T prep
gwas.train <- big_univLogReg(G,tPheno[ind.train],ind.train = ind.train)
gwas_gc <- snp_gc(gwas.train)
p.value.all <- predict(gwas.train, log10 = FALSE)
thrs <- sort(c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100))))
lpS <- -predict(gwas_gc)

# p-value grid:
gwas_gc$lpS <- lpS
rownames(gwas_gc) <- obj.SNP$map[,2]

#################################################################################################################################################
# Stacked C+T (indep)
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS, exclude = which(is.na(lpS)),ncores = NCORES)
multi_PRS <- snp_grid_PRS(G, all_keep, gwas_gc$estim, lpS = lpS, ind.row = ind.train,ncores=NCORES)
final_mod <- snp_grid_stacking(multi_PRS, tPheno[ind.train], ncores = NCORES, K = 4)
new_beta <- final_mod$beta.G
ind <- which(new_beta != 0)
pred <- final_mod$intercept + big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)
i.SCT.auc.a <- round(AUCBoot(pred, tPheno[ind.test]),4)[[1]]
i.SCT.auc.b <- round(AUCBoot(pred, tPheno[ind.test]),4)[[2]]
i.SCT.auc.c <- round(AUCBoot(pred, tPheno[ind.test]),4)[[3]]
i.SCT.n <- length(ind)
i.SCT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#################################################################################################################################################
# Traditional C+T (indep) [also known as maxCT in Prive paper]
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
  tidyr::unnest(cols = "thr.lp")
s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
        single_PRS <- X[, ind]
        bigstatsr::AUC(single_PRS, y.train)
        }, ind = 1:s, s = s, y.train = tPheno[ind.train],a.combine = 'c', block.size = 1, ncores = NCORES)
max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
pred <- snp_PRS(G, gwas_gc$estim[ind.keep], ind.test = ind.test, ind.keep = ind.keep,lpS.keep = lpS[ind.keep], thr.list = max_prs$thr.lp)
i.CT.auc.a <- round(AUCBoot(pred,tPheno[ind.test]),4)[[1]]
i.CT.auc.b <- round(AUCBoot(pred,tPheno[ind.test]),4)[[2]]
i.CT.auc.c <- round(AUCBoot(pred,tPheno[ind.test]),4)[[3]]
i.CT.n <- length(ind.keep)
i.CT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#################################################################################################################################################
# PCA-PRS indep : For CT and p-value thresold used in Stacking
multi_PRS <- snp_grid_PRS(G, all_keep, gwas_gc$estim, lpS = lpS, ind.row = ind.test,ncores=NCORES)
corr_matrix <- cor(t(cbind(multi_PRS[])))
pca.prs <- princomp(corr_matrix)$loadings[,1]
i.pca.ct.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
i.pca.ct.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
i.pca.ct.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
i.pca.ct.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
i.pca.ct.n <- dim(G)[2]

#################################################################################################################################################
# PCA-PRS indep: Only p-value and fixing R2 to 0.2
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS, exclude = which(is.na(lpS)),ncores = NCORES,grid.thr.r2 = c(0.2),grid.base.size=c(50))
multi_PRS <- snp_grid_PRS(G, all_keep, gwas_gc$estim, lpS = lpS, ind.row = ind.test,ncores=NCORES)
corr_matrix <- cor(t(cbind(multi_PRS[])))
pca.prs <- princomp(corr_matrix)$loadings[,1]
i.pca.ct.p.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
i.pca.ct.p.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
i.pca.ct.p.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
i.pca.ct.p.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
i.pca.ct.p.n <- dim(G)[2]

#################################################################################################################################################
# Penalized linear regression (indep)
logit <- big_spLogReg(X = G, y01.train = tPheno[ind.train],alphas = c(1,0.5,0.05,0.001),ind.train = ind.train,ncores=NCORES)
pred1 <- predict(logit, X = G, ind.row = ind.test)

plr.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4)
plr.auc.a <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[1]]
plr.auc.b <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[2]]
plr.auc.c <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[3]]
plr.n <- summary(logit)$nb_var[which(summary(logit)$validation_loss == min(summary(logit)$validation_loss))]

#################################################################################################################################################
# SIS + PLR (SS)
t <- lpS
t.ind <- sort(abs(t), decreasing=T, index.return=T, na.last = TRUE)$ix
screen.ind.sis <-  t.ind[1:length(ind.train)-1]
set.seed(seed)
K.fold <- sample(rep_len(1:5, length(ind.train)))

X.sis <- FBM(n <- dim(G)[1], length(screen.ind.sis), init = G[, screen.ind.sis])
mod.sis <- big_spLogReg(X.sis, tPheno[ind.train], K=5, ind.train = ind.train, ind.sets = K.fold,ncores=NCORES)
pred.sis <- predict(mod.sis, X.sis, ind.test)

sis.lasso.auc.a <- round(AUCBoot(pred.sis,tPheno[ind.test]),4)[[1]] 
sis.lasso.auc.b <- round(AUCBoot(pred.sis,tPheno[ind.test]),4)[[2]]
sis.lasso.auc.c <- round(AUCBoot(pred.sis,tPheno[ind.test]),4)[[3]]
sis.lasso.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred.sis,family = "binomial"))$R2, 4)
sis.lasso.n <- summary(mod.sis)$nb_var

#################################################################################################################################################
# Stochastic Gradient Descent to see if this line can be fit better
sgd <- function(
  par,                       # parameter estimates
  X,                         # model matrix
  y,                         # target variable
  stepsize = 1e-2,           # the learning rate; suggest 1e-3 for non-adagrad methods
  type = 'adagrad',          # one of adagrad, rmsprop, adam or nadam
  average = FALSE,           # a variation of the approach
  ...                        # arguments to pass to an updating function, e.g. gamma in rmsprop
){
  
  # initialize
  beta = par
  names(beta) = colnames(X)
  betamat = matrix(0, nrow(X), ncol = length(beta))      # Collect all estimates
  v    = rep(0, length(beta))                    # gradient variance (sum of squares)
  m    = rep(0, length(beta))                    # average of gradients for n/adam
  eps  = 1e-8                                    # a smoothing term to avoid division by zero
  grad_old = rep(0, length(beta))
  
  update_ff <- function(type, ...) {
    
    # if stepsize_tau > 0, a check on the LR at early iterations
    adagrad <- function(grad, stepsize_tau = 0) {
      v <<- v + grad^2  
      
      stepsize/(stepsize_tau + sqrt(v + eps)) * grad
    }
    
    rmsprop <- function(grad, grad_old, gamma = .9) {
      v = gamma * grad_old^2 + (1 - gamma) * grad^2
      
      stepsize / sqrt(v + eps) * grad
    }
    
    adam <- function(grad, b1 = .9, b2 = .999) {
      m <<- b1 * m + (1 - b1) * grad
      v <<- b2 * v + (1 - b2) * grad^2
      
      if (type == 'adam')
        # dividing v and m by 1 - b*^i is the 'bias correction'
        stepsize/(sqrt(v / (1 - b2^i)) + eps) *  (m / (1 - b1^i))
      else 
        # nadam
        stepsize/(sqrt(v / (1 - b2^i)) + eps) *  (b1 * m  +  (1 - b1)/(1 - b1^i) * grad)
    }
    
    switch(
      type,
      adagrad = function(grad, ...) adagrad(grad, ...),
      rmsprop = function(grad, ...) rmsprop(grad, grad_old, ...),
      adam    = function(grad, ...) adam(grad, ...),
      nadam   = function(grad, ...) adam(grad, ...)
    )
  }

  update = update_ff(type, ...)
  
  for (i in 1:nrow(X)) {
    Xi   = X[i, , drop = FALSE]
    yi   = y[i]
    LP   = Xi %*% beta                           # matrix operations not necessary, 
    grad = t(Xi) %*% (LP - yi)                   # but makes consistent with standard gd func
    
    # update
    beta = beta - update(grad, ...)
    
    if (average & i > 1) {
      beta = beta - 1/i * (betamat[i - 1, ] - beta)   # a variation
    } 
    
    betamat[i,] = beta
    grad_old = grad
  }
  
  LP = X %*% beta
  lastloss = crossprod(LP - y)
  
  list(
    par = beta,                               # final estimates
    par_chain = betamat,                      # estimates at each iteration
    RMSE = sqrt(sum(lastloss)/nrow(X)),
    fitted = LP
  )
}

x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]
init = summary(mod.sis)$beta[1][[1]]
names(init) <- obj.SNP$map[,2][screen.ind.sis]
fit_gd = sgd(init,X = x[ind.train,screen.ind.sis],y = tPheno[ind.train],stepsize=1e-3,type="adam")
pred <- big_prodVec(G, fit_gd$par, ind.row = ind.test, ind.col = screen.ind.sis)
sis.lasso.sgd.auc.a <- round(AUCBoot(pred, tPheno[ind.test]),4)[[1]]
sis.lasso.sig.auc.b <- round(AUCBoot(pred, tPheno[ind.test]),4)[[2]]
sis.lasso.sgd.auc.c <- round(AUCBoot(pred, tPheno[ind.test]),4)[[3]]
sis.lasso.sgd.n <- summary(mod.sis)$nb_var
sis.lasso.sgd.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#################################################################################################################################################
# SIS + RF

x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]

p <- length(screen.ind.sis)
hyper_grid <- expand.grid(
  mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,sqrt(p)*5,0.1 * p,0.5 * p,p),
  node_size  = c(1,20,40,80,100),
  num.trees = c(250, 500, 1000),
  OOB_RMSE   = 0
)

print("normal method")
# normal method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,screen.ind.sis])     
rfNR2.2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.2a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.2b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.2c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

print("ic method")
# impurity_corrected method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,screen.ind.sis])
rfNR2.1 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.1a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.1b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.1c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

sis.rf.auc.a <- c(rfAUC.1a,rfAUC.2a)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rf.auc.b <- c(rfAUC.1b,rfAUC.2b)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rf.auc.c <- c(rfAUC.1c,rfAUC.2c)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rf.nr2 <- max(rfNR2.1 ,rfNR2.2 )
sis.rf.n <- length(screen.ind.sis)

# importance values
rf2 <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=1,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
rf3 <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=2,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
p <- importance_pvalues(rf)
p2 <- importance_pvalues(rf2)
p3 <- importance_pvalues(rf3)
rfImpnr <- NULL
rfImpauc1 <- NULL
rfImpauc2 <- NULL
rfImpauc3 <- NULL
rfImpn <- NULL

# RF Lasso Method and RF Importance Value Method
p <- cbind(p,1:dim(p)[1])
p2 <- cbind(p2,1:dim(p2)[1])
p3 <- cbind(p3,1:dim(p3)[1])

# Picking Best p-values for searching (p-values are inversely proportional to important measures)
for (i in seq(0.201,1.01,0.1)){
	idx1 <- (p[p[,2] <= i,])[,3]
        idx2 <- (p2[p2[,2] <= i,])[,3]
        idx3 <- (p3[p3[,2] <= i,])[,3]
	idx <- intersect(idx1,intersect(idx2,idx3))
	val <- rowMeans(cbind(p[,1],p2[,1],p3[,1]))
        rfMat <- x[ind.test,screen.ind.sis]
        val <- (val+0.0001)/sum(val)
        rfPRS <- rfMat[,idx] %*% val[idx]
        rfImpnr <- rbind(rfImpnr,round(NagelkerkeR2(glm(tPheno[ind.test]~rfPRS,family = "binomial"))$R2, 4))
        rfImpauc1 <- rbind(rfImpauc1,round(AUCBoot(rfPRS,tPheno[ind.test])[[1]],4))
	rfImpauc2 <- rbind(rfImpauc2,round(AUCBoot(rfPRS,tPheno[ind.test])[[2]],4))
	rfImpauc3 <- rbind(rfImpauc3,round(AUCBoot(rfPRS,tPheno[ind.test])[[3]],4))
	rfImpn <- rbind(rfImpn,length(idx))
}

sis.rfimp.auc.a <- max(rfImpauc1)
sis.rfimp.auc.b <- max(rfImpauc2)
sis.rfimp.auc.c <- max(rfImpauc3)
sis.rfimp.nr2 <- max(rfImpnr)
sis.rfimp.n <- rfImpn[which(rfImpnr == max(rfImpnr))]

#################################################################################################################################################
# SIS + RFprob (SS)
# Input the betas/beta_se as a value of how often the snp should be selected [using split.select.weights]
print("Starting RFprob")
weight <- (t[screen.ind.sis] - min( t[screen.ind.sis]) )/ (max( t[screen.ind.sis]) - min( t[screen.ind.sis]))
p <- length(screen.ind.sis)
hyper_grid <- expand.grid(
  mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,sqrt(p)*5,0.1 * p,0.5 * p),
  node_size  = c(1,20,40),
  num.trees = c(250, 500),
  OOB_RMSE   = 0
)

# normal method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),split.select.weights = weight,seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,screen.ind.sis])
rfNR2.2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.2a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.2b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.2c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

# impurity_corrected method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,split.select.weights=weight,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,screen.ind.sis]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",
                num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,screen.ind.sis])
rfNR2.1 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.1a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.1b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.1c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

sis.rfprior.auc.a <- c(rfAUC.1a,rfAUC.2a)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rfprior.auc.b <- c(rfAUC.1b,rfAUC.2b)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rfprior.auc.c <- c(rfAUC.1c,rfAUC.2c)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
sis.rfprior.nr2 <- max(rfNR2.1 ,rfNR2.2 )
sis.rfprior.n <- length(screen.ind.sis)

#################################################################################################################################################
# RFs for different clumping thresholds: relaxed
print("Clumping threshold prep:realaxed")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]

rfnr2 <- NULL
rfauc1 <- NULL
rfauc2 <- NULL
rfauc3 <- NULL
rfn <- NULL

count = 0
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS,grid.thr.r2=c(0.3, 0.6, 0.9, 0.95),grid.base.size = c(100),exclude = which(is.na(lpS)),ncores = NCORES)

for (idx in all_keep$`21`){
        print("All_keep:")
        print(count)
        count = count + 1
	p <- length(idx)
                if (p>3){
                        hyper_grid <- expand.grid(
                                mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,0.1 * p,0.5 * p),
                                num.trees = c(250, 500),
                                OOB_RMSE   = 0
                                )
                        for(i in 1:nrow(hyper_grid)){
                                if(hyper_grid$mtry[i] < length(idx)){
                                        rf <- ranger(x=x[ind.train,idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
                                        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
                                        }
                                else{ hyper_grid$OOB_RMSE[i] <- 100 }
                                }
                        i =  which.min(hyper_grid$OOB_RMSE)
                        rf <- ranger(x=x[ind.train,idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
                        pred.Train <- predict(rf,x[ind.train,idx])$predictions[,2]
                        pred.Test <- predict(rf,x[ind.test,idx])$predictions[,2]
                        rfnr2 <- rbind(rfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred.Test,family = "binomial"))$R2, 4))
                        rfauc1 <- rbind(rfauc1,round(AUCBoot(pred.Test,tPheno[ind.test])[[1]],4))
                        rfauc2 <- rbind(rfauc2,round(AUCBoot(pred.Test,tPheno[ind.test])[[2]],4))
                        rfauc3 <- rbind(rfauc3,round(AUCBoot(pred.Test,tPheno[ind.test])[[3]],4))
                	rfn <- rbind(rfn,p)
		}
}

id = which.max(rfnr2)
v.rf.nr2 <- rfnr2[id]
v.rf.auc.a <- rfauc1[id]
v.rf.auc.b <- rfauc2[id]
v.rf.auc.c <- rfauc3[id]
v.rf.n <- rfn[id]

#################################################################################################################################################
# RFs for different clumping thresholds: rigorous

rfnr2 <- NULL
rfauc1 <- NULL
rfauc2 <- NULL
rfauc3 <- NULL
rfn <- NULL

count = 0
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS,exclude = which(is.na(lpS)),ncores = NCORES)

for (idx in all_keep$`21`){
        print("All_keep:")
        print(count)
        count = count + 1
        p <- length(idx)
                if (p>3){
                        hyper_grid <- expand.grid(
                                mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,0.1 * p,0.5 * p),
                                num.trees = c(250, 500),
                                OOB_RMSE   = 0
                                )
                        for(i in 1:nrow(hyper_grid)){ 
                                if(hyper_grid$mtry[i] < length(idx)){
                                        rf <- ranger(x=x[ind.train,idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i]) 
                                        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
                                        }
                                else{ hyper_grid$OOB_RMSE[i] <- 100 }
                                }
                        i =  which.min(hyper_grid$OOB_RMSE)
                        rf <- ranger(x=x[ind.train,idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
                        pred.Train <- predict(rf,x[ind.train,idx])$predictions[,2]
                        pred.Test <- predict(rf,x[ind.test,idx])$predictions[,2]
                        rfnr2 <- rbind(rfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred.Test,family = "binomial"))$R2, 4))
                        rfauc1 <- rbind(rfauc1,round(AUCBoot(pred.Test,tPheno[ind.test])[[1]],4))
                        rfauc2 <- rbind(rfauc2,round(AUCBoot(pred.Test,tPheno[ind.test])[[2]],4))
                        rfauc3 <- rbind(rfauc3,round(AUCBoot(pred.Test,tPheno[ind.test])[[3]],4))
                        rfn <- rbind(rfn,p)
                }
}

id = which.max(rfnr2)
vs.rf.nr2 <- rfnr2[id]
vs.rf.auc.a <- rfauc1[id]
vs.rf.auc.b <- rfauc2[id]
vs.rf.auc.c <- rfauc3[id]
vs.rf.n <- rfn[id]

#################################################################################################################################################
# RFs for different CT threshold Prep with p-values
print("CT threshold prep")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]

pred.Train.cp <- data.frame(matrix(ncol = 800, nrow = 0))
pred.Test.cp <- data.frame(matrix(ncol = 200, nrow = 0))

rfnr2 <- NULL
rfauc1 <- NULL
rfauc2 <- NULL
rfauc3 <- NULL
rfn <- NULL
rfc <- NULL
rfp <- NULL
rfi <- NULL

count = 0
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS,grid.base.size = c(50, 100, 200),exclude = which(is.na(lpS)),ncores = NCORES)

for (idx in all_keep$`21`){
	print("All_keep:")
	print(count)
        count = count + 1
	p.val <- lpS[idx]
	grid.lpS.thr = 0.9999 * seq_log(max(0.1, min(p.val, na.rm = TRUE)),max(p.val, na.rm = TRUE), 5)
	for (pcut in grid.lpS.thr){
		a <- gwas_gc[idx,]
		names <- rownames(a[a$lpS <= pcut,])
		names.idx <- match(names,obj.SNP$map[,2])
		p <- length(names)
		if (p>3){
			hyper_grid <- expand.grid(
  				mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,0.1 * p,0.5 * p),
  				num.trees = c(250, 500),
  				OOB_RMSE   = 0
				)
			for(i in 1:nrow(hyper_grid)){
				if(hyper_grid$mtry[i] < length(names.idx)){
        				rf <- ranger(x=x[ind.train,names.idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        				hyper_grid$OOB_RMSE[i] <- rf$prediction.error
					}
				else{ hyper_grid$OOB_RMSE[i] <- 100 }
				}
			i =  which.min(hyper_grid$OOB_RMSE)
			rfi <- rbind(rfi,i)
			rf <- ranger(x=x[ind.train,names.idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
			pred.Train <- predict(rf,x[ind.train,names.idx])$predictions[,2]
			pred.Test <- predict(rf,x[ind.test,names.idx])$predictions[,2]
			pred.Train.cp <- rbind(pred.Train.cp,pred.Train)
			pred.Test.cp <- rbind(pred.Test.cp,pred.Test)
			rfnr2 <- rbind(rfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred.Test,family = "binomial"))$R2, 4))
			rfauc1 <- rbind(rfauc1,round(AUCBoot(pred.Test,tPheno[ind.test])[[1]],4))
			rfauc2 <- rbind(rfauc2,round(AUCBoot(pred.Test,tPheno[ind.test])[[2]],4))
			rfauc3 <- rbind(rfauc3,round(AUCBoot(pred.Test,tPheno[ind.test])[[3]],4))
			rfn <- rbind(rfn,p)
			rfc <- rbind(rfc,count)
			rfp <- rbind(rfp,pcut)
		}
	}
}
colnames(pred.Test.cp) <- obj.SNP$fam[,1][ind.test]
colnames(pred.Train.cp) <- obj.SNP$fam[,1][ind.train]

id = which.max(rfnr2)
rf.nr2 <- rfnr2[id]
rf.auc.a <- rfauc1[id]
rf.auc.b <- rfauc2[id]
rf.auc.c <- rfauc3[id]
rf.n <- rfn[id]

rf.c.idx <- rfc[id]
rf.p <- rfp[id]

#################################################################################################################################################
# PCA: CT + RF (SS)
print("Starting CT+RF")
corr_matrix <- cor(pred.Test.cp)
pca.prs <- princomp(corr_matrix)$loadings[,1]
pca.rf.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
pca.rf.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
pca.rf.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
pca.rf.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
z <- length(all_keep$`21`)
pca.rf.n <- length(all_keep$`21`[[z]])

#################################################################################################################################################
# Stacked: CT + RF (SS)
logit <- big_spLogReg(X = as_FBM(t(pred.Train.cp)), y01.train = tPheno[ind.train],alphas=c(1,0.5,0.05,0.001),ncores=NCORES)
pred1 <- predict(logit,X=as_FBM(t(pred.Test.cp)))
stack.rf.auc.a <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[1]]
stack.rf.auc.b <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[2]]
stack.rf.auc.c <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[3]]
stack.rf.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4)
stack.rf.n <- length(all_keep$`21`[[z]])

#################################################################################################################################################
# RF (indep) and RFLR (indep)
# impurity_corrected method
print("Starting impurity corrected method")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]

p <- dim(x)[2]
hyper_grid <- expand.grid(
  mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,0.1 * p,0.5*p),
  node_size  = c(1,10,20,40),
  num.trees = c(250, 500),
  OOB_RMSE   = 0
)
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
	print(i)
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,])
rfic.NR2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfic.AUC.a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfic.AUC.b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfic.AUC.c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)
rfic.n <- p

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# importance values
print("Running RFLR")
rf2 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=1,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
rf3 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=2,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
p <- importance_pvalues(rf)
p2 <- importance_pvalues(rf2)
p3 <- importance_pvalues(rf3)
rfImpnr <- NULL
rflrnr <- NULL
rfImpauc1 <- NULL
rflrauc1 <- NULL
rfImpauc2 <- NULL
rflrauc2 <- NULL
rfImpauc3 <- NULL
rflrauc3 <- NULL
rfImpn <- NULL
rflrn <- NULL
sgd.auc.a<- NULL
sgd.auc.b <- NULL
sgd.auc.c <- NULL
sgd.nr2 <- NULL
sgd.n <- NULL
seq.data <- NULL

# RF Lasso Method and RF Importance Value Method
p <- cbind(p,1:dim(p)[1])
p2 <- cbind(p2,1:dim(p2)[1])
p3 <- cbind(p3,1:dim(p3)[1])

# Picking Best p-values for searching (p-values are inversely proportional to important measures)
for (i in seq(0.201,1.01,0.1)){
        idx1 <- (p[p[,2] <= i,])[,3]
        idx2 <- (p2[p2[,2] <= i,])[,3]
        idx3 <- (p3[p3[,2] <= i,])[,3]
        idx <- intersect(idx1,intersect(idx2,idx3))
	seq.data <- rbind(seq.data,i)
        val <- rowMeans(cbind(p[,1],p2[,1],p3[,1]))
        rfMat <- x[ind.test,]
        val <- (val+0.0001)/sum(val)
        rfPRS <- rfMat[,idx] %*% val[idx]
        rfImpnr <- rbind(rfImpnr,round(NagelkerkeR2(glm(tPheno[ind.test]~rfPRS,family = "binomial"))$R2, 4))
        rfImpauc1 <- rbind(rfImpauc1,round(AUCBoot(rfPRS,tPheno[ind.test])[[1]],4))
        rfImpauc2 <- rbind(rfImpauc2,round(AUCBoot(rfPRS,tPheno[ind.test])[[2]],4))
        rfImpauc3 <- rbind(rfImpauc3,round(AUCBoot(rfPRS,tPheno[ind.test])[[3]],4))
        rfImpn <- rbind(rfImpn,length(idx))

	logit <- big_spLogReg(X = G, y01.train = tPheno[ind.train],alphas = c(1,0.5,0.05,0.001),ind.train = ind.train,ind.col=idx[idx<=dim(as.matrix(x[ind.train,]))[2]])
        pred1 <- predict(logit, X = G, ind.row = ind.test ,ind.col=idx[idx<=dim(as.matrix(x[ind.train,]))[2]])
        rflrnr <- rbind(rflrnr,round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4))
	rflrauc1 <- rbind(rflrauc1,round(AUCBoot(pred1,tPheno[ind.test])[[1]],4))
        rflrauc2 <- rbind(rflrauc2,round(AUCBoot(pred1,tPheno[ind.test])[[2]],4))
        rflrauc3 <- rbind(rflrauc3,round(AUCBoot(pred1,tPheno[ind.test])[[3]],4))	
	rflrn <- rbind(rflrn,summary(logit)$nb_var[which.min(summary(logit)$validation_loss)])

	init = summary(logit)$beta[1][[1]]
	names(init) <- obj.SNP$map[,2][idx]
	fit_gd = sgd(init,X = x[ind.train,idx],y = tPheno[ind.train],stepsize=1e-3,type="adam")
	pred <- big_prodVec(G, fit_gd$par, ind.row = ind.test, ind.col = idx)
	sgd.auc.a <- rbind(sgd.auc.a,round(AUCBoot(pred, tPheno[ind.test]),4)[[1]])
	sgd.auc.b <- rbind(sgd.auc.b,round(AUCBoot(pred, tPheno[ind.test]),4)[[2]])
	sgd.auc.c <- rbind(sgd.auc.c,round(AUCBoot(pred, tPheno[ind.test]),4)[[3]])
	sgd.n <- rbind(sgd.n,length(idx))
	sgd.nr2 <- rbind(sgd.nr2,round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4))
}

idx = which.max(sgd.nr2)
rflr.sgd.auc.a <- sgd.auc.a[idx]
rflr.sgd.auc.b <- sgd.auc.b[idx]
rflr.sgd.auc.c <- sgd.auc.c[idx] 
rflr.sgd.nr2 <- sgd.nr2[idx] 
rflr.sgd.n <- sgd.n[idx]

idx = which(rfImpnr == max(rfImpnr))
rfimp.auc.a <- rfImpauc1[idx]
rfimp.auc.b <- rfImpauc2[idx]
rfimp.auc.c <- rfImpauc3[idx]
rfimp.nr2 <- max(rfImpnr)
rfimp.n <- rfImpn[idx]

idx = which(rflrnr == max(rflrnr))
rflr.auc.a <- rflrauc1[idx]
rflr.auc.b <- rflrauc2[idx]
rflr.auc.c <- rflrauc3[idx]
rflr.nr2 <- max(rflrnr)
rflr.n <- rflrn[idx]

imp.i <- i

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RF-LR implemented in bipolar disorder paper (li-chung chuang)
print("Running OG RFLR")
rf4 <- ranger(x=x[ind.train,],seed = 0, y=as.factor(tPheno[ind.train]),probability=TRUE,classification=TRUE,num.trees=500,mtry=sqrt(dim(x)[2])*3,importance="permutation")
imp1 <- rf2$variable.importance[rf2$variable.importance>summary(rf2$variable.importance)[3]]
imp2 <- rf4$variable.importance[rf4$variable.importance>summary(rf4$variable.importance)[3]]
imp <- intersect(names(imp1),names(imp2))
imp.idx <- match( imp,obj.SNP$map[,2])
x <- G[]
colnames(x) <- obj.SNP$map[,2]

logit <- big_spLogReg(X = G, y01.train = tPheno[ind.train],alphas = c(1,0.5,0.05,0.001),ind.train = ind.train,ind.col=imp.idx)
pred1 <- predict(logit, X = G, ind.row = ind.test ,ind.col=imp.idx)

rflr.og.auc.a <- round(AUCBoot(pred1,tPheno[ind.test])[[1]],4)
rflr.og.auc.b <- round(AUCBoot(pred1,tPheno[ind.test])[[2]],4)
rflr.og.auc.c <- round(AUCBoot(pred1,tPheno[ind.test])[[3]],4)
rflr.og.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4)
rflr.og.n <- summary(logit)$nb_var

init = summary(logit)$beta[1][[1]]
names(init) <- obj.SNP$map[,2][imp.idx]
fit_gd = sgd(init,X = x[ind.train,imp.idx],y = tPheno[ind.train],stepsize=1e-3,type="adam")
pred <- big_prodVec(G, fit_gd$par, ind.row = ind.test, ind.col = imp.idx)
rflr.og.sgd.auc.a <- round(AUCBoot(pred, tPheno[ind.test]),4)[[1]]
rflr.og.sgd.auc.b <- round(AUCBoot(pred, tPheno[ind.test]),4)[[2]]
rflr.og.sgd.auc.c <- round(AUCBoot(pred, tPheno[ind.test]),4)[[3]]
rflr.og.sgd.n <- length(imp)
rflr.og.sgd.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RF-RF
print("running RF-RF")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]
cut <- seq(0.201,1.01,0.1)

rfrfnr2 <- NULL
rfrfauc.a <- NULL
rfrfauc.b <- NULL
rfrfauc.c <- NULL
rfrfn <- NULL
rfPred <- NULL

for (j in seq(0.201,1.01,0.1)){
	print(j)
	cutSnps <- rownames(p[p[,2] < j,])
	hyper_grid <- expand.grid(
  		mtry = c(sqrt(length(cutSnps))/2,sqrt(length(cutSnps)),sqrt(length(cutSnps))*2,length(cutSnps)*0.1,length(cutSnps)*0.5),
  		node_size  = c(1,20,40),
  		num.trees = c(250, 500),
  		OOB_RMSE   = 0
	)
	cutSnpsRF <- cutSnps
	for(i in 1:nrow(hyper_grid)){
        	rf <- ranger(x=x[ind.train,colnames(x) %in% cutSnps],y=tPheno[ind.train],seed=876,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        	hyper_grid$OOB_RMSE[i] <- rf$prediction.error
	}
	i =  which.min(hyper_grid$OOB_RMSE)
	rf <- ranger(x=x[ind.train,colnames(x) %in% cutSnps],y=tPheno[ind.train],seed=876,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
	pred <- predict(rf,x[ind.test,colnames(x) %in% cutSnps])
	rfPred <- rbind(rfPred,pred)
	rfrfnr2 <- rbind(rfrfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4))
	rfrfauc.a <- rbind(rfrfauc.a,round(AUCBoot(pred$predictions[,1],tPheno[ind.test]),4)[[1]])
	rfrfauc.b <- rbind(rfrfauc.b,round(AUCBoot(pred$predictions[,1],tPheno[ind.test]),4)[[2]])
	rfrfauc.c <- rbind(rfrfauc.c,round(AUCBoot(pred$predictions[,1],tPheno[ind.test]),4)[[3]])
	rfrfn <- rbind(rfrfn,length(cutSnps))	
}

idx <- which.max(rfrfnr2)
rfrf.nr2 <- rfrfnr2[idx]
rfrf.auc.a <- rfrfauc.a[idx]
rfrf.auc.b <- rfrfauc.b[idx]
rfrf.auc.c <- rfrfauc.c[idx]
rfrf.n <- rfrfn[idx]

#################################################################################################################################################
# RFprior
# Input the betas/beta_se as a value of how often the snp should be selected [using split.select.weights]
print("Starting RFbeta")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]
weight <- (t - min(t))/ (max(t) - min(t))
p <- dim(G)[2]
hyper_grid <- expand.grid(
  mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,sqrt(p)*5,0.1 * p),
  node_size  = c(1,20,40),
  num.trees = c(250, 500),
  OOB_RMSE   = 0
)
print("normal method")
# normal method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),split.select.weights = weight,seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,])
rfNR2.2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.2a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.2b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.2c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

print("ic method")
# impurity_corrected method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,split.select.weights=weight,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,])
rfNR2.1 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.1a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.1b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.1c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

rfprior.auc.a <- c(rfAUC.1a,rfAUC.2a)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
rfprior.auc.b <- c(rfAUC.1b,rfAUC.2b)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
rfprior.auc.c <- c(rfAUC.1c,rfAUC.2c)[which(c(rfNR2.1 ,rfNR2.2) == max(rfNR2.1 ,rfNR2.2))]
rfprior.nr2 <- max(rfNR2.1 ,rfNR2.2 )
rfprior.n <- p

print("Finished running")
#################################################################################################################################################
# Output results

write.table(cbind(modelN,causalN,seed,split,strength,i.SCT.nr2,i.CT.nr2,i.pca.ct.nr2,i.pca.ct.p.nr2,plr.nr2,sis.lasso.nr2,sis.lasso.sgd.nr2,sis.rf.nr2,sis.rfimp.nr2,sis.rfprior.nr2,v.rf.nr2,vs.rf.nr2,rf.nr2,pca.rf.nr2,stack.rf.nr2,rfimp.nr2,rflr.nr2,rflr.sgd.nr2,rflr.og.nr2,rflr.og.sgd.nr2,rfrf.nr2,rfprior.nr2,rf.c.idx,rf.p,imp.i),file="/home3/DrYao/FINALest/Simulation/NR2Indep",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/NR2Indep"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,i.SCT.auc.a,i.CT.auc.a,i.pca.ct.auc.a,i.pca.ct.p.auc.a,plr.auc.a,sis.lasso.auc.a,sis.lasso.sgd.auc.a,sis.rf.auc.a,sis.rfimp.auc.a,sis.rfprior.auc.a,v.rf.auc.a,vs.rf.auc.a,rf.auc.a,pca.rf.auc.a,stack.rf.auc.a,rfimp.auc.a,rflr.auc.a,rflr.sgd.auc.a,rflr.og.auc.a,rflr.og.sgd.auc.a,rfrf.auc.a,rfprior.auc.a),file="/home3/DrYao/FINALest/Simulation/AUCs.aIndep",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.aIndep"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,i.SCT.n,i.CT.n,i.pca.ct.n,i.pca.ct.p.n,plr.n,sis.lasso.n,sis.lasso.sgd.n,sis.rf.n,sis.rfimp.n,sis.rfprior.n,v.rf.n,vs.rf.n,rf.n,pca.rf.n,stack.rf.n,rfimp.n,rflr.n,rflr.sgd.n,rflr.og.n,rflr.og.sgd.n,rfrf.n,rfprior.n),file="/home3/DrYao/FINALest/Simulation/NsIndep",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/NsIndep"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,i.SCT.auc.b,i.CT.auc.b,i.pca.ct.auc.b,i.pca.ct.p.auc.b,plr.auc.b,sis.lasso.auc.b,sis.lasso.sgd.auc.b,sis.rf.auc.b,sis.rfimp.auc.b,sis.rfprior.auc.b,v.rf.auc.b,vs.rf.auc.b,rf.auc.b,pca.rf.auc.b,stack.rf.auc.b,rfimp.auc.b,rflr.auc.b,rflr.sgd.auc.b,rflr.og.auc.b,rflr.og.sgd.auc.b,rfrf.auc.b,rfprior.auc.b),file="/home3/DrYao/FINALest/Simulation/AUCs.bIndep",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.bIndep"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,i.SCT.auc.c,i.CT.auc.c,i.pca.ct.auc.c,i.pca.ct.p.auc.c,plr.auc.c,sis.lasso.auc.c,sis.lasso.sgd.auc.c,sis.rf.auc.c,sis.rfimp.auc.c,sis.rfprior.auc.c,v.rf.auc.c,vs.rf.auc.c,rf.auc.c,pca.rf.auc.c,stack.rf.auc.c,rfimp.auc.c,rflr.auc.c,rflr.sgd.auc.c,rflr.og.auc.c,rflr.og.sgd.auc.c,rfrf.auc.c,rfprior.auc.c),file="/home3/DrYao/FINALest/Simulation/AUCs.cIndep",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.cIndep"),append=T)

# Remove tmp files
system("cd ../")
system(paste("rm -r logs",paste(seed,modelN,sep="_"),sep=""))
system("rm -r tmp-data")
file.remove(paste0(tmpfile, ".rds"))


