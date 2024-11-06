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
suppressPackageStartupMessages(library(tidyverse))

# Load Data
tmpfile <- tempfile()
bedfile = "/home3/DrYao/PATCH_LIGHTs/LIGHTS/May2022/dataQC/FinalLightsData.bed"
snp_readBed(bedfile, backingfile = tmpfile)

obj.SNP <- snp_attach(paste0(tmpfile, ".rds"))
G   <- obj.SNP$genotypes
CHR <- obj.SNP$map$chromosome
POS <- obj.SNP$map$physical.pos

# Load indices
ind.train <- get(load("/home3/DrYao/FINALest/Atopy/splits/ind.train1"))
ind.test <- get(load("/home3/DrYao/FINALest/Atopy/splits/ind.test1"))

# Covariate file
cov <- read.table("/home3/DrYao/PATCH_LIGHTs/LIGHTS/May2022/dataQC/Cov.atopy",header=TRUE,sep="\t",stringsAsFactors=FALSE)
Cov <- cov[cov$Test_ID %in% obj.SNP$fam[,1],]
pca <- read.table("/home3/DrYao/PATCH_LIGHTs/LIGHTS/May2022/dataQC/PCA",header=TRUE,sep=" ",stringsAsFactors=FALSE)
rownames(pca) <- pca$V1
pca <- pca[,3:22]

TrainCov <- Cov[ind.train,]
TestCov <- Cov[ind.test,]
AtopyTrainCov <- cbind(as.matrix(pca[ind.train,1:10]),TrainCov$age,TrainCov$sex)
AtopyTestCov <- cbind(as.matrix(pca[ind.test,1:10]),TestCov$age,TestCov$sex)
tPheno <- Cov$pheno
colnames(AtopyTrainCov) <- c(colnames(pca[,1:10]),"age","sex")
colnames(AtopyTestCov) <- c(colnames(pca[,1:10]),"age","sex")
NCORES = 4
#########################################################################################################################
# C+T prep
print("Prepping for C+T")
gwas.train <- big_univLogReg(G,tPheno[ind.train],ind.train = ind.train,covar.train = AtopyTrainCov, ncores=NCORES)
gwas_gc <- snp_gc(gwas.train)

thrs <- sort(c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100))))
lpS <- -predict(gwas_gc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stacked C+T
print("Stacked C+T starting...")
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS, exclude = which(is.na(lpS)),ncores = NCORES)
multi_PRS <- snp_grid_PRS(G, all_keep, gwas_gc$estim, lpS = lpS, ind.row = ind.train,ncores=NCORES)
final_mod <- snp_grid_stacking(multi_PRS, tPheno[ind.train], ncores = NCORES, K = 4,covar.train=AtopyTrainCov,pf.covar=(rep(0,12)))
new_beta <- final_mod$beta.G
ind <- which(new_beta != 0)
pred <- final_mod$intercept + big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind) + (as.matrix(AtopyTestCov) %*% final_mod$beta.covar)

SCT.auc.a <- round(AUCBoot(pred, tPheno[ind.test]),4)[[1]]
SCT.auc.b <- round(AUCBoot(pred, tPheno[ind.test]),4)[[2]]
SCT.auc.c <- round(AUCBoot(pred, tPheno[ind.test]),4)[[3]]
SCT.n <- length(ind)
SCT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)
print(SCT.nr2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Traditional C+T
print("Starting Traditional C+T")
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
  tidyr::unnest(cols = "thr.lp")
s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
        single_PRS <- rowSums(X[, ind + s * (0:21)])  
        bigstatsr::AUC(single_PRS, y.train)
        }, ind = 1:s, s = s, y.train = tPheno[ind.train],a.combine = 'c', block.size = 1, ncores = NCORES)

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
pred <- snp_PRS(G, gwas_gc$estim[ind.keep], ind.test = ind.test, ind.keep = ind.keep,lpS.keep = lpS[ind.keep], thr.list = max_prs$thr.lp)
round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)
# Calculating coefficients for covariates
dt <- data.frame(y=tPheno[ind.train],AtopyTrainCov)
lmfit <- lm(y~., data=dt)
pred <- pred + AtopyTestCov %*% lmfit$coefficients[2:13]

CT.auc.a <- round(AUCBoot(pred,tPheno[ind.test]),4)[[1]]
CT.auc.b <- round(AUCBoot(pred,tPheno[ind.test]),4)[[2]]
CT.auc.c <- round(AUCBoot(pred,tPheno[ind.test]),4)[[3]]
CT.n <- length(ind.keep)
CT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)
print(CT.nr2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA PRS
print("Starting Traditional PRS")

multi_PRS <- snp_grid_PRS(G, all_keep, gwas_gc$estim, lpS = lpS, ind.row = ind.test,ncores=NCORES)
multi_PRS <- sweep(multi_PRS[],1, (AtopyTestCov %*% lmfit$coefficients[2:13]),"+")

corr_matrix <- cor(t(multi_PRS))
pca.prs <- princomp(corr_matrix)$loadings[,1]
pca.ct.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
pca.ct.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
pca.ct.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
pca.ct.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
pca.ct.n <- dim(G)[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Penalized Linear Regression
print("Starting PLR")
# load subset values
logit <- big_spLogReg(X = G, y01.train = tPheno[ind.train],alphas = c(1,0.5,0.05,0.001),ind.train = ind.train,ncores=NCORES,covar.train=AtopyTrainCov)
pred1 <- predict(logit, X = G, ind.row = ind.test,covar.row=AtopyTestCov)

plr.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4)
plr.auc.a <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[1]]
plr.auc.b <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[2]]
plr.auc.c <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[3]]
plr.n <- summary(logit)$nb_var[which(summary(logit)$validation_loss == min(summary(logit)$validation_loss))]

print("Writing Results")
write.table(cbind("Split1",SCT.auc.a,SCT.auc.b,SCT.auc.c,SCT.n,SCT.nr2,CT.auc.a,CT.auc.b,CT.auc.c,CT.n,CT.nr2,pca.ct.auc.a,pca.ct.auc.b,pca.ct.auc.c,pca.ct.nr2,pca.ct.n,plr.auc.a,plr.auc.b,plr.auc.c,plr.nr2,plr.n),file = "CT-PLRouts",quote=FALSE,row.names=FALSE,col.names= !file.exists("CT-PLRouts"),append=T)

file.remove(paste0(tmpfile, ".rds"))



