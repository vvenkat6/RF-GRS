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
Sys.sleep(seed * 5 * modelN * strength)
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
lpS = -log10(sumstats$P)
NCORES = 4
gwas_gc <- sumstats
gwas_gc$lpS <- -log10(sumstats$P)
rownames(gwas_gc) <- gwas_gc$SNP
#################################################################################################################################################
# STacked C+T with summary statistics

# computes sets of variants resulting from the clumping procedure that is applied repeatedly with different values of hyper-parameters
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = -log10(sumstats$P), exclude = which(is.na(-log10(sumstats$P))),ncores = NCORES)

# Computes C+T scores for all the SNP sets in all_keep
multi_PRS <- snp_grid_PRS(G, all_keep, sumstats$BETA, lpS = -log10(sumstats$P), ind.row = ind.train,ncores=NCORES)

# Penalized regression used to learn optimal combination of C+T scores
final_mod <- snp_grid_stacking(multi_PRS, tPheno[ind.train], ncores = NCORES, K = 4)
new_beta <- final_mod$beta.G
ind <- which(new_beta != 0)
pred <- final_mod$intercept + big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)
SCT.auc.a <- round(AUCBoot(pred, tPheno[ind.test]),4)[[1]]
SCT.auc.b <- round(AUCBoot(pred, tPheno[ind.test]),4)[[2]]
SCT.auc.c <- round(AUCBoot(pred, tPheno[ind.test]),4)[[3]]
SCT.n <- length(ind)
SCT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#################################################################################################################################################
# Traditional C+T (ss) [also known as maxCT in Prive paper]
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
pred <- snp_PRS(G, sumstats$BETA[ind.keep], ind.test = ind.test, ind.keep = ind.keep,lpS.keep = sumstats$lpval[ind.keep], thr.list = max_prs$thr.lp)
ss.CT.auc.a <- round(AUCBoot(pred,tPheno[ind.test]),4)[[1]]
ss.CT.auc.b <- round(AUCBoot(pred,tPheno[ind.test]),4)[[2]]
ss.CT.auc.c <- round(AUCBoot(pred,tPheno[ind.test]),4)[[3]]
ss.CT.n <- length(ind.keep)
ss.CT.nr2 <- round(NagelkerkeR2(glm(as.factor(tPheno[ind.test])~pred,family="binomial"))$R2,4)

#################################################################################################################################################
# PCA-PRS with summary statistics: R2 and p-value
multi_PRS <- snp_grid_PRS(G, all_keep, sumstats$BETA, lpS = -log10(sumstats$P), ind.row = ind.test,ncores=NCORES)
corr_matrix <- cor(t(cbind(multi_PRS[])))
pca.prs <- princomp(corr_matrix)$loadings[,1]
ss.pca.ct.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
ss.pca.ct.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
ss.pca.ct.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
ss.pca.ct.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
ss.pca.ct.n <- dim(G)[2]

#################################################################################################################################################
# PCA-PRS with summary statistics: Only p-value
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = -log10(sumstats$P), exclude = which(is.na(-log10(sumstats$P))),ncores = NCORES,grid.thr.r2 = c(0.1),grid.base.size=c(50))
multi_PRS <- snp_grid_PRS(G, all_keep, sumstats$BETA, lpS = -log10(sumstats$P), ind.row = ind.test,ncores=NCORES)
corr_matrix <- cor(t(cbind(multi_PRS[])))
pca.prs <- princomp(corr_matrix)$loadings[,1]
ss.pca.ct.p.auc.a <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[1]]
ss.pca.ct.p.auc.b <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[2]]
ss.pca.ct.p.auc.c <- round(AUCBoot(pca.prs,tPheno[ind.test]),4)[[3]]
ss.pca.ct.p.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pca.prs,family = "binomial"))$R2,4)
ss.pca.ct.p.n <- dim(G)[2]

#################################################################################################################################################
# Lassosum/Ldpred prep
print("Starting Lassosum/Ldpred prep")
dat <- as.data.frame(cbind(CHR,POS))
colnames(dat) <- c("chr","pos")
POS.hg19 <- snp_modifyBuild(info_snp=dat, liftOver = "/home3/vvenka23/FNC-Lasso/liftover/liftOver",from = "hg38",to = "hg19")$pos
#remove.idx <- RSID[which(is.na(POS.hg19))]
system("mkdir tmp-data")
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")),add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
POS2 <- snp_asGeneticPos(CHR, POS.hg19, dir = "tmp-data", ncores = 4)
chr = 21
ind.chr <- which(sumstats$CHR == chr)
ord <- order(POS2[ind.chr])
corr0 <- snp_cor(
    G,
    ind.col = ind.chr[ord],
    ncores = 4,
    infos.pos = POS2[ind.chr[ord]],
    size = 3 / 1000
  )

ld <- Matrix::colSums(corr0^2)
corr <- as_SFBM(corr0)
df_beta <- sumstats[ind.chr[ord], ]

df_beta[is.na(df_beta)] <- 0
colnames(df_beta) <- c("CHR","SNP","BP","A1","TEST","n_eff","beta","beta_se","L50","U50","STAT","P")
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2,
                                 sample_size = 10654, blocks = NULL)))
ldsc_h2_est <- ldsc[["h2"]]
if (ldsc_h2_est <0){ldsc_h2_est = 0.01}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lassosum (SS)

beta_lassosum2 <- snp_lassosum2(corr, df_beta, nlambda = 10, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")
pred_grid2 <- big_prodMat(G, beta_lassosum2)
params2$score <- apply(pred_grid2[ind.train, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
    summary(glm(tPheno[ind.train] ~ x, family = "binomial"))$coef["x", 3]
  })
best_grid_lassosum2 <- params2 %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_lassosum2[, .]
pred_lassosum <- big_prodVec(G, best_grid_lassosum2, ind.row = ind.test)

lassosum.auc.a <- round(AUCBoot(pred_lassosum,tPheno[ind.test]),4)[[1]]
lassosum.auc.b <- round(AUCBoot(pred_lassosum,tPheno[ind.test]),4)[[2]]
lassosum.auc.c <- round(AUCBoot(pred_lassosum,tPheno[ind.test]),4)[[3]]
lassosum.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred_lassosum,family = "binomial"))$R2,4)
lassosum.n <- sum(best_grid_lassosum2!=0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LdPred (SS)
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(h2_seq <- round(ldsc_h2_est * c(0.3, 0.7, 1, 1.4), 4))
grid.param <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid)
pred_grid[is.na(pred_grid)] <- 0
grid.param$score <- big_univLogReg(as_FBM(pred_grid[ind.train, ]), tPheno[ind.train])$score
grid.param %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(desc(score)) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)
best_grid_sp <- grid.param %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]
pred_sp <- big_prodVec(G, best_grid_sp, ind.row = ind.test)

ld.auc.a <- round(AUCBoot(pred_sp,tPheno[ind.test]),4)[[1]]
ld.auc.b <- round(AUCBoot(pred_sp,tPheno[ind.test]),4)[[2]]
ld.auc.c <- round(AUCBoot(pred_sp,tPheno[ind.test]),4)[[3]]
ld.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred_sp,family = "binomial"))$R2,4)
ld.n <- sum(best_grid_sp!=0)

#################################################################################################################################################
# SIS + Lasso (SS)
t <- sumstats$BETA/sumstats$SE
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
# RFProb (ss) and RFLR (ss)

weight <- (t - min(t) )/ (max(t) - min(t))
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
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),split.select.weights = weight,seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
	print(i)
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,split.select.weights = weight,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
pred <- predict(rf,x[ind.test,])
rfic.NR2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfic.AUC.a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfic.AUC.b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfic.AUC.c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)
rfic.n <- p

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# importance values
print("Running RFLR")
rf2 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=1,verbose=FALSE,probability=TRUE,classification=TRUE,split.select.weights = weight,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
rf3 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=2,verbose=FALSE,probability=TRUE,classification=TRUE,split.select.weights = weight,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
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
}

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
rf4 <- ranger(x=x[ind.train,],seed = 0,split.select.weights = weight, y=as.factor(tPheno[ind.train]),probability=TRUE,classification=TRUE,num.trees=500,mtry=sqrt(dim(x)[2])*3,importance="permutation")
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

#################################################################################################################################################
# RFprior : C+T
# RFs for different CT threshold Prep with p-values
print("CT threshold prep")
x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]
weight <- (t - min(t) )/ (max(t) - min(t))

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
                                        rf <- ranger(x=x[ind.train,names.idx],y=tPheno[ind.train],seed=0,split.select.weights = weight[names.idx],verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
                                        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
                                        }
                                else{ hyper_grid$OOB_RMSE[i] <- 100 }
                                }
                        i =  which.min(hyper_grid$OOB_RMSE)
                        rfi <- rbind(rfi,i)
                        rf <- ranger(x=x[ind.train,names.idx],y=tPheno[ind.train],seed=0,verbose=FALSE,probability=TRUE,split.select.weights = weight[names.idx],classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i])
                        pred.Train <- predict(rf,x[ind.train,names.idx])$predictions[,2]
                        pred.Test <- predict(rf,x[ind.test,names.idx])$predictions[,2]
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

id = which.max(rfnr2)
rfprior.nr2 <- rfnr2[id]
rfprior.auc.a <- rfauc1[id]
rfprior.auc.b <- rfauc2[id]
rfprior.auc.c <- rfauc3[id]
rfprior.n <- rfn[id]


print("Finished running")
#################################################################################################################################################
# Output results

write.table(cbind(modelN,causalN,seed,split,strength,ldsc_h2_est,SCT.nr2,ss.CT.nr2,ss.pca.ct.nr2,ss.pca.ct.p.nr2,lassosum.nr2,ld.nr2,sis.lasso.nr2,sis.rf.nr2,sis.rfimp.nr2,sis.rfprior.nr2,v.rf.nr2,vs.rf.nr2,rf.nr2,pca.rf.nr2,stack.rf.nr2,rfic.NR2,rfimp.nr2,rflr.nr2,rflr.og.nr2,rfprior.nr2),file="/home3/DrYao/FINALest/Simulation/NR2SS",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/NR2SS"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,ldsc_h2_est,SCT.auc.a,ss.CT.auc.a,ss.pca.ct.auc.a,ss.pca.ct.p.auc.a,lassosum.auc.a,ld.auc.a,sis.lasso.auc.a,sis.rf.auc.a,sis.rfimp.auc.a,sis.rfprior.auc.a,v.rf.auc.a,vs.rf.auc.a,rf.auc.a,pca.rf.auc.a,stack.rf.auc.a,rfic.AUC.a,rfimp.auc.a,rflr.auc.a,rflr.og.auc.a,rfprior.auc.a),file="/home3/DrYao/FINALest/Simulation/AUCs.aSS",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.aSS"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,ldsc_h2_est,SCT.n,ss.CT.n,ss.pca.ct.n,ss.pca.ct.p.n,lassosum.n,ld.n,sis.lasso.n,sis.rf.n,sis.rfimp.n,sis.rfprior.n,v.rf.n,vs.rf.n,rf.n,pca.rf.n,stack.rf.n,rfic.n,rfimp.n,rflr.n,rflr.og.n,rfprior.n),file="/home3/DrYao/FINALest/Simulation/NsSS",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/NsSS"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,ldsc_h2_est,SCT.auc.b,ss.CT.auc.b,ss.pca.ct.auc.b,ss.pca.ct.p.auc.b,lassosum.auc.b,ld.auc.b,sis.lasso.auc.b,sis.rf.auc.b,sis.rfimp.auc.b,sis.rfprior.auc.b,v.rf.auc.b,vs.rf.auc.b,rf.auc.b,pca.rf.auc.b,stack.rf.auc.b,rfic.AUC.b,rfimp.auc.b,rflr.auc.b,rflr.og.auc.b,rfprior.auc.b),file="/home3/DrYao/FINALest/Simulation/AUCs.bSS",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.bSS"),append=T)

write.table(cbind(modelN,causalN,seed,split,strength,ldsc_h2_est,SCT.auc.c,ss.CT.auc.c,ss.pca.ct.auc.c,ss.pca.ct.p.auc.c,lassosum.auc.c,ld.auc.c,sis.lasso.auc.c,sis.rf.auc.c,sis.rfimp.auc.c,sis.rfprior.auc.c,v.rf.auc.c,vs.rf.auc.c,rf.auc.c,pca.rf.auc.c,stack.rf.auc.c,rfic.AUC.c,rfimp.auc.c,rflr.auc.c,rflr.og.auc.c,rfprior.auc.c),file="/home3/DrYao/FINALest/Simulation/AUCs.cSS",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/AUCs.cSS"),append=T)

# Remove tmp files
system("cd ../")
system(paste("rm -r logs",paste(seed,modelN,sep="_"),sep=""))
system("rm -r tmp-data")
file.remove(paste0(tmpfile, ".rds"))


