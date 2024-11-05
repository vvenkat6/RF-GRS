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
        grid.lpS.thr = 0.9999 * seq_log(max(0.1, min(p.val, na.rm = TRUE)),max(p.val, na.rm = TRUE), 7)
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
rf.nr22 <- rfnr2[id]
rf.auc.a2 <- rfauc1[id]
rf.auc.b2 <- rfauc2[id]
rf.auc.c2 <- rfauc3[id]
rf.n2 <- rfn[id]

rf.c.idx2 <- rfc[id]
rf.p2 <- rfp[id]
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
        grid.lpS.thr = c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
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
rf.nr23 <- rfnr2[id]
rf.auc.a3 <- rfauc1[id]
rf.auc.b3 <- rfauc2[id]
rf.auc.c3 <- rfauc3[id]
rf.n3 <- rfn[id]

rf.c.idx3 <- rfc[id]
rf.p3 <- rfp[id]
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
                                mtry = c(sqrt(p)),
                                num.trees = c(500),
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
rf.nr24 <- rfnr2[id]
rf.auc.a4 <- rfauc1[id]
rf.auc.b4 <- rfauc2[id]
rf.auc.c4 <- rfauc3[id]
rf.n3 <- rfn[id]

rf.c.idx4 <- rfc[id]
rf.p4 <- rfp[id]


print("Finished running")
#################################################################################################################################################
# Output results

write.table(cbind(modelN,causalN,seed,split,strength,rf.nr2,rf.nr22,rf.nr23,rf.nr24),file="/home3/DrYao/FINALest/Simulation/NR2rfs",quote=FALSE,row.names=FALSE,col.names= !file.exists("/home3/DrYao/FINALest/Simulation/NR2rfs"),append=T)

# Remove tmp files
system("cd ../")
system(paste("rm -r logs",paste(seed,modelN,sep="_"),sep=""))
system("rm -r tmp-data")
file.remove(paste0(tmpfile, ".rds"))


