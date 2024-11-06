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
bedfile = "/home3/DrYao/projects/LIGHTS/May2022/dataQC/FinalLightsData.bed"
snp_readBed(bedfile, backingfile = tmpfile)

obj.SNP <- snp_attach(paste0(tmpfile, ".rds"))
G   <- obj.SNP$genotypes
CHR <- obj.SNP$map$chromosome
POS <- obj.SNP$map$physical.pos

# Load indices
ind.train <- get(load("/home3/DrYao/FINALest/Atopy/splits/ind.train1"))
ind.test <- get(load("/home3/DrYao/FINALest/Atopy/splits/ind.test1"))
split = "Split1"

# Covariate file
cov <- read.table("/home3/DrYao/projects/LIGHTS/May2022/dataQC/Cov.atopy",header=TRUE,sep="\t",stringsAsFactors=FALSE)
Cov <- cov[cov$Test_ID %in% obj.SNP$fam[,1],]
pca <- read.table("/home3/DrYao/projects/LIGHTS/May2022/dataQC/PCA",header=TRUE,sep=" ",stringsAsFactors=FALSE)
rownames(pca) <- pca$V1
pca <- pca[,3:22]

TrainCov <- Cov[ind.train,]
TestCov <- Cov[ind.test,]
AtopyCov <- cbind(as.matrix(pca[,1:10]),Cov$age,Cov$sex)
colnames(AtopyCov) <- c(colnames(pca[,1:10]),"age","sex")
AtopyTrainCov <- cbind(as.matrix(pca[ind.train,1:10]),TrainCov$age,TrainCov$sex)
AtopyTestCov <- cbind(as.matrix(pca[ind.test,1:10]),TestCov$age,TestCov$sex)
tPheno <- Cov$pheno
colnames(AtopyTrainCov) <- c(colnames(pca[,1:10]),"age","sex")
colnames(AtopyTestCov) <- c(colnames(pca[,1:10]),"age","sex")
NCORES = 8
#########################################################################################################################
# Clumping prep for RF

print("Prepping for C+T")
gwas.train <- big_univLogReg(G,tPheno[ind.train],ind.train = ind.train,covar.train = AtopyTrainCov, ncores=NCORES)
gwas_gc <- snp_gc(gwas.train)
p.value.all <- predict(gwas.train, log10 = FALSE)
thrs <- sort(c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100))))
lpS <- -predict(gwas_gc)

grid.thr.r2 = c(0.001,0.01, 0.1, 0.5)
i = 1
all_keep = list()

for (val in grid.thr.r2){
        all_keep[[i]] = snp_clumping(G, infos.chr = CHR,ind.row = ind.train,thr.r2 = val, S=abs(gwas_gc$score),size = 250,infos.pos = POS,ncores=NCORES)
        i = i + 1
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Main RF

x <- G[]
colnames(x) <- obj.SNP$map[,2]
rownames(x) <- obj.SNP$fam[,1]
x <- cbind(x,AtopyCov)

rfnr2 <- NULL
rfauc1 <- NULL
rfauc2 <- NULL
rfauc3 <- NULL
rfn <- NULL
count = 0
pred.Train <- NULL
pred.Test <- NULL

rfnr2.1 <- NULL
rfauc1.1 <- NULL
rfauc2.1 <- NULL
rfauc3.1 <- NULL
rfn.1 <- NULL

rfnr2.2 <- NULL
rfauc1.2 <- NULL
rfauc2.2 <- NULL
rfauc3.2 <- NULL
rfn.2 <- NULL

for (idx in all_keep){
        print("All keep index:")
        print(count)
	count = count + 1
        p <- length(idx)
	p.val <- lpS[idx]
	grid.lpS.thr = 0.9999 * seq_log(max(0.1, min(p.val, na.rm = TRUE)),max(p.val, na.rm = TRUE), 5)
	for (pcut in grid.lpS.thr){
                a <- gwas_gc[idx,]
                names <- rownames(a[p.val >= pcut,])
		p <- length(names)
		if (p>3){
        		hyper_grid <- expand.grid(
                	mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,sqrt(p)*5,0.1 * p,0.5 * p),
                	node_size = c(100,300),
                	num.trees = c(2000,5000),
                	OOB_RMSE   = 0
                	)
			idx.fil <- c(names,which(colnames(x) =="PC1"),which(colnames(x) =="PC2"),which(colnames(x) =="PC3"),which(colnames(x) =="PC4"),which(colnames(x) =="PC5"),which(colnames(x) =="PC6"),which(colnames(x) =="PC7"),which(colnames(x) =="PC8"),which(colnames(x) =="PC9"),which(colnames(x) =="PC10"),which(colnames(x) =="age"),which(colnames(x) =="sex"))
        		a <- table(as.factor(tPheno[ind.train]))[[1]]/(table(as.factor(tPheno[ind.train]))[[1]] + table(as.factor(tPheno[ind.train]))[[2]])
        		b <- table(as.factor(tPheno[ind.train]))[[2]]/(table(as.factor(tPheno[ind.train]))[[1]] + table(as.factor(tPheno[ind.train]))[[2]])
			for(i in 1:nrow(hyper_grid)){
				if(hyper_grid$mtry[i] < length(names)){
                			rf <- ranger(x=x[ind.train,as.numeric(idx.fil)],y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b),always.split.variables=colnames(AtopyCov),importance="impurity_corrected")
                			hyper_grid$OOB_RMSE[i] <- rf$prediction.error
                		} else{ hyper_grid$OOB_RMSE[i] <- 100 }
			}
        	i =  which.min(hyper_grid$OOB_RMSE)
        	rf <- ranger(x=x[ind.train,as.numeric(idx.fil)],y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b),always.split.variables=colnames(AtopyCov),importance="impurity_corrected")
        	pred.Train <- rbind(pred.Train,predict(rf,x[ind.train,as.numeric(idx.fil)])$predictions[,2])
        	pred.Test <- rbind(pred.Test,predict(rf,x[ind.test,as.numeric(idx.fil)])$predictions[,2])
        	pred <- predict(rf,x[ind.test,as.numeric(idx.fil)])$predictions[,2]
        	rfnr2 <- rbind(rfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred,family = "binomial"))$R2, 4))
        	rfauc1 <- rbind(rfauc1,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
        	rfauc2 <- rbind(rfauc2,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
        	rfauc3 <- rbind(rfauc3,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
        	rfn <- rbind(rfn,p)
        	print(rfnr2)
		}
		if(count == 3){
			# calculation for trf
			rfnr2.1 <- rbind(rfnr2.1,round(NagelkerkeR2(glm(tPheno[ind.test]~pred,family = "binomial"))$R2, 4))
                	rfauc1.1 <- rbind(rfauc1.1,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                	rfauc2.1 <- rbind(rfauc2.1,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                	rfauc3.1 <- rbind(rfauc3.1,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                	rfn.1 <- rbind(rfn.1,p)	
		}
		if(pcut == grid.lpS.thr[length(grid.lpS.thr)]){
			#calculation for crf
                        rfnr2.2 <- rbind(rfnr2.2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred,family = "binomial"))$R2, 4))
                        rfauc1.2 <- rbind(rfauc1.2,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                        rfauc2.2 <- rbind(rfauc2.2,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                        rfauc3.2 <- rbind(rfauc3.2,round(AUCBoot(pred,tPheno[ind.test])[[1]],4))
                        rfn.2 <- rbind(rfn.2,p)
                }
		
	}
}

idx2 <- which.max(rfnr2)
ctrf.nr2 <- rfnr2[idx2]
ctrf.auc.a <- rfauc1[idx2]
ctrf.auc.b <- rfauc2[idx2]
ctrf.auc.c <- rfauc3[idx2]
ctrf.n <- rfn[idx2]

idx2 <- which.max(rfnr2.1)
trf.nr2 <- rfnr2.1[idx2]
trf.auc.a <- rfauc1.1[idx2]
trf.auc.b <- rfauc2.1[idx2]
trf.auc.c <- rfauc3.1[idx2]
trf.n <- rfn.1[idx2]

idx2 <- which.max(rfnr2.2)
crf.nr2 <- rfnr2.2[idx2]
crf.auc.a <- rfauc1.2[idx2]
crf.auc.b <- rfauc2.2[idx2]
crf.auc.c <- rfauc3.2[idx2]
crf.n <- rfn.2[idx2]

print(ctrf.nr2,trf.nr2,crf.nr2)
# Stacked CTRF
logit <- big_spLogReg(X = as_FBM(t(pred.Train)), y01.train = tPheno[ind.train],alphas=c(1,0.5,0.05,0.001),ncores=NCORES)
pred1 <- predict(logit,X=as_FBM(t(pred.Test)))
stack.rf.auc.a <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[1]]
stack.rf.auc.b <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[2]]
stack.rf.auc.c <- round(AUCBoot(pred1,tPheno[ind.test]),4)[[3]]
stack.rf.nr2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred1,family = "binomial"))$R2, 4)
stack.rf.n <- dim(x)[2]
print(ctrf.nr2,trf.nr2,crf.nr2,stack.rf.nr2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RFLR
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected")
rf2 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=1,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected")
rf3 <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=2,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected")
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

        logit <- big_spLogReg(X = G, y01.train = tPheno[ind.train],alphas = c(1,0.5,0.05,0.001),ind.train = ind.train,ind.col=idx[idx<=dim(G)[2]])
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Weighted RF
# RFprior
# Input the betas/beta_se as a value of how often the snp should be selected [using split.select.weights]
t <- lpS
print("Starting RFbeta")
weight <- (t - min(t))/ (max(t) - min(t))
weight <- c(weight,rep(1,12))
p <- dim(x)[2]
hyper_grid <- expand.grid(
	mtry = c(sqrt(p)/2,sqrt(p),sqrt(p)*2,sqrt(p)*5,0.1 * p,0.5 * p),
	node_size = c(100,300),
	num.trees = c(2000,5000),
	OOB_RMSE   = 0
)

print("normal method")
# normal method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),split.select.weights = weight,seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b))
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b))
pred <- predict(rf,x[ind.test,])
rfNR2.2 <- round(NagelkerkeR2(glm(tPheno[ind.test]~pred$predictions[,2],family = "binomial"))$R2, 4)
rfAUC.2a <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[1]],4)
rfAUC.2b <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[2]],4)
rfAUC.2c <- round(AUCBoot(pred$predictions[,2],tPheno[ind.test])[[3]],4)

print("ic method")
# impurity_corrected method
for(i in 1:nrow(hyper_grid)){
        rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,split.select.weights=weight,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b))
        hyper_grid$OOB_RMSE[i] <- rf$prediction.error
}
i =  which.min(hyper_grid$OOB_RMSE)
rf <- ranger(x=as.matrix(x[ind.train,]),y=as.factor(tPheno[ind.train]),seed=0,verbose=FALSE,probability=TRUE,classification=TRUE,importance="impurity_corrected",num.trees=hyper_grid$num.trees[i],mtry=hyper_grid$mtry[i],min.node.size=hyper_grid$node_size[i],class.weights=c(a,b))
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# output variables
write.table(cbind(split,ctrf.nr2,crf.nr2,crf.nr2,stack.rf.nr2,rfimp.nr2,rflr.nr2,rfprior.nr2),file = "RFAtopyOuts.nr2",quote=FALSE,row.names=FALSE,col.names= !file.exists("RFAtopyOuts.nr2"),append=T)
write.table(cbind(split,ctrf.auc.a,crf.auc.a,crf.auc.a,stack.rf.auc.a,rfimp.auc.a,rflr.auc.a,rfprior.auc.a ),file = "RFAtopyOuts.auc.a",quote=FALSE,row.names=FALSE,col.names= !file.exists("RFAtopyOuts.auc.a"),append=T)
write.table(cbind(split,ctrf.auc.b,crf.auc.b,crf.auc.b,stack.rf.auc.b,rfimp.auc.b,rflr.auc.b,rfprior.auc.b),file = "RFAtopyOuts.auc.b",quote=FALSE,row.names=FALSE,col.names= !file.exists("RFAtopyOuts.auc.b"),append=T)
write.table(cbind(split,ctrf.auc.c,crf.auc.c,crf.auc.c,stack.rf.auc.c,rfimp.auc.c,rflr.auc.c,rfprior.auc.c),file = "RFAtopyOuts.auc.c",quote=FALSE,row.names=FALSE,col.names= !file.exists("RFAtopyOuts.auc.c"),append=T)
write.table(cbind(split,ctrf.n,crf.n,crf.n,stack.rf.n,rfimp.n,rflr.n,rfprior.n),file = "RFAtopyOuts.n",quote=FALSE,row.names=FALSE,col.names= !file.exists("RFAtopyOuts.n"),append=T)





