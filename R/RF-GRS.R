#' Perform SNP Analysis
#'
#' This function calculates genetic risk score using random forest when there is no base data available
#'
#' @param bedfile Path to the BED file.
#' @param seed Seed for random number generation.
#' @param r2_grid R2 grid for clumping.
#' @param sumstats_file Path to the summary statistics file.
#' @return A list containing the analysis results.
#' @export

snp_analysis <- function(bedfile = NULL, seed = 123, r2_grid = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95), sumstats_file = NULL, NCORES = 4) {

  set.seed(seed)  

  if (is.null(bedfile)) {
  	print("Error: Please provide bedfile.")
	quit(status = 1)
  }

  tmpfile <- tempfile()
  snp_readBed(bedfile, backingfile = tmpfile)
  obj.SNP <- snp_attach(paste0(tmpfile, '.res'))
  G <- obj.SNP$genotypes
  CHR <- obj.SNP$map$chromosome
  POS <- obj.SNP$map$physical.pos

  if (!is.null(sumstats_file)) {
  	# Process summary statistics file
  	print("Processing summary statistics file...")

	if ("P" %in% colnames(sumstats_file)) {
      		print("Found 'P' column in summary statistics file")
      		# Process using 'P' column
		lpS <- -log10(sumstats_file$P)
    	} else if ("lpS" %in% colnames(sumstats_file)) {
      		print("Found 'lpS' column in summary statistics file")
    		lpS <- sumstats_file$lpS
    	} else {
      	log_print("Error: Summary statistics file must contain a column named 'P' for p-values or 'lpS' for -log10(P)")
      	quit(status = 1)
    	}
	train.ratio <- 0.8
        ind.train <- sample(nrow(G), round(train.ratio*nrow(G)))
        ind.test <- setdiff(rows_along(G), ind.train)

  } else {

  	print("No summary statistics file provided, splitting given data into training and test set to compute GRS...")

  	# Split data into train and validation set varying with seed
  	train.ratio <- 0.8
  	ind.train <- sample(nrow(G), round(train.ratio*nrow(G)))
  	ind.test <- setdiff(rows_along(G), ind.train)

  	# C+T prep
  	gwas.train <- big_univLogReg(G,tPheno[ind.train],ind.train = ind.train)
  	gwas_gc <- snp_gc(gwas.train)
  	p.value.all <- predict(gwas.train, log10 = FALSE)
  	thrs <- sort(c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100))))
  	lpS <- -predict(gwas_gc)

  	# p-value grid:
  	gwas_gc$lpS <- lpS
  	rownames(gwas_gc) <- obj.SNP$map[,2]
  }

  x <- G[]
  colnames(x) <- obj.SNP$map[,2]
  rownames(x) <- obj.SNP$fam[,1]

  pred.Train.cp <- data.frame(matrix(ncol = length(ind.train), nrow = 0))
  pred.Test.cp <- data.frame(matrix(ncol = length(ind.test), nrow = 0))
  rfnr2 <- NULL

  count = 0

  # Run clumping 
  all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train, lpS = lpS,grid.base.size = c(50, 100, 200),exclude = which(is.na(lpS)),ncores = NCORES)

  for (idx in all_keep){
    print('All_keep:')
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
          num.trees = c(5000),
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
        rfnr2 <- rbind(rfnr2,round(NagelkerkeR2(glm(tPheno[ind.test]~pred.Test,family = 'binomial'))$R2, 4))
      }
    }
  }
  colnames(pred.Test.cp) <- obj.SNP$fam[,1][ind.test]
  colnames(pred.Train.cp) <- obj.SNP$fam[,1][ind.train]

  id = which.max(rfnr2)
  rf.nr2 <- rfnr2[id]

  return(pred.Test[id,])
}
