#' Clumping and Thresholding Random Forest (ctRF)
#'
#' @param G Genotype matrix.
#' @param CHR Chromosome data.
#' @param POS Position data.
#' @param train_indices Training indices.
#' @param test_indices Test indices.
#' @param p_values P-values for SNPs.
#' @param covariates Covariates like age and sex.
#' @param grid_type Type of grid for hyperparameter tuning: 'light', 'medium', or 'intensive'.
#' @param n_cores Number of cores for parallel processing.
#'
#' @return A list containing the best model and predictions.
#' @export

ctrf <- function(G, CHR, POS, train_indices, test_indices, p_values, covariates, grid_type = "medium", n_cores = 4) {
  # Preprocessing: Clumping & Thresholding logic based on LD and p-values
  
  clumped_indices <- snp_clumping(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.01)  # Adjust according to grid_type
  
  # Grid setup: based on grid_type (intensive, medium, light)
  hyper_grid <- switch(grid_type,
    "intensive" = expand.grid(mtry = c(sqrt(p)/2, sqrt(p), sqrt(p)*2), num.trees = c(5000, 10000)),
    "medium" = expand.grid(mtry = c(sqrt(p), sqrt(p)*2), num.trees = c(1000, 5000)),
    "light" = expand.grid(mtry = c(sqrt(p)), num.trees = c(1000))
  )
  
  # Random Forest with cross-validation and hyperparameter tuning
  best_model <- NULL
  best_rmse <- Inf
  
  for (i in 1:nrow(hyper_grid)) {
    model <- ranger(x = G[train_indices, clumped_indices], y = as.factor(covariates$pheno[train_indices]),
                    mtry = hyper_grid$mtry[i], num.trees = hyper_grid$num.trees[i], importance = "impurity_corrected",
                    classification = TRUE, probability = TRUE)
    
    if (model$prediction.error < best_rmse) {
      best_rmse <- model$prediction.error
      best_model <- model
    }
  }
  
  # Return the best model and predictions
  predictions <- predict(best_model, G[test_indices, clumped_indices])$predictions[,2]
  return(list(model = best_model, predictions = predictions))
}
