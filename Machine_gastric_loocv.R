library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(linkET)



##################################
#### preparation ####
##################################


data<-readRDS('subset_input_ML_dataframes_after_univariate.RDS') # obtain from "ML_input_preparation.R"
data <- lapply(data,function(x){
  x <- x[x$time != 0, ]
  
  x[,-c(1:2)] <- scale(x[,-c(1:2)])
  return(x)})
data$tcga<-na.omit(data$tcga)
##################################
#### preparation ####
##################################

# result <- data.frame()
 est_dd <- data$tcga
# 
# val_data_list <- data
 selected_gene <- colnames(est_dd)[-c(1:2)]
# val_dd_list <- lapply(val_data_list,function(x){x[,c('time','status',selected_gene),drop = FALSE]})
# 



rf_nodesize <- 5
seed <- 1


convert_to_numeric <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col) || is.character(col)) {
      as.numeric(as.character(col))  # Convert factor/character to numeric
    } else {
      as.numeric(col)  # Ensure all other columns are numeric
    }
  })
  return(df)
}
# # Apply the conversion function to your training and validation data
# est_dd <- convert_to_numeric(data$tcga)


val_dd_list <- lapply(data, function(df) {
  convert_to_numeric(df)  # Convert all columns to numeric
})

##################################
#### 1.RSF ####
##################################
result <- data.frame()

# Perform dataset-level LOOCV
rs_list <- lapply(seq_along(data), function(i) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Train RSF model
  fit <- rfsrc(Surv(time, status) ~ ., data = est_dd,
               ntree = 1000, nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = TRUE, proximity = TRUE, forest = TRUE, seed = seed)
  
  # Predict risk scores for the validation dataset
  cbind(val_dd[, 1:2], RS = predict(fit, newdata = val_dd)$predicted)
})

# Calculate C-index for each dataset
cc <- data.frame(Cindex = sapply(rs_list, function(x) {
  as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])
})) %>%
  rownames_to_column('ID')

# Add model name
cc$ID <- names(data)  # Match dataset names
cc$Model <- 'RSF'

# Append results
result <- rbind(result, cc)

##################################
#### 2.LASSO ####
##################################


# Perform LOOCV
rs_list <- lapply(seq_along(data), function(i) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Prepare data for glmnet
  est_dd$time[est_dd$time <= 0] <- 0.1  # Replace non-positive times with a small positive value
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Features (excluding time and status)
  y_est <- Surv(est_dd$time, est_dd$status)  # Surv object for training
  
  # Fit LASSO model with cross-validation
  set.seed(seed)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  fit <- glmnet(x_est, y_est, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
  
  # Predict risk scores for the validation dataset
  x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation features
  preds <- predict(fit, newx = x_val, s = "lambda.min", type = "link")  # Predicted risk score
  cbind(val_dd[, 1:2], RS = as.vector(preds))  # Combine with time and status
})

# Calculate C-index for each dataset
cc <- data.frame(Cindex = sapply(rs_list, function(x) {
  as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])
})) %>%
  rownames_to_column('ID')

# Add model identifier
cc$ID <- names(data)  # Match dataset names
cc$Model <- 'LASSO'

# Combine with the result object
result <- rbind(result, cc)
##################################
#### 3.RIDGE ####
##################################

# Fit Ridge (alpha = 0 for Ridge, use cv.glmnet for cross-validation)
set.seed(seed)
rs_list <- lapply(seq_along(data), function(i) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Prepare data for glmnet
  est_dd$time[est_dd$time <= 0] <- 0.1  # Replace non-positive times with a small positive value
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Features (excluding time and status)
  y_est <- Surv(est_dd$time, est_dd$status)  # Surv object for training
  
  # Fit Ridge regression model with cross-validation
  set.seed(seed)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 0, nfolds = 10)
  fit <- glmnet(x_est, y_est, family = "cox", alpha = 0, lambda = cvfit$lambda.min)
  
  # Predict risk scores for the validation dataset
  x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation features
  preds <- predict(fit, newx = x_val, s = "lambda.min", type = "link")  # Predicted risk score
  cbind(val_dd[, 1:2], RS = as.vector(preds))  # Combine with time and status
})

# Calculate C-index for each dataset
cc <- data.frame(Cindex = sapply(rs_list, function(x) {
  as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])
})) %>%
  rownames_to_column('ID')

# Add model identifier
cc$ID <- names(data)  # Match dataset names
cc$Model <- 'Ridge'

# Combine with the result object
result <- rbind(result, cc)

##################################
#### 4.Enet ####
##################################

for (alpha in seq(0, 1, 0.1)) {  # Iterate over different alpha values
  rs_list <- lapply(seq_along(data), function(i) {
    # Leave one dataset out for validation
    val_dd <- data[[i]]
    train_dd_list <- data[-i]  # All other datasets for training
    
    # Combine remaining datasets into a single training set
    est_dd <- do.call(rbind, train_dd_list)
    
    # Convert all columns to numeric
    est_dd <- convert_to_numeric(est_dd)
    val_dd <- convert_to_numeric(val_dd)
    
    # Prepare training data
    est_dd$time[est_dd$time <= 0] <- 0.1  # Replace non-positive times with a small positive value
    x_train <- as.matrix(est_dd[, selected_gene])  # Features
    y_train <- Surv(est_dd$time, est_dd$status)  # Surv object for training
    
    # Fit Elastic Net model with cross-validation
    set.seed(seed)
    fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = alpha, nfolds = 10)
    
    # Predict risk scores for the validation dataset
    x_val <- as.matrix(val_dd[, selected_gene])  # Validation features
    preds <- predict(fit, newx = x_val, s = fit$lambda.min, type = "link")  # Predicted risk scores
    cbind(val_dd[, 1:2], RS = as.vector(preds))  # Combine with time and status
  })
  
  # Calculate C-index for each dataset
  cc <- data.frame(Cindex = sapply(rs_list, function(x) {
    as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])
  })) %>%
    rownames_to_column('ID')
  
  # Add model identifier
  cc$ID <- names(data)  # Match dataset names
  cc$Model <- paste0('Enet[α=', alpha, ']')  # Identify model with specific alpha value
  
  # Combine with the result object
  result <- rbind(result, cc)
}

##################################
#### 5.StepCox ####
##################################

for (direction in c("both", "backward", "forward")) {  # Iterate over different stepwise directions
  rs_list <- lapply(seq_along(data), function(i) {
    # Leave one dataset out for validation
    val_dd <- data[[i]]
    train_dd_list <- data[-i]  # All other datasets for training
    
    # Combine remaining datasets into a single training set
    est_dd <- do.call(rbind, train_dd_list)
    
    # Convert all columns to numeric
    est_dd <- convert_to_numeric(est_dd)
    val_dd <- convert_to_numeric(val_dd)
    
    # Prepare the Cox model using stepwise selection
    fit <- step(coxph(Surv(time, status) ~ ., data = est_dd), direction = direction)
    
    # Predict risk scores for the validation dataset
    rs <- predict(fit, newdata = val_dd, type = 'risk')  # Predicted risk scores
    cbind(val_dd[, 1:2], RS = rs)  # Combine with time and status
  })
  
  # Calculate C-index for each dataset
  cc <- data.frame(Cindex = sapply(rs_list, function(x) {
    as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])
  })) %>%
    rownames_to_column('ID')
  
  # Add model identifier
  cc$ID <- names(data)  # Match dataset names
  cc$Model <- paste0('StepCox[', direction, ']')  # Identify model with specific direction
  
  # Combine with the result object
  result <- rbind(result, cc)
}



##################################
#### 6.CoxBoost ####
##################################

set.seed(seed)
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # CoxBoost penalty selection
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                              trace = TRUE, start.penalty = 500, parallel = TRUE)
  
  # Cross-validation to find the optimal penalty
  cv.res <- cv.CoxBoost(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                        maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
  
  # Fit the CoxBoost model using the optimal penalty and step number
  fit <- CoxBoost(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  
  # Apply the fitted CoxBoost model to validation datasets
  rs <- cbind(val_dd[,1:2], RS = as.numeric(predict(fit, newdata = val_dd[,-c(1, 2)],
                                                    newtime = val_dd[, 1], newstatus = val_dd[, 2], type = "lp")))
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  
  # Add model identifier
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'CoxBoost'
  
  # Combine with the result object
  result <- rbind(result, cc)
}

##################################
#### 7.plsRcox####
##################################

set.seed(seed)
# Perform LOOCV for plsRcox
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Fit plsRcox model using cross-validation to select optimal number of components (nt)
  set.seed(seed)
  cv.plsRcox.res <- cv.plsRcox(list(x = est_dd[, selected_gene], time = est_dd$time, status = est_dd$status),
                               nt = 10, verbose = FALSE)
  
  # Fit the plsRcox model with optimal number of components (nt)
  fit <- plsRcox(est_dd[, selected_gene], time = est_dd$time, event = est_dd$status, 
                 nt = as.numeric(cv.plsRcox.res[5]))
  
  # Apply the fitted plsRcox model to validation datasets
  rs <- cbind(val_dd[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = val_dd[, -c(1, 2)])))
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  
  # Add model identifier
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'plsRcox'
  
  # Combine with the result object
  result <- rbind(result, cc)
}

##################################
#### 8.superpc####
##################################

for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Prepare data for SuperPC
  data1 <- list(x = t(est_dd[, -c(1, 2)]), y = est_dd$time, censoring.status = est_dd$status, 
                featurenames = colnames(est_dd)[-c(1, 2)])
  
  # Train SuperPC model
  set.seed(seed)
  fit <- superpc.train(data = data1, type = 'survival', s0.perc = 0.5)
  
  # Perform cross-validation for SuperPC to determine best threshold
  cv.fit <- superpc.cv(fit, data1, n.threshold = 20, n.fold = 10, n.components = 3, 
                       min.features = 5, max.features = nrow(data1$x), compute.fullcv = TRUE, 
                       compute.preval = TRUE)
  
  # Apply SuperPC model to validation datasets
  test <- list(x = t(val_dd[, -c(1, 2)]), y = val_dd$time, censoring.status = val_dd$status, 
               featurenames = colnames(val_dd)[-c(1, 2)])
  ff <- superpc.predict(fit, data1, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(val_dd[, 1:2], RS = rr)
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rr2))$concordance[1])) %>%
    rownames_to_column('ID')
  
  # Add model identifier
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'SuperPC'
  
  # Combine with the result object
  result <- rbind(result, cc)
}


##################################
#### 9.GBM ####
##################################

set.seed(seed)

for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Train GBM model
  fit <- gbm(formula = Surv(time, status) ~ ., data = est_dd, distribution = 'coxph',
             n.trees = 10000, interaction.depth = 3, n.minobsinnode = 10, 
             shrinkage = 0.001, cv.folds = 10, n.cores = 6)
  
  # Find the optimal number of trees based on minimum CV error
  best <- which.min(fit$cv.error)
  
  # Fit GBM model again with the optimal number of trees
  set.seed(seed)
  fit <- gbm(formula = Surv(time, status) ~ ., data = est_dd, distribution = 'coxph',
             n.trees = best, interaction.depth = 3, n.minobsinnode = 10,
             shrinkage = 0.001, cv.folds = 10, n.cores = 8)
  
  # Apply the fitted model to the single held-out validation dataset
  preds <- predict(fit, val_dd, n.trees = best, type = 'link')
  val_dd_with_rs <- cbind(val_dd[, 1:2], RS = as.numeric(preds))  # Add risk scores
  
  # Calculate C-index for the single validation dataset
  cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd_with_rs))$concordance[1])
  
  # Create a result row for this iteration
  cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = 'GBM')
  
  # Combine with the result object
  result <- rbind(result, cc)
}


##################################
#### 10.survivalsvm ####
##################################

# Perform LOOCV for survival-SVM
set.seed(seed)

for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  
  # Train the Survival-SVM model
  fit <- survivalsvm(Surv(time, status) ~ ., data = est_dd, gamma.mu = 1)
  
  # Apply the fitted model to the single held-out validation dataset
  preds <- as.numeric(predict(fit, val_dd)$predicted)
  val_dd_with_rs <- cbind(val_dd[, 1:2], RS = preds)  # Add risk scores
  
  # Calculate C-index for the single validation dataset
  cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd_with_rs))$concordance[1])
  
  # Create a result row for this iteration
  cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = 'survival-SVM')
  
  # Combine with the result object
  result <- rbind(result, cc)
}

##################################
#### 11.LASSO + RSF ####
##################################

# Fit LASSO (Elastic Net) to select important features
set.seed(seed)
# Perform LOOCV for LASSO + RSF model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  
  # Fit LASSO model
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Remove OS.time and OS columns
  y_est <- Surv(est_dd$time, est_dd$status)  # Surv object for training
  
  # Apply LASSO
  set.seed(seed)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  
  # Get the selected features
  rid <- rownames(coef(cvfit, s = "lambda.min"))[which(coef(cvfit, s = "lambda.min") != 0)]  # Selected features
  
  # Subset the data with the selected features
  est_dd2 <- est_dd[, c('time', 'status', rid)]
  val_dd2 <- val_dd[, c('time', 'status', rid)]
  
  # Fit RSF model using selected features from LASSO
  set.seed(seed)
  fit <- rfsrc(Surv(time, status) ~ ., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, 
               splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)
  
  # Predict using the RSF model on the validation dataset
  rs <- cbind(val_dd2[, 1:2], RS = predict(fit, newdata = val_dd2)$predicted)
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'LASSO + RSF'
  
  # Combine with the result object
  result <- rbind(result, cc)
}



################
### 12. Ridge + RSF ####
######################

# Set seed for reproducibility

set.seed(seed)
# Perform LOOCV for Ridge + RSF model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  # Fit Ridge regression (alpha = 0)
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Remove OS.time and OS columns
  y_est <- Surv(est_dd$time, est_dd$status)  # Create Surv object for training
  
  # Apply Ridge regression
  set.seed(seed)
  cvfit_ridge <- cv.glmnet(x_est, y_est, family = "cox", alpha = 0, nfolds = 10)
  
  # Get the selected features
  rid <- rownames(coef(cvfit_ridge, s = "lambda.min"))[which(coef(cvfit_ridge, s = "lambda.min") != 0)]  # Selected features
  
  # Subset the data with the selected features
  est_dd2 <- est_dd[, c('time', 'status', rid)]
  val_dd2 <- val_dd[, c('time', 'status', rid)]
  
  # Fit RSF model using selected features from Ridge
  set.seed(seed)
  fit_rsf <- rfsrc(Surv(time, status) ~ ., data = est_dd2,
                   ntree = 1000, nodesize = rf_nodesize,  ## adjust this value if needed
                   splitrule = 'logrank', importance = TRUE, 
                   proximity = TRUE, forest = TRUE)
  
  # Predict using the RSF model on the validation dataset
  rs <- cbind(val_dd2[, 1:2], RS = predict(fit_rsf, newdata = val_dd2)$predicted)
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'Ridge + RSF'
  
  # Combine with the result object
  result <- rbind(result, cc)
}


##################################
#### 13. RSF + Enet ####
##################################

set.seed(seed)
# Perform LOOCV for RSF + Enet (with varying alpha) model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  for (alpha in seq(0, 1, 0.1)) {
    # Fit Elastic Net (Enet) model with selected alpha
    x1 <- as.matrix(est_dd[, selected_gene])  # Use selected genes for features
    x2 <- as.matrix(Surv(est_dd$time, est_dd$status))  # Create Surv object for training
    
    set.seed(seed)
    fit <- cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
    
    # Get the selected features from Enet
    rid <- rownames(coef(fit, s = "lambda.min"))[which(coef(fit, s = "lambda.min") != 0)]
    est_dd2 <- est_dd[, c('time', 'status', rid)]
    val_dd2 <- val_dd[, c('time', 'status', rid)]
    
    # Fit RSF model using selected features from Enet
    set.seed(seed)
    fit_rsf <- rfsrc(Surv(time, status) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize,
                     splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)
    
    # Predict using the RSF model on the validation dataset
    rs <- cbind(val_dd2[, 1:2], RS = predict(fit_rsf, newdata = val_dd2)$predicted)
    
    # Calculate C-index for the validation dataset
    cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
    cc$ID <- names(data)[i]  # Match dataset names
    cc$Model <- paste0('RSF + Enet', '[α=', alpha, ']')
    
    # Combine with the result object
    result <- rbind(result, cc)
  }
}

##################################
#### 14. RSF + StepCox ####
##################################

set.seed(seed)
# Perform LOOCV for RSF + StepCox (with different directions) model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  for (direction in c("both", "backward", "forward")) {
    # Fit StepCox model for feature selection
    fit <- step(coxph(Surv(time, status) ~ ., est_dd), direction = direction)
    
    # Get the selected features from StepCox
    rid <- names(coef(fit))[-1]  # Exclude intercept
    est_dd2 <- est_dd[, c('time', 'status', rid)]
    val_dd2 <- val_dd[, c('time', 'status', rid)]
    
    # Fit RSF model using selected features from StepCox
    set.seed(seed)
    fit_rsf <- rfsrc(Surv(time, status) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize,
                     splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)
    
    # Predict using the RSF model on the validation dataset
    rs <- cbind(val_dd2[, 1:2], RS = predict(fit_rsf, newdata = val_dd2)$predicted)
    
    # Calculate C-index for the validation dataset
    cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
    cc$ID <- names(data)[i]  # Match dataset names
    cc$Model <- paste0('RSF + StepCox', '[', direction, ']')
    
    # Combine with the result object
    result <- rbind(result, cc)
  }
}

##################################
#### 15. RSF + CoxBoost ####
##################################

set.seed(seed)
# Perform LOOCV for RSF + CoxBoost model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  # Apply CoxBoost for feature selection
  pen <- optimCoxBoostPenalty(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                              trace = TRUE, start.penalty = 500, parallel = TRUE)
  cv.res <- cv.CoxBoost(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                        maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
  fit_cboost <- CoxBoost(est_dd[,'time'], est_dd[,'status'], as.matrix(est_dd[,-c(1, 2)]),
                         stepno = cv.res$optimal.step, penalty = pen$penalty)
  
  # Get the selected features from CoxBoost
  selected_features <- colnames(est_dd)[-c(1, 2)][cv.res$selected.variables]
  if (length(selected_features) > 0) {
  
  # Subset the data with selected features
  est_dd2 <- est_dd[, c('time', 'status', selected_features)]
  val_dd2 <- val_dd[, c('time', 'status', selected_features)]
  
  # Fit RSF model using selected features from CoxBoost
  set.seed(seed)
  fit_rsf <- rfsrc(Surv(time, status) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize,
                   splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)
  
  # Predict using the RSF model on the validation dataset
  rs <- cbind(val_dd2[, 1:2], RS = predict(fit_rsf, newdata = val_dd2)$predicted)
  
  # Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
  cc$ID <- names(data)[i]  # Match dataset names
  cc$Model <- 'RSF + CoxBoost'
  
  # Combine with the result object
  result <- rbind(result, cc)
  }
  else {
    message("No selected features from plsRcox, skipping this fold.")
  }
}


##################################
#### 16. RSF + plsRcox ####
##################################

set.seed(seed)
# Perform LOOCV for RSF + plsRcox model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1
  
  # Apply plsRcox for feature selection
  cv.plsRcox.res <- cv.plsRcox(list(x = est_dd[, selected_gene], time = est_dd$time, status = est_dd$status),
                               nt = 10, verbose = FALSE)
  fit_plsrc <- plsRcox(est_dd[, selected_gene], time = est_dd$time, event = est_dd$status, 
                       nt = as.numeric(cv.plsRcox.res[5]))
  
  # Get the selected features from plsRcox
  selected_features <- names(est_dd)[-c(1, 2)][cv.plsRcox.res$selected]
  
  # Ensure that selected features are not empty
  if (length(selected_features) > 0) {
    # Subset the data with selected features
    est_dd2 <- est_dd[, c('time', 'status', selected_features)]
    val_dd2 <- val_dd[, c('time', 'status', selected_features)]
    
    # Fit RSF model using selected features from plsRcox
    set.seed(seed)
    fit_rsf <- rfsrc(Surv(time, status) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize,
                     splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)
    
    # Predict using the RSF model on the validation dataset
    rs <- cbind(val_dd2[, 1:2], RS = predict(fit_rsf, newdata = val_dd2)$predicted)
    
    # Calculate C-index for the validation dataset
    cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1]))
    cc$ID <- names(data)[i]  # Match dataset names
    cc$Model <- 'RSF + plsRcox'
    
    # Combine with the result object
    result <- rbind(result, cc)
  } else {
    message("No selected features from plsRcox, skipping this fold.")
  }
}



##################################
#### 17. RSF + SuperPC ####
##################################
set.seed(seed)

# Perform LOOCV for RSF + SuperPC model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure no zero or negative survival times
  
  # Prepare the data for SuperPC
  data2 <- list(
    x = t(est_dd[, -c(1, 2)]), 
    y = est_dd$time, 
    censoring.status = est_dd$status, 
    featurenames = colnames(est_dd)[-c(1, 2)]
  )
  
  # Train the SuperPC model
  fit_superpc <- superpc.train(data = data2, type = 'survival', s0.perc = 0.5)
  
  # Perform cross-validation for SuperPC
  cv.fit <- superpc.cv(
    fit_superpc, data2, 
    n.threshold = 20, n.fold = 10, 
    n.components = 3, 
    min.features = 5, 
    max.features = nrow(data2$x), 
    compute.fullcv = TRUE, 
    compute.preval = TRUE
  )
  
  # Select the optimal threshold
  optimal_threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])]
  
  # Predict using the SuperPC model on the single validation dataset
  test <- list(
    x = t(val_dd[, -c(1, 2)]), 
    y = val_dd$time, 
    censoring.status = val_dd$status, 
    featurenames = colnames(val_dd)[-c(1, 2)]
  )
  pred <- superpc.predict(fit_superpc, data2, test, threshold = optimal_threshold, n.components = 1)
  val_dd_with_rs <- cbind(val_dd[, 1:2], RS = as.numeric(pred$v.pred))  # Add risk scores
  
  # Calculate C-index for the validation dataset
  cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd_with_rs))$concordance[1])
  
  # Create a result row for this iteration
  cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = 'RSF + SuperPC')
  
  # Combine with the result object
  result <- rbind(result, cc)
}



##################################
#### 18. RSF + GBM ####
##################################

set.seed(seed)

# Perform LOOCV for RSF + GBM model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure no zero or negative survival times
  
  # Train GBM model on the training dataset
  set.seed(seed)
  fit_gbm <- gbm(
    formula = Surv(time, status) ~ ., 
    data = est_dd, 
    distribution = "coxph", 
    n.trees = 10000, 
    interaction.depth = 3, 
    n.minobsinnode = 10, 
    shrinkage = 0.001, 
    cv.folds = 10, 
    n.cores = 6
  )
  
  # Find the number of trees with minimum CV error
  best <- which.min(fit_gbm$cv.error)
  
  # Retrain GBM with the optimal number of trees
  set.seed(seed)
  fit_gbm <- gbm(
    formula = Surv(time, status) ~ ., 
    data = est_dd, 
    distribution = "coxph", 
    n.trees = best, 
    interaction.depth = 3, 
    n.minobsinnode = 10, 
    shrinkage = 0.001, 
    cv.folds = 10, 
    n.cores = 8
  )
  
  # Predict using the GBM model on the single validation dataset
  val_dd_with_rs <- cbind(
    val_dd[, 1:2], 
    RS = as.numeric(predict(fit_gbm, val_dd, n.trees = best, type = "link"))
  )
  
  # Calculate C-index for the validation dataset
  cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd_with_rs))$concordance[1])
  
  # Create a result row for this iteration
  cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = 'RSF + GBM')
  
  # Combine with the result object
  result <- rbind(result, cc)
}




##################################
#### 19. RSF + survival-SVM ####
##################################

set.seed(seed)

# Perform LOOCV for RSF + Survival-SVM model
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine remaining datasets into a single training set
  est_dd <- do.call(rbind, train_dd_list)
  
  # Convert all columns to numeric for both training and validation
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure no zero or negative survival times
  
  # Fit RSF model on the training dataset
  set.seed(seed)
  fit_rsf <- rfsrc(
    Surv(time, status) ~ ., 
    data = est_dd, 
    ntree = 1000, 
    nodesize = rf_nodesize, 
    splitrule = "logrank", 
    importance = TRUE, 
    proximity = TRUE, 
    forest = TRUE
  )
  
  # Predict risk scores using the RSF model on the single validation dataset
  val_dd_with_rsf <- cbind(
    val_dd[, 1:2], 
    RS = predict(fit_rsf, newdata = val_dd)$predicted
  )
  
  # Fit Survival-SVM using the RSF risk scores
  fit_svm <- survivalsvm(Surv(time, status) ~ RS, data = val_dd_with_rsf, gamma.mu = 1)
  
  # Calculate C-index for the validation dataset
  cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd_with_rsf))$concordance[1])
  
  # Create a result row for this iteration
  cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = 'RSF + Survival-SVM')
  
  # Combine with the result object
  result <- rbind(result, cc)
}


#################33
## 20. lasso + survival-svm
################


###############
## 21. lasso + stepcox
################

set.seed(seed)


# Perform LOOCV for LASSO + StepCox
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Prepare matrix for LASSO
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Create survival object
  
  # Fit LASSO and get selected features
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  fit_lasso <- glmnet(x_est, y_est, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
  
  # Extract non-zero coefficients
  rid_lasso <- rownames(as.matrix(coef(fit_lasso, s = cvfit$lambda.min)))[as.matrix(coef(fit_lasso, s = cvfit$lambda.min)) != 0]
  rid_lasso <- rid_lasso[rid_lasso != "(Intercept)"]  # Remove intercept if present
  
  if (length(rid_lasso) > 0) {
    # Subset training and validation data with LASSO-selected features
    est_dd2 <- est_dd[, c("time", "status", rid_lasso), drop = FALSE]
    val_dd2 <- val_dd[, c("time", "status", rid_lasso), drop = FALSE]
    
    # Iterate over stepwise directions
    for (direction in c("both", "backward", "forward")) {
      # Fit StepCox on LASSO-selected features
      fit_step <- step(coxph(Surv(time, status) ~ ., est_dd2), direction = direction)
      
      # Fit Survival-SVM using the features from StepCox
      fit_svm <- survivalsvm(Surv(time, status) ~ ., data = est_dd2, gamma.mu = 1)
      
      # Predict risk scores for validation data
      val_dd2$RS <- as.numeric(predict(fit_svm, newdata = val_dd2)$predicted)
      
      # Calculate C-index for validation data
      cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, val_dd2))$concordance[1])
      
      # Store results
      cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = paste0("LASSO + StepCox [", direction, "]"))
      result <- rbind(result, cc)
    }
  } else {
    # Log if LASSO selects no features
    message("No LASSO features selected for fold: ", names(data)[i])
  }
}



#################
## 22. lasso + coxboost
################

# Fit LASSO model
set.seed(seed)


# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Prepare matrix for LASSO
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Create survival object
  
  # Fit LASSO and get selected features
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  fit_lasso <- glmnet(x_est, y_est, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
  
  # Extract non-zero coefficients
  rid_lasso <- rownames(as.matrix(coef(fit_lasso, s = cvfit$lambda.min)))[as.matrix(coef(fit_lasso, s = cvfit$lambda.min)) != 0]
  rid_lasso <- rid_lasso[rid_lasso != "(Intercept)"]  # Remove intercept if present
  
  if (length(rid_lasso) > 0) {
    # Subset training and validation data with LASSO-selected features
    est_dd2 <- est_dd[, c("time", "status", rid_lasso), drop = FALSE]
    val_dd2 <- val_dd[, c("time", "status", rid_lasso), drop = FALSE]
    
    # Fit CoxBoost on LASSO-selected features
    pen <- optimCoxBoostPenalty(est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]), trace = TRUE, start.penalty = 500, parallel = TRUE)
    cv.res <- cv.CoxBoost(est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
    fit_coxboost <- CoxBoost(est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, penalty = pen$penalty)
    
    # Predict risk scores using CoxBoost
    rs <- cbind(val_dd2[, 1:2], RS = as.numeric(predict(fit_coxboost, newdata = as.matrix(val_dd2[, -c(1, 2)]), newtime = val_dd2$time, newstatus = val_dd2$status, type = "lp")))
    
    # Calculate C-index for validation dataset
    cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])
    
    # Store results
    cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = "LASSO + CoxBoost")
    result <- rbind(result, cc)
  } else {
    # Log if LASSO selects no features
    message("No LASSO features selected for fold: ", names(data)[i])
  }
}

#################
## 23. lasso + plsRcox
################

set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit LASSO model on training data
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  rid <- rownames(as.matrix(coef(cvfit, s = "lambda.min")))[as.matrix(coef(cvfit, s = "lambda.min")) != 0]
  rid <- rid[rid != "(Intercept)"]  # Remove intercept if present
  
  if (length(rid) > 0) {
    # Step 2: Subset training and validation data with LASSO-selected features
    if (all(rid %in% colnames(est_dd))) {
      est_dd2 <- est_dd[, c("time", "status", rid), drop = FALSE]
      val_dd2 <- val_dd[, c("time", "status", rid), drop = FALSE]
      
      # Step 3: Fit plsRcox model using selected features from LASSO
      set.seed(seed)
      cv.plsRcox.res <- cv.plsRcox(
        list(x = est_dd2[, -c(1, 2)], time = est_dd2$time, status = est_dd2$status), 
        nt = 10, verbose = FALSE
      )
      fit <- plsRcox(est_dd2[, -c(1, 2)], time = est_dd2$time, event = est_dd2$status, nt = as.numeric(cv.plsRcox.res[5]))
      
      # Step 4: Predict using the plsRcox model on validation dataset
      rs <- cbind(val_dd2[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = val_dd2[, -c(1, 2)])))
      
      # Step 5: Calculate C-index for the validation dataset
      cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])
      
      # Store results
      cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = "LASSO + plsRcox")
      result <- rbind(result, cc)
    } else {
      message("Selected features not found in training dataset for fold: ", names(data)[i])
    }
  } else {
    # Log if LASSO selects no features
    message("No LASSO features selected for fold: ", names(data)[i])
  }
}



######################
### 24. lasso + superPC 
#######################

# Step 1: Fit LASSO and select features
set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit LASSO model on training data
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  rid <- rownames(as.matrix(coef(cvfit, s = "lambda.min")))[as.matrix(coef(cvfit, s = "lambda.min")) != 0]
  rid <- rid[rid != "(Intercept)"]  # Remove intercept if present
  
  if (length(rid) > 0) {
    # Step 2: Subset training and validation data with LASSO-selected features
    est_dd2 <- est_dd[, c("time", "status", rid), drop = FALSE]
    val_dd2 <- val_dd[, c("time", "status", rid), drop = FALSE]
    
    # Step 3: Prepare data for superPC
    data3 <- list(
      x = t(est_dd2[, -c(1, 2)]), 
      y = est_dd2$time, 
      censoring.status = est_dd2$status, 
      featurenames = colnames(est_dd2)[-c(1, 2)]
    )
    
    # Fit superPC model
    fit <- superpc.train(data = data3, type = 'survival', s0.perc = 0.5)
    
    # Cross-validation for superPC
    cv.fit <- superpc.cv(
      fit, data3, n.threshold = 20, 
      n.fold = 10, 
      n.components = 3, 
      min.features = 5, 
      max.features = nrow(data3$x), 
      compute.fullcv = TRUE, 
      compute.preval = TRUE
    )
    
    # Step 4: Predict risk scores using superPC on the validation dataset
    test <- list(
      x = t(val_dd2[, -c(1, 2)]), 
      y = val_dd2$time, 
      censoring.status = val_dd2$status, 
      featurenames = colnames(val_dd2)[-c(1, 2)]
    )
    
    # Predict risk scores
    ff <- superpc.predict(
      fit, data3, test, 
      threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], 
      n.components = 1
    )
    rr <- as.numeric(ff$v.pred)
    
    # Combine with time and status
    rs <- cbind(val_dd2[, 1:2], RS = rr)
    
    # Step 5: Calculate C-index for validation dataset
    cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])
    
    # Store results
    cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = "LASSO + SuperPC")
    result <- rbind(result, cc)
  } else {
    # Log if LASSO selects no features
    message("No LASSO features selected for fold: ", names(data)[i])
  }
}



#####################
##### 25. LASSO + GBM ###
#####################

# Step 1: Fit LASSO and select features
set.seed(seed)


# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit LASSO model on training data
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  rid <- rownames(as.matrix(coef(cvfit, s = "lambda.min")))[as.matrix(coef(cvfit, s = "lambda.min")) != 0]
  rid <- rid[rid != "(Intercept)"]  # Remove intercept if present
  
  if (length(rid) > 0) {
    # Step 2: Subset training and validation data with LASSO-selected features
    if (all(rid %in% colnames(est_dd))) {
      est_dd2 <- est_dd[, c("time", "status", rid), drop = FALSE]
      val_dd2 <- val_dd[, c("time", "status", rid), drop = FALSE]
      
      # Step 3: Fit GBM model using selected features from LASSO
      set.seed(seed)
      fit <- gbm(
        formula = Surv(time, status) ~ ., data = est_dd2, distribution = 'coxph',
        n.trees = 10000,
        interaction.depth = 3,
        n.minobsinnode = 10,
        shrinkage = 0.001,
        cv.folds = 10, n.cores = 6
      )
      
      # Find index for number of trees with minimum CV error
      best <- which.min(fit$cv.error)
      
      # Refit GBM with the optimal number of trees
      set.seed(seed)
      fit <- gbm(
        formula = Surv(time, status) ~ ., data = est_dd2, distribution = 'coxph',
        n.trees = best,
        interaction.depth = 3,
        n.minobsinnode = 10,
        shrinkage = 0.001,
        cv.folds = 10, n.cores = 8
      )
      
      # Step 4: Predict using the GBM model on validation dataset
      rs <- cbind(val_dd2[, 1:2], RS = as.numeric(predict(fit, val_dd2, n.trees = best, type = 'link')))
      
      # Step 5: Calculate C-index for the validation dataset
      cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])
      
      # Store results
      cc <- data.frame(ID = names(data)[i], Cindex = cindex, Model = "LASSO + GBM")
      result <- rbind(result, cc)
    } else {
      message("Selected features not found in training dataset for fold: ", names(data)[i])
    }
  } else {
    # Log if LASSO selects no features
    message("No LASSO features selected for fold: ", names(data)[i])
  }
}



############################
### 26. Ridge + stepcox ####
############################
# Ridge regression is performed by setting alpha = 0 in glmnet
# Ridge regression
set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit Ridge model on training data
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  # Perform Ridge regression (alpha = 0)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 0, nfolds = 10)
  
  # Use the best lambda from cross-validation to fit Ridge model
  fit_ridge <- glmnet(x_est, y_est, family = "cox", alpha = 0, lambda = cvfit$lambda.min)
  
  # Step 2: Predict risk scores using Ridge model on validation data
  x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation data excluding time and status
  preds_ridge <- predict(fit_ridge, newx = x_val, s = "lambda.min", type = "link")
  rs_ridge <- cbind(val_dd[, 1:2], RS = as.vector(preds_ridge))  # Combine with OS.time and OS
  
  # Step 3: Apply Stepwise Cox to the original dataset (training data)
  for (direction in c("both", "backward", "forward")) {
    # Fit Stepwise Cox model on training data
    fit_stepcox <- step(coxph(Surv(time, status) ~ ., est_dd), direction = direction)
    
    # Predict risk scores using the Stepwise Cox model on validation data
    rs_stepcox <- cbind(val_dd[, 1:2], RS = predict(fit_stepcox, type = 'risk', newdata = val_dd))
    
    # Calculate C-index for the validation dataset
    cindex_stepcox <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs_stepcox))$concordance[1])
    
    # Create the result dataframe with consistent columns
    cc <- data.frame(ID = names(data)[i], Cindex = cindex_stepcox, Model = paste0('Ridge + StepCox[', direction, ']'))
    
    # Combine with the result object
    result <- rbind(result, cc)
  }
}






############################
### 27. Ridge + coxboost ####
############################

# Ridge regression
set.seed(seed)


# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit Ridge model on training data
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  y_est <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  # Perform Ridge regression (alpha = 0)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 0, nfolds = 10)
  
  # Use the best lambda from cross-validation to fit Ridge model
  fit_ridge <- glmnet(x_est, y_est, family = "cox", alpha = 0, lambda = cvfit$lambda.min)
  
  # Step 2: Predict risk scores using Ridge model on validation data
  x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation data excluding time and status
  preds_ridge <- predict(fit_ridge, newx = x_val, s = "lambda.min", type = "link")
  rs_ridge <- cbind(val_dd[, 1:2], RS = as.vector(preds_ridge))  # Combine with OS.time and OS
  
  # Step 3: Fit CoxBoost model on training data
  set.seed(seed)
  fit_coxboost <- CoxBoost(
    time = est_dd$time,
    status = est_dd$status,
    x = as.matrix(est_dd[, -c(1, 2)]),
    stepno = 100
  )
  
  # Step 4: Predict risk scores using CoxBoost model on validation data
  preds_coxboost <- predict(fit_coxboost, newdata = x_val, type = "lp")
  rs_coxboost <- cbind(val_dd[, 1:2], RS = as.vector(preds_coxboost))  # Combine with OS.time and OS
  
  # Step 5: Calculate C-index for the validation dataset using CoxBoost predictions
  cindex_coxboost <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs_coxboost))$concordance[1])
  
  # Step 6: Add results to the result data frame
  cc <- data.frame(ID = names(data)[i], Cindex = cindex_coxboost, Model = 'Ridge + CoxBoost')
  result <- rbind(result, cc)
}




#############################
### 28. ENet + stepcox  #####
###                     #####
#############################
set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Prepare training data
  x1 <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  x2 <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  # Loop over alpha values for Elastic Net
  for (alpha in seq(0, 1, 0.1)) {
    set.seed(seed)
    # Fit Enet model with the current alpha
    fit_enet <- cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
    
    # Predict risk scores using Enet on the validation dataset
    x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation data excluding time and status
    preds_enet <- predict(fit_enet, newx = x_val, s = fit_enet$lambda.min, type = "link")
    rs_enet <- cbind(val_dd[, 1:2], RS = as.vector(preds_enet))  # Combine with time and status
    
    # Fit Stepwise Cox on the training data
    for (direction in c("both", "backward", "forward")) {
      fit_stepcox <- step(coxph(Surv(time, status) ~ ., est_dd), direction = direction, trace = 0)
      
      # Predict risk scores using Stepwise Cox on the validation dataset
      rs_stepcox <- cbind(val_dd[, 1:2], RS = predict(fit_stepcox, type = 'risk', newdata = val_dd))
      
      # Calculate the concordance index (C-index)
      cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs_stepcox))$concordance[1])
      
      # Add results to the result data frame
      cc <- data.frame(ID = names(data)[i],
                       Cindex = cindex,
                       Model = paste0('Enet[α=', alpha, '] + StepCox[', direction, ']'))
      result <- rbind(result, cc)
    }
  }
}





##########################
#### 29. Enet + coxboost
##########################

set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Prepare training data
  x1 <- as.matrix(est_dd[, -c(1, 2)])  # Exclude time and status
  x2 <- Surv(est_dd$time, est_dd$status)  # Survival object
  
  # Loop over alpha values for Elastic Net
  for (alpha in seq(0, 1, 0.1)) {
    set.seed(seed)
    # Fit Enet model with the current alpha
    fit_enet <- cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
    
    # Predict risk scores using Enet on the validation dataset
    x_val <- as.matrix(val_dd[, -c(1, 2)])  # Validation data excluding time and status
    preds_enet <- predict(fit_enet, newx = x_val, s = fit_enet$lambda.min, type = "link")
    rs_enet <- cbind(val_dd[, 1:2], RS = as.vector(preds_enet))  # Combine with time and status
    
    # Fit CoxBoost on the training data
    fit_coxboost <- CoxBoost(time = est_dd$time, status = est_dd$status, 
                             x = as.matrix(est_dd[, -c(1, 2)]), stepno = 100)
    
    # Predict risk scores using CoxBoost on the validation dataset
    preds_coxboost <- predict(fit_coxboost, newdata = x_val, type = "lp")
    rs_coxboost <- cbind(val_dd[, 1:2], RS = as.vector(preds_coxboost))
    
    # Calculate the concordance index (C-index)
    cindex <- as.numeric(summary(coxph(Surv(time, status) ~ RS, rs_coxboost))$concordance[1])
    
    # Add results to the result data frame
    cc <- data.frame(ID = names(data)[i],
                     Cindex = cindex,
                     Model = paste0('Enet[α=', alpha, '] + CoxBoost'))
    result <- rbind(result, cc)
  }
}


#####################
#####  30. stepcox + coxboost

##########################3
set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric and valid
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Iterate over stepwise directions
  for (direction in c("both", "backward", "forward")) {
    # Step 1: Apply Stepwise Cox to the combined training data
    fit_stepcox <- step(coxph(Surv(time, status) ~ ., data = est_dd), direction = direction, trace = 0)
    
    # Step 2: Get selected variables from Stepwise Cox
    selected_vars <- names(coef(fit_stepcox))
    
    # Subset training and validation datasets with selected variables
    train_dd2 <- est_dd[, c("time", "status", selected_vars)]
    val_dd2 <- val_dd[, c("time", "status", selected_vars), drop = FALSE]
    
    # Step 3: Apply CoxBoost to the training data
    fit_coxboost <- CoxBoost(
      time = train_dd2$time, 
      status = train_dd2$status, 
      x = as.matrix(train_dd2[, -c(1, 2)]), 
      stepno = 100
    )
    
    # Step 4: Predict risk scores for the validation dataset
    val_x <- as.matrix(val_dd2[, -c(1, 2), drop = FALSE])  # Exclude time and status
    val_dd2$RS <- as.numeric(predict(fit_coxboost, newdata = val_x, type = "lp"))  # Risk score
    
    # Step 5: Calculate C-index for the validation dataset
    c_index <- as.numeric(summary(coxph(Surv(time, status) ~ RS, data = val_dd2))$concordance[1])
    
    # Store results
    result <- rbind(result, data.frame(
      ID = names(data)[i], 
      Cindex = c_index, 
      Model = paste0("StepCox[", direction, "] + CoxBoost")
    ))
  }
}



#####################
#####  31. stepcox + plsRcox

##########################

set.seed(seed)


# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric and valid
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Iterate over stepwise directions
  for (direction in c("both", "backward", "forward")) {
    # Step 1: Apply Stepwise Cox to the combined training data
    fit_stepcox <- step(coxph(Surv(time, status) ~ ., data = est_dd), direction = direction, trace = 0)
    
    # Step 2: Get selected variables from Stepwise Cox
    selected_vars <- names(coef(fit_stepcox))
    
    # Subset training and validation datasets with selected variables
    train_dd2 <- est_dd[, c("time", "status", selected_vars)]
    val_dd2 <- val_dd[, c("time", "status", selected_vars), drop = FALSE]
    
    # Step 3: Cross-validated plsRcox to determine the optimal number of components
    cv.plsRcox.res <- cv.plsRcox(
      list(x = train_dd2[, -c(1, 2)], time = train_dd2$time, status = train_dd2$status), 
      nt = 10, verbose = FALSE
    )
    nt_optimal <- as.numeric(cv.plsRcox.res[5])  # Optimal number of components
    
    # Step 4: Fit plsRcox with the optimal number of components
    fit_plsRcox <- plsRcox(
      train_dd2[, -c(1, 2)], 
      time = train_dd2$time, 
      event = train_dd2$status, 
      nt = nt_optimal
    )
    
    # Step 5: Predict risk scores for the validation dataset
    val_x <- as.matrix(val_dd2[, -c(1, 2), drop = FALSE])  # Exclude time and status
    val_dd2$RS <- as.numeric(predict(fit_plsRcox, type = "lp", newdata = val_x))  # Risk score
    
    # Step 6: Calculate C-index for the validation dataset
    c_index <- as.numeric(summary(coxph(Surv(time, status) ~ RS, data = val_dd2))$concordance[1])
    
    # Store results
    result <- rbind(result, data.frame(
      ID = names(data)[i], 
      Cindex = c_index, 
      Model = paste0("StepCox[", direction, "] + plsRcox")
    ))
  }
}



  #######################
##### 32, Stepcox + superPC
#############################









##################################
#### 34.StepCox+ gbm         ####
##################################


set.seed(seed)

# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric and valid
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Iterate over stepwise directions
  for (direction in c("both", "backward", "forward")) {
    # Step 1: Apply Stepwise Cox to the combined training data
    fit_stepcox <- step(coxph(Surv(time, status) ~ ., data = est_dd), direction = direction, trace = 0)
    
    # Step 2: Get selected variables from Stepwise Cox
    selected_vars <- names(coef(fit_stepcox))
    
    # Skip if no variables are selected
    if (length(selected_vars) == 0) {
      next
    }
    
    # Subset training and validation datasets with selected variables
    train_dd2 <- est_dd[, c("time", "status", selected_vars), drop = FALSE]
    val_dd2 <- val_dd[, c("time", "status", selected_vars), drop = FALSE]
    
    # Skip if validation data has no rows or columns
    if (nrow(val_dd2) == 0 || ncol(val_dd2) <= 2) {
      next
    }
    
    # Step 3: Train the GBM model on the training data
    set.seed(seed)
    
    # Initial GBM training to determine the best number of trees
    fit_gbm <- gbm(
      formula = Surv(time, status) ~ ., 
      data = train_dd2, 
      distribution = 'coxph', 
      n.trees = 10000, 
      interaction.depth = 3, 
      n.minobsinnode = 10, 
      shrinkage = 0.001, 
      cv.folds = 10, 
      n.cores = 6,
      verbose = FALSE
    )
    
    # Find the index for the number of trees with minimum CV error
    best_trees <- which.min(fit_gbm$cv.error)
    
    # Retrain the GBM model using the best number of trees
    set.seed(seed)
    fit_gbm <- gbm(
      formula = Surv(time, status) ~ ., 
      data = train_dd2, 
      distribution = 'coxph', 
      n.trees = best_trees, 
      interaction.depth = 3, 
      n.minobsinnode = 10, 
      shrinkage = 0.001, 
      cv.folds = 10, 
      n.cores = 6,
      verbose = FALSE
    )
    
    # Step 4: Predict risk scores for the validation dataset
    preds <- predict(fit_gbm, val_dd2, n.trees = best_trees, type = 'link')
    val_dd2$RS <- as.numeric(preds)
    
    # Step 5: Calculate C-index for the validation dataset
    c_index <- as.numeric(summary(coxph(Surv(time, status) ~ RS, data = val_dd2))$concordance[1])
    
    # Store results
    result <- rbind(result, data.frame(
      ID = names(data)[i], 
      Cindex = c_index, 
      Model = paste0("StepCox[", direction, "] + GBM")
    ))
  }
}




##################################
#### 35.StepCox+survival-SVM ####
##################################

for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric and valid
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Iterate over stepwise directions
  for (direction in c("both", "backward", "forward")) {
    # Step 1: Apply Stepwise Cox to the combined training data
    fit_stepcox <- step(coxph(Surv(time, status) ~ ., data = est_dd), direction = direction, trace = 0)
    
    # Step 2: Get selected variables from Stepwise Cox
    selected_vars <- names(coef(fit_stepcox))
    
    # Skip if no variables are selected
    if (length(selected_vars) == 0) {
      next
    }
    
    # Subset training and validation datasets with selected variables
    train_dd2 <- est_dd[, c("time", "status", selected_vars), drop = FALSE]
    val_dd2 <- val_dd[, c("time", "status", selected_vars), drop = FALSE]
    
    # Skip if validation data has no rows or columns
    if (nrow(val_dd2) == 0 || ncol(val_dd2) <= 2) {
      next
    }
    
    # Step 3: Train the survival-SVM model on the training data
    set.seed(seed)
    fit_survivalsvm <- survivalsvm(
      Surv(time, status) ~ ., 
      data = train_dd2, 
      gamma.mu = 1, 
      type = "regression"
    )
    
    # Step 4: Predict risk scores for the validation dataset
    preds <- predict(fit_survivalsvm, val_dd2)$predicted
    val_dd2$RS <- as.numeric(preds)
    
    # Step 5: Calculate C-index for the validation dataset
    c_index <- as.numeric(summary(coxph(Surv(time, status) ~ RS, data = val_dd2))$concordance[1])
    
    # Store results
    result <- rbind(result, data.frame(
      ID = names(data)[i], 
      Cindex = c_index, 
      Model = paste0("StepCox[", direction, "] + survival-SVM")
    ))
  }
}

##################################
#### 34.Coxboost + plsRcox ####
##################################
# Fit CoxBoost Model
# Perform LOOCV
library(survival)
library(CoxBoost)
library(plsRcox)
library(glmnet)


# Perform LOOCV
for (i in seq_along(data)) {
  # Leave one dataset out for validation
  val_dd <- data[[i]]
  train_dd_list <- data[-i]  # All other datasets for training
  
  # Combine training datasets into one
  est_dd <- do.call(rbind, train_dd_list)
  
  # Ensure data is numeric and valid
  est_dd <- convert_to_numeric(est_dd)
  val_dd <- convert_to_numeric(val_dd)
  est_dd$time[est_dd$time <= 0] <- 0.1  # Ensure valid survival times
  
  # Step 1: Fit the LASSO model for feature selection
  set.seed(seed)
  x_est <- as.matrix(est_dd[, -c(1, 2)])  # Exclude 'time' and 'status'
  y_est <- Surv(est_dd$time, est_dd$status)
  cvfit <- cv.glmnet(x_est, y_est, family = "cox", alpha = 1, nfolds = 10)
  
  # Selected features
  rid <- rownames(coef(cvfit, s = "lambda.min"))[which(coef(cvfit, s = "lambda.min") != 0)]
  
  # Skip iteration if no variables are selected
  if (length(rid) == 0) next
  
  # Step 2: Subset the data with selected features
  est_dd2 <- est_dd[, c("time", "status", rid), drop = FALSE]
  val_dd2 <- val_dd[, c("time", "status", rid), drop = FALSE]
  
  # Skip if validation data or training data has invalid dimensions
  if (ncol(est_dd2) <= 2 || ncol(val_dd2) <= 2) next
  
  # Step 3: Fit the plsRcox model
  set.seed(seed)
  cv.plsRcox.res <- cv.plsRcox(
    list(x = est_dd2[, -c(1, 2)], time = est_dd2$time, status = est_dd2$status),
    nt = 10, verbose = FALSE
  )
  fit_plsrc <- plsRcox(
    est_dd2[, -c(1, 2)], 
    time = est_dd2$time, 
    event = est_dd2$status, 
    nt = as.numeric(cv.plsRcox.res[5])
  )
  
  # Step 4: Fit the CoxBoost model
  pen <- optimCoxBoostPenalty(
    est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]),
    trace = TRUE, start.penalty = 500, parallel = TRUE
  )
  cv.res <- cv.CoxBoost(
    est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]),
    maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
  )
  fit_cx <- CoxBoost(
    est_dd2$time, est_dd2$status, as.matrix(est_dd2[, -c(1, 2)]),
    stepno = cv.res$optimal.step, penalty = pen$penalty
  )
  
  # Step 5: Predict risk scores for the validation dataset
  rs_combined <- tryCatch({
    pred_cx <- as.numeric(predict(fit_cx, newdata = val_dd2[, -c(1, 2)], newtime = val_dd2$time, newstatus = val_dd2$status, type = "lp"))
    pred_pls <- as.numeric(predict(fit_plsrc, type = "lp", newdata = val_dd2[, -c(1, 2)]))
    RS <- (pred_cx + pred_pls) / 2  # Combine predictions (e.g., averaging)
    cbind(val_dd2[, 1:2], RS)  # Combine with time and status
  }, error = function(e) NULL)  # Skip in case of errors
  
  # Skip if predictions are invalid
  if (is.null(rs_combined)) next
  
  # Step 6: Calculate C-index for the validation dataset
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(time, status) ~ RS, data = rs_combined))$concordance[1]))
  cc$ID <- names(data)[i]
  cc$Model <- 'CoxBoost + plsRcox'
  
  # Combine with the result object
  result <- rbind(result, cc)
}








##### plot heatmap of final

result2 <- result
result2$Model <- gsub('α','a',result2$Model)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(cowplot)
library(gridExtra)
range(result2$Cindex)


result2 %>%
  
  ggplot(aes(Cindex, reorder(Model, Cindex))) +
  geom_bar(width = 0.7, stat = 'summary', fun = 'mean', fill = 'lightblue') +
  geom_text(aes(label = round(Cindex, 2)), 
              # Adjust text position based on the value
            
            position = position_stack(vjust = 0.7), color = 'gray10') + 
  scale_x_continuous(limits = c(0, 0.7)) + # Adding values inside the bars
  theme_classic() +
  labs(y = NULL)

dd <- result2%>%
  #filter(ID!='TCGA')%>%
  group_by(Model)%>%
  summarise(Cindex=mean(Cindex))

bar_plot <- dd %>%
  ggplot(aes(Cindex, reorder(Model, Cindex))) +
  geom_bar(width = 0.8, stat = 'identity', fill = 'lightblue') +
  geom_text(aes(label = round(Cindex, 2)), size = 3.5,
            position = position_stack(vjust = 0.7), color = 'gray10') + 
  scale_x_continuous(limits = c(0, 0.8)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(), 
   # axis.text.y = element_blank(),# Remove y-axis title
    # Remove y-axis text
    axis.ticks.y = element_blank()   # Remove y-axis ticks
  ) +
  labs(y = NULL)  # Remove default y-axis
bar_plot

# Extract the order of Model from the ggplot
model_order <- dd$Model[order(-dd$Cindex)]  # Order based on Cindex in descending order

dd2 <- pivot_wider(result2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()

# Create dd2 and dd3 as before
dd3 <- column_to_rownames(dd2, var = 'Model')

# Reorder dd3 based on the extracted model order
dd3 <- as.matrix(dd3[model_order, , drop = FALSE])  # Use drop=FALSE to avoid dropping rows if levels don't match


# Load necessary library
library(pheatmap)
library(colorRamp2)
library(RColorBrewer)
groups <- c("TCGA","GSE26253", "GSE15459")  # Example group information

annotation_col <- data.frame(
  Group = groups  # Custom grouping information
)
rownames(annotation_col) <- colnames(dd3)  # Set the rownames to match the column names of dd3
ann_colors = list(
  Group = c(
    TCGA = "#a6cee3",      # Color for TCGA
    GSE26253 = "#b2df8a",  # Color for GSE84437
    GSE15459 = "#fb9a99"  # Color for GSE84433
    #GSE84433 = "orange3"   # Color for GSE84426
  )
)
# Create a vector for gaps after every column
gaps_col <- seq(1, ncol(dd3) - 1)  # Creates a sequence of indices for gaps
custom_colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100) # Example using a blue gradient
breaks <- seq(0.3, 0.88, length.out = 101)  # 101 breaks from 0 to 1

# Create the heatmap with column annotation
pheatmap_grob <- pheatmap(dd3,
         labels_row = rownames(dd3),
         display_numbers = T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = F,
         cluster_rows  = F,
         border_color = 'gray30',
         show_colnames = FALSE,
         show_rownames = FALSE,
         
         gaps_col = gaps_col,
         color = custom_colors,
         breaks = breaks)

pheatmap_grob


###### choose the best model
####  choose the features

x1 <- as.matrix(est_dd[, selected_gene])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))

set.seed(seed)
fit = cv.glmnet(x1, x2, family = "cox", alpha = 0.1, nfolds = 10)

rid <- rownames(coef(fit, s = "lambda.min"))[which(coef(fit, s = "lambda.min") != 0)]  # Selected features
est_dd2 <- est_dd[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x) { x[, c('time', 'status', rid)] })

# Fit RSF model using selected features from Enet
set.seed(seed)
fit_rsf <- rfsrc(Surv(time,status) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize,
                 splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE)

  
  




# Extract variable importance
# Extract the importance scores from the RSF model
importance_scores <- fit_rsf$importance

# Create a data frame of features and their corresponding importance scores
important_features <- data.frame(
  Feature = names(importance_scores),
  Importance = importance_scores
)

# View the important features
print(important_features)

# Normalize importance scores to sum to 1 (optional)
importance_weights <- importance_scores / sum(importance_scores)

# Extract only the positive importance scores
positive_genes <- names(importance_scores[importance_scores > 0])  #### this is the final gene list of our paper. Use this for single cell validation



positive_weights <- importance_scores[positive_genes]

data<-readRDS('subset_input_ML_dataframes_after_univariate.RDS') # obtain from "ML_input_preparation.R"


# Calculate risk scores and add the 'Risk' column
# Construct the formula for the Cox model using the positive genes
# Loop through each dataset in the 'data' list
risk_scores <- lapply(data, function(dataset) {
  # Construct the formula for the Cox model using the positive genes
  formula <- as.formula(paste("Surv(as.numeric(time), as.numeric(status)) ~", paste(positive_genes, collapse = "+")))
  
  # Fit the Cox model to the current dataset
  fit_multi <- coxph(formula, data = dataset)
  
  # Extract the coefficients from the fitted Cox model
  coefficients <- coef(fit_multi)
  
  # Select the features corresponding to positive genes
  features <- as.matrix(dataset[, positive_genes])
  
  # Calculate the risk score for each sample using the Cox model coefficients
  risk_score <- apply(features, 1, function(row) sum(row * coefficients))
  
  # Determine the median risk score for this dataset
  median_risk <- median(risk_score)
  
  # Assign 'high risk' or 'low risk' based on the median
  risk_label <- ifelse(risk_score > median_risk, "high risk", "low risk")
  
  # Add 'RiskScore' and 'Risk' columns to the original data
  cbind(dataset[, 1:2], RiskScore = risk_score, Risk = risk_label)
})

# Now `risk_scores` will contain the risk scores and risk categories for each sample in each dataset


# Now `risk_scores` will contain the risk scores and risk categories for each sample in each dataset

# Check one example output
print(risk_scores[[1]])


# Load required package
# Load required packages
library(ggpubr)
library(dplyr)

# Extract importance scores and convert to a data frame
importance_df <- data.frame(Feature = names(importance_scores), 
                            Importance = importance_scores)

# Sort features by importance value
importance_df <- importance_df %>%
  arrange(Importance) %>%
  mutate(Direction = ifelse(Importance > 0, "Positive", "Negative"))  # Add direction

# Ensure 'Direction' is treated as a factor
importance_df$Direction <- factor(importance_df$Direction, levels = c("Positive", "Negative"))

# Create the lollipop plot
ggdotchart(importance_df, 
           x = "Feature", 
           y = "Importance", 
           color = "Direction",  # Use the new Direction column for color
           sorting = "descending",
           rotate = TRUE, 
           add = "segments",  # Add lines for lollipop effect
           add.params = list(color = "lightgray", size = 0.5), 
           dot.size = 5, 
           ggtheme = theme_classic2()) +
  scale_color_manual(values = c("Positive" = "forestgreen", "Negative" = "red")) +
  theme(legend.text = element_text(size = 15),      # Adjust the size of legend text
        axis.text = element_text(size = 20),        # Adjust the size of axis labels
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 15,face = 'bold'),# Hide x-axis tick labels
        # Adjust the size of axis titles
        strip.text = element_text(size = 15),       # Adjust the size of strip (facet) labels
        plot.caption = element_text(size = 15),     # Adjust the size of the plot caption
        axis.text.x = element_text(size = 15),      # Adjust the size of x-axis tick labels
        axis.text.y = element_text(size = 12),      # Adjust the size of y-axis tick labels
        strip.text.x = element_text(size = 15)) + 
  labs(title = "RSF + Stepcox[forward]", 
       x = "Feature", 
       y = "Importance")

#############################################
#############################################
#####                    ####################
##### SURVIVAL ANALYSIS  ####################
#####                    ####################
#############################################
library(survival)
library(glmnet)
library("survminer")
surv_data<-risk_scores$tcga
surv_data$time<-as.numeric(surv_data$time)
surv_data$status<-as.numeric(surv_data$status)
# find optimal cutpoint for risk score
surv_cut <- surv_cutpoint(surv_data, time = 'time', event = 'status',variables = 'RiskScore')
cutoff<-surv_cut$cutpoint$cutpoint
surv_data$risk<-ifelse(surv_data$RiskScore <= cutoff, "low risk", "high risk")

survival_fit <- survfit(Surv(time, status) ~ Risk, data = surv_data)
ggsurvplot(survival_fit, data = surv_data, pval = TRUE,
           pval.cExpression.Levelrd = c(50, 1),
           pval.size=4,legend="top",
           palette = c("Red","green4"),
           risk.table = F,font.main = c(12, "bold"),
           legend.title = "group",
           font.legend = c(14, "bold"),
           #title='EIF2 signaling (CGGA)',
           font.x = c(14, "bold"),
           font.y = c(14, "bold"),
           legend.labs = c("High risk","Low risk"),
           font.tickslab = c(12, "plain"),xlab="Time")


#Extract the survival probabilities and event status for survivalROC:
library(survivalROC)
library(timeROC)
ROC <- timeROC(T=surv_data$time,   
               delta=surv_data$status,   
               marker=surv_data$RiskScore,   
               cause=1,               
               weighting="marginal",   
               times=c(1*12, 3*12,5*12),       
               iid=T,
               ROC=T)
# Set the size of axis labels and ticks
x_axis_size <- 1.5  # Change this to increase x-axis size
y_axis_size <- 1.5  # Change this to increase y-axis size

# Plot the ROC curve for different time points
plot(ROC, 
     time=1*12, 
     col="red", 
     lwd=2,
     title = "")

plot(ROC,
     time=3*12, 
     col="blue", 
     add=TRUE, 
     lwd=2)

plot(ROC,
     time=5*12, 
     col="orange", 
     add=TRUE, 
     lwd=2)

# Add the legend
legend("bottomright",
       c(paste0("AUC at 1 year: ", round(ROC[["AUC"]][1], 2)), 
         paste0("AUC at 3 year: ", round(ROC[["AUC"]][2], 2)), 
         paste0("AUC at 5 year: ", round(ROC[["AUC"]][3], 2))),
       col=c("red", "blue", "orange"),
       lty=1, 
       lwd=1, 
       bty = "n")


##### Cach 2
time_ROC_df <- data.frame(
  TP_1year = ROC$TP[, 1],
  FP_1year = ROC$FP[, 1],
  TP_3year = ROC$TP[, 2],
  FP_3year = ROC$FP[, 2],
  TP_5year = ROC$TP[, 3],
  FP_5year = ROC$FP[, 3])

library(ggplot2)

ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.2f", ROC$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.2f", ROC$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.2f", ROC$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

################## save Risk score list for further clinical related analysis

#saveRDS(risk_scores,file = 'risk_score_list_for_clinical.rds')  ## This is used for forestplot and nomogram, refer to "clinical_analysis_stad.R"

