##' @title  Synthetic Control Groups -- Functions
##' @author Pascal Amiet
##' @date   12.08.2023


# Packages ----

.pkg <- c('dplyr', 'glmnet', 'mgcv', 'rpart', 'randomForest', 'mgcv', 'earth',
          'progress','tidyr','cubature','sfsmisc','ks','e1071','modelr','caret',
          'tidyverse')

for (.p in .pkg) {
  if (!require(.p, character.only = TRUE)) {
    install.packages(.p)
  }
}

# Colors ----

# Define the color palette for plots

.c <- c('aquamarine2','aquamarine3','aquamarine4','skyblue4','skyblue3',
        'skyblue2','violetred4','violetred3','violetred2','mediumpurple2',
        'mediumpurple3','mediumpurple4')

.colorShuffle <- c(2, 5, 8, 10, 3, 6, 9, 12, 1, 4, 7, 11)

.c2 <- .c[.colorShuffle]

# Synthetic Control Estimation ----

##' @description:
##' This function computes the synthetic control group given some
##' input data and a given method.
##' @param data: The first input is the data matrix that is to be
##' evaluated. It is important that the rows are the periods and
##' the columns are the different groups/covariates.
##' @param treated: Here, we define the column which will be
##' treated, i.e. the treatment group. This is either a string or
##' the column number.
##' @param t: Here we indicate the treatment period. This has to
##' correspond to the specific index.
##' @param method: The method defines which statistical model is
##' to be used for the computations. One can choose between `lasso`,
##' `ols`, `gam`.
##' @param significance: Indicates whether or not to calculate the significance
##' level of the estimate. The default is `significance = TRUE`.
##' @return The function returns a list with the estimation of the
##' synthetic control group, as well as some error statistics, significance,
##' and the absolute treatment effect.

synthCtrl <- function(data, treated, t, method, significance = TRUE) {
  
  # Error handling 1
  if (! (method %in% c("ols", "lasso", "gam"))) {
    cat("\n\nPlease choose a valid method.\n\n")
  } else if (! (t %in% rownames(data))) {
    cat("\n\nPlease choose a valid intervention period.\n\n")
  } else if (! (treated %in% names(data))) {
    cat("\n\nPlease choose a valid treatment group.\n\n")
  } else {

    # Train/test split
    .train <- data[1:which(rownames(data) == t), ]
    .test  <- data[(which(rownames(data) == t) + 1):nrow(data), ]

    .y_train <- .train %>% dplyr::pull(treated)
    .x_train <- .train %>%
        dplyr::select(-dplyr::all_of(treated)) %>%
        data.matrix()

    .y_test <- .test %>% dplyr::pull(treated)
    .x_test <- .test %>%
        dplyr::select(-dplyr::all_of(treated)) %>%
        data.matrix
    
    # Error handling 2
    if (nrow(.x_train) < ncol(.x_train) &
        ! (method %in% c('lasso'))) {
      cat("\n\nThe matrix does not have full rank. Please choose a method that accounts for that, i.e. a shrinkage estimator.\n\n")
    } else {

      # Fitting
      if (method == "lasso") {
        
          .cv <- glmnet::cv.glmnet(.x_train, .y_train, alpha = 1,
                                        intercept = FALSE)
          .lambda <- .cv$lambda.min
          .model <- glmnet::glmnet(.x_train, .y_train, alpha = 1,
                                   lambda = .lambda)
          .wgt <- glmnet::glmnet(.x_train, .y_train, alpha = 1,
                                 lambda = .lambda) %>%
                                 coef()
          .wgt <- data.matrix(.wgt) %>%
            as.data.frame() %>%
            dplyr::rename(coef = s0) %>%
            dplyr::slice(2:dplyr::n()) %>%
            t()
          
      } else if (method == "ols") {
        
        .model <- lm(.y_train ~ .x_train)
        
      } else if (method == "gam") {
        
        .model <- gam(.y_train ~ .x_train)
        
      }
      
      # Adjust data to be predicted for different methods
      if (method %in% c("ols")) {
        .topredict <- as.data.frame(rbind(.x_train, .x_test))
      } else {
        .topredict <- rbind(.x_train, .x_test)
      }
  
      # Predict
      .synth <- predict(.model, .topredict) %>%
          as.data.frame() %>%
          dplyr::pull(1)
      
      names(.synth) <- rownames(data)
  
      # Stats
      .stats <- getStats(.x_train, .y_train, .x_test, .y_test, .synth)
      
      # Treatment effect
      .eff <- treatEff(c(.y_train, .y_test), .synth, t)
      
      # Significance
      if (significance == TRUE) {
        .sig <- testStat(data, treated, t, method, plot = FALSE) 
      }
  
      # Return
      .ls <- vector(mode = "list")
      .ls$synth <- .synth
      .ls$stats <- .stats
      if (significance == TRUE) .ls$significance <- .sig
      .ls$effect <- .eff
      return(.ls)
      
    } # Error handling 2
  } # Error handling 1
}


##' @description
##' This function calculates all the necessary numbers to compare methods.
##' @param xtrain: The training covariates.
##' @param ytrain: The training response variable.
##' @param xtest: The testing covariates.
##' @param ytest: The testing response variable.
##' @param synth: The predictions for the synthetic control group.
##' @return The function returns a matrix containing MSE, RMSE, MAE, and MAPE
##' for both in-sample predictions and out-of-sample predictions.

getStats <- function(xtrain, ytrain, xtest, ytest, synth) {
  
  .stats <- matrix(NA, nrow = 2, ncol = 5)
  rownames(.stats) <- c("Prediction", "In-Sample")
  colnames(.stats) <- c("MSE", "RMSE", "MAE", "MAPE", "n")
  .stats[, "n"] <- c(nrow(xtest), nrow(xtrain))
  .stats[, "MSE"] <- c(mean((ytest - synth[(length(ytrain) + 1):length(synth)])^2), # nolint: line_length_linter.
                       mean((ytrain - synth[1:length(ytrain)])^2))
  .stats[, "RMSE"] <- sqrt(.stats[, "MSE"])
  .stats[, "MAE"] <- c(mean(abs(ytest - synth[(length(ytrain) + 1):length(synth)])), # nolint: line_length_linter.
                       mean(abs(ytrain - synth[1:length(ytrain)])))
  .stats[, "MAPE"] <- c(mean(abs(ytest - synth[(length(ytrain) + 1):length(synth)]) / ytest), # nolint: line_length_linter.
                        mean(abs(ytrain - synth[1:length(ytrain)]) / ytrain)) # nolint: line_length_linter.
  .stats <- round(.stats, digits = 6)
  
  return(.stats)
  
}


##' @description
##' This function computes the test statistic, i.e. the significance level, of
##' for the synthetic control group. It does so by using a placebo test as 
##' suggested by the literature.
##' @param data: The first input is the data matrix that is to be
##' evaluated. It is important that the rows are the periods and
##' the columns are the different groups/covariates.
##' @param treated: Here, we define the column which will be
##' treated, i.e. the treatment group. This is either a string or
##' the column number.
##' @param t: Here we indicate the treatment period. This has to
##' correspond to the specific index.
##' @param method: The method defines which statistical model is
##' to be used for the computations. One can choose between `lasso`,
##' `ols`, `gam`.
##' @param plot: If true, it will return a histogram of the estimated gaps
##' to illustrate the significance. The default is that `plot = TRUE`.
##' @return The function returns the significance level provided by the placebo
##' test.

testStat <- function(data, treated, t, method, plot = TRUE) {
  # Gap for the actual synthetic control group
  .gap <- synthCtrl(data, treated, t, method, significance = FALSE)$stats[1,1]
  .gaps <- c()
  .others <- names(data)[-which(names(data) == treated)]
  for (.n in .others) {
    .g <- synthCtrl(data, .n, t, method, significance = FALSE)$stats[1,1]
    .gaps <- c(.gaps, .g)
  }
  # Plot
  if (plot == TRUE) {
    hist(c(.gaps, .gap), main = "Histogram of MSPEs", xlab = "MSPE", 
         probability = TRUE)
    abline(v = .gap, lty = "dashed") 
  }
  .test <- sum(c(.gaps, .gap) >= .gap) / length(c(.gaps, .g))
  return(.test)
}


# Treatment Effect ----

##' @description
##' This function calculates the treatment effect given the treatment and the
##' (synthetic) control group.
##' @param treatment: Time series of the treatment group.
##' @param control: Time series of the (synthetic) control group.
##' @param t: Treatment period.
##' @return The function returns the treatment effect.

treatEff <- function(treatment, control, t) {
  t <- as.character(t)
  if (length(treatment) != length(control)) {
    ##' Error handling
    cat("\n\nThe two time series are of different length.\n\n")
  } else if (length(treatment) == 2 & length(control) == 2) {
    ##' Closed form solution
    .eff = (treatment[2] - treatment[1]) - (control[2] - control[1])
    return(.eff)
  } else {
    if (is.null(names(control))) {
      .timeFE <- rep(as.character(1:length(control)), 2)
    } else {
      .timeFE <- rep(names(control), 2)
    }
    .unitFE <- c(rep("treatment", length(treatment)), 
                 rep("control", length(control)))
    .i <- as.numeric(ifelse(is.null(names(control)), t, which(names(control) == t)))
    .treat <- c(rep(0, (.i-1)), rep(1, length(treatment) - .i + 1), 
                rep(0, length(control)))
    .X <- data.frame(outcome = c(treatment, control), time = .timeFE, 
                    unit = .unitFE, treatment = .treat)
    .eff <- lm(outcome ~ treatment + as.factor(time) + as.factor(unit), data = .X) %>% 
      coef() %>% .[[2]]
    return(.eff)
  }
}


# Statistical Methods ----

##' @description
##' This function fits a bagged tree model on the data.
##' @param data: The data to be fitted.
##' @param y: The variable of interested.
##' @param B: The number of trees to be fitted, by default `B = 20`.
##' @param contr: Controls for the splits when estimating the single trees,
##' given as the standard `rpart.control(...)` element from the `rpart` package.
##' @return The function returns the fitted bagged tree model.

baggedTree <- function(data, y, B=20, 
                       contr=rpart.control(minsplit=2, cp=0, minbucket=1)){
  # Create list for storage
  .l <- list()
  for (b in 1:B){
    # Bootstrap
    .s <- data[sample(nrow(data), nrow(data), replace=T), ]
    # Train model
    .f <- as.formula(paste(y,'~ .'))
    .t <- rpart(.f, data = .s, control = contr)
    .l[[b]] <- .t
  }
  # Return list
  return(.l)
}


##' @description
##' This function is used for prediction of the bagged tree model.
##' @param model: The bagged tree model that one gets as an output from the 
##' `baggedTree()` function.
##' @param newdata: The data (covariates) used to predict.
##' @return The function returns a vector of predicted variables.

predict.baggedTree <- function(model, newdata){
  .p <- sapply(model, function(model) predict(model, newdata))
  .m <- apply(.p, 1, mean)
  return(.m)
}


##' @description
##' This function fits a boosted tree model on the data.
##' @param data: The data to be fitted.
##' @param y: The variable of interest.
##' @param B: The number of trees to be fitted, by default `B = 20`.
##' @param rate: The learning rate, by default `rate = 0.1`.
##' @param info: Whether or not information should be displayed during training,
##' by default `infor = TRUE`.
##' @param contr: Controls for the splits when estimating the single trees,
##' given as the standard `rpart.control(...)` element from the `rpart` package.
##' @return The function returns the fitted boosted tree model.

boostedTree <- function(data, y, B=20, rate=0.1, info=T,
                        contr=rpart.control(minsplit=2, cp=0, minbucket=1)){
  # Create list for storage
  .l <- list()
  # The residuals are initially equal to the data
  .res <- data
  attr(.l, 'rate') <- rate
  for (b in 1:B){
    .f <- as.formula(paste(y,'~ .'))
    .t <- rpart(.f, data=.res, control=contr)
    .l[[b]] <- .t
    .res[[y]] <- .res[[y]] - (rate * predict(.t, .res))
    if (info==T) print(paste(b,'/',B,'-- RMSE =',sqrt(mean((.res[[y]])^2))))
  }
  # Return list
  return(.l)
}


##' @description
##' This function is used for prediction of the boosted tree model.
##' @param model: The bagged tree model that one gets as an output from the 
##' `boostedTree()` function.
##' @param newdata: The data (covariates) used to predict.
##' @return The function returns a vector of predicted variables.

predict.boostedTree <- function(model, newdata){
  .r <- attributes(model)$rate
  .p <- sapply(model, function(model) .r * predict(model, newdata))
  .s <- apply(.p, 1, sum)
  return(.s)
}


##' @description
##' This function uses distance measures to estimate the weights of different 
##' control groups in a non-parametric way.
##' @param X: The data frame used for fitting the weights.
##' @param k: The variable of interest as a string.
##' @param id: The variable identifying the different groups as a string.
##' @param id.interest: The group we want to create a counterfactual for as a
##' string.
##' @param periods: The variable name as a string of the different periods.
##' @param intervention: The period in which the intervention took place, as an
##' integer.
##' @param controls: The control variables as strings in a vector.
##' @param distance: The distance measure to compute the disimilarity between
##' observations. The options are: `Norm`, `FDS`.
##' @param method: The method of selecting a single weight per control group.
##' This can be either `'Aggregate'` or `'Norm'`, where aggregating means 
##' taking the mean. By default, `method = 'Norm'`.
##' @param norm.wgt: Boolean value indicating whether or not to normalize the
##' estiamted weights, such that the sum of all weights is equal to 1. The 
##' default is set to `norm.wgt = FALSE`.
##' @param true_cf: A vector containing the true counterfactual, if that is
##' known. This will allow to calculate the MSPE. This does not necessairily 
##' have to be provided.
##' @param standardize: A boolean value indicating whether or not to standardize
##' the control variables, i.e. setting mean to 0 and standard deviation to 1.
##' By default, this is not done, i.e. `standardize = FALSE`.
##' @return The function returns a list. The list contains the dissimilarity 
##' matrix, the estimated weights, predictions, as well as performance measures.

dissimilarityMatrix <- function(X, k, id, id.interest, periods = NULL,
                                intervention = NULL, controls = NULL, 
                                distance = 'Norm', method = 'Norm', norm.wgt = F, 
                                true_cf = NULL, standardize = F){
  
  # Controls
  if (is.null(controls)){ controls <- c('x_ctrl'); X <- X %>% mutate(x_ctrl=X[[k]]) }
  
  # Standardize 
  if (standardize){
    X <- X %>% mutate_at(controls, ~(scale(.) %>% as.vector))
  }
  
  # Periods and intervention
  if (is.null(periods)){
    .n <- nrow(X[X[[id]] == id.interest,]); periods = 'x_period'
    X <- X %>% group_by_(.dots = id) %>% mutate(x_period = (1:.n)) %>% 
      as.data.frame()
  }
  if (is.null(intervention)) intervention <- max(X[X[[id]] == id.interest,periods])
  
  # Distance measure
  if (distance == 'Euclidean'){
    .d <- function(u, v){ norm(u - v, type = '2') }
  } else if (distance == 'Manhattan'){
    .d <- function(u, v){ norm(u - v, type = '1') }
  } else if (distance == 'Correlation'){
    .d <- function(u, v){ cor(u, v) }
  }
  
  # Calculate weights
  .units <- unique(X[[id]])[-which(unique(X[[id]]) == id.interest)]
  .mat <- matrix(NA, length(controls), intervention - min(X[[periods]]) + 1)
  rownames(.mat) <- .units; colnames(.mat) <- min(X[[periods]]):intervention
  
  for (.u in .units){
    .IN <- X[X[[id]] == id.interest,]
    .CO <- X[X[[id]] == .u,]
    
    if (!is.null(dim(.IN[,controls]))){
      
      .in <- apply(.IN[,controls], 1, function(x) norm(x ,type = '2'))
      .co <- apply(.CO[,controls], 1, function(x) norm(x ,type = '2'))
      for (.i in 1:length(.in)){
        .dist <- .d(.in[.i], .co[.i])
        if (.i == 1){ .vec <- c(.dist) } else { .vec <- c(.vec, .dist) }
      }
      .mat[toString(.u), ] <- .vec[1:ncol(.mat)]
        
    } else {
      
      .in <- .IN[,controls]; .co <- .CO[,controls]
      for (.i in 1:length(.in)){
        .dist <- .d(.in[.i], .co[.i])
        if (.i == 1){ .vec <- c(.dist) } else { .vec <- c(.vec, .dist) }
      }
      .mat[toString(.u), ] <- .vec[1:ncol(.mat)]
      
    }
  }
  
  # Weights
  if (distance != 'Correlation'){
    if (method == 'Aggregate') .wgt <- apply(.mat^(-1), 1, mean)
    if (method == 'Norm') .wgt <- apply(.mat^(-1), 1, function(x) norm(x, type = '2'))
  } else {
    if (method == 'Aggregate') .wgt <- apply(.mat, 1, mean)
    if (method == 'Norm') .wgt <- apply(.mat, 1, function(x) norm(x, type = '2'))
  }
  
  # Normalize weights
  if (norm.wgt) .wgt <- .wgt / sum(.wgt)
  
  # Separate predictors and outcome variable
  .X <- X[X[[id]] != id.interest, c(id,periods,k)] %>% 
    pivot_wider(names_from = id, values_from = k) %>% 
    select(-all_of(c(periods))) %>% 
    as.matrix()
  
  .Y <- X[X[[id]] == id.interest, k]
  
  # Perfromance metrics
  .preds <- as.vector(.X %*% .wgt); names(.preds) <- unique(X[[periods]])
  .idx <- which(d.in[[periods]] > intervention)
  if (intervention < max(X[X[[id]] == id.interest,periods])){
    .mse <- mean((.preds[-.idx] - .IN[[k]][-.idx])^2)
  } else {
    .mse <- mean((.preds - .IN[[k]])^2)
  }
  if (intervention < max(X[X[[id]] == id.interest,periods]) &
      !is.null(true_cf)){
    .mspe <- mean((.preds[.idx] - true_cf[.idx])^2)
  } else {
    .mspe <- NULL
  }
  
  # Return results
  .l <- list(method = method, matrix = .mat, predictions = .preds, wgt = .wgt, 
             mse = .mse)
  if(!is.null(.mspe)) .l$mspe <- .mspe
  return(.l)
  
}


##' @description
##' This function uses kernels to estimate the weights of different control
##' groups in a non-parametric way. It uses only univariate kernels. 
##' @param X: The data frame used for fitting the weights.
##' @param k: The variable of interest as a string.
##' @param id: The variable identifying the different groups as a string.
##' @param id.interest: The group we want to create a counterfactual for as a
##' string.
##' @param periods: The variable name as a string of the different periods.
##' @param intervention: The period in which the intervention took place, as an
##' integer.
##' @param controls: The control variables as strings in a vector.
##' @param method: The method of selecting a single weight per control group.
##' This can be either `'Aggregate'` or `'Norm'`, where aggregating means 
##' taking the mean. By default, `method = 'Norm'`.
##' @param h: The bandwidth used for the kernel. If set to `NULL`, which is the
##' default, then the function will choose the optimal bandwidth based on 
##' cross validation.
##' @param kernel: The kernel used. The default is `kernel = "Gaussian"`, but
##' one can also choose to instead use `Epanechnikov`, `Triangular`, or 
##' `Tri-cube`.
##' @param norm.wgt: Boolean value indicating whether or not to normalize the
##' estiamted weights, such that the sum of all weights is equal to 1. The 
##' default is set to `norm.wgt = FALSE`.
##' @param true_cf: A vector containing the true counterfactual, if that is
##' known. This will allow to calculate the MSPE. This does not necessairily 
##' have to be provided.
##' @param standardize: A boolean value indicating whether or not to standardize
##' the control variables, i.e. setting mean to 0 and standard deviation to 1.
##' By default, this is not done, i.e. `standardize = FALSE`.
##' @return The function returns a list. The first element of that list is some
##' specifications about the method. The second element are the estimated 
##' weights.

uniKernel <- function(X, k, id, id.interest, periods = NULL, 
                      intervention = NULL, controls = NULL, h = NULL, 
                      kernel = 'Gaussian', method = 'Norm', norm.wgt = F, 
                      true_cf = NULL, standardize = F){
  
  # Controls
  if (is.null(controls)){ controls <- c('x_ctrl'); X <- X %>% mutate(x_ctrl=X[[k]]) }
  
  # Standardize 
  if (standardize){
    X <- X %>% mutate_at(controls, ~(scale(.) %>% as.vector))
  }
  
  # Periods and intervention
  if (is.null(periods)){
    .n <- nrow(X[X[[id]] == id.interest,]); periods = 'x_period'
    X <- X %>% group_by_(.dots = id) %>% mutate(x_period = (1:.n)) %>% 
      as.data.frame()
  }
  if (is.null(intervention)) intervention <- max(X[X[[id]] == id.interest,periods])
  
  # Optimal bandwidth estimation
  if (is.null(h)){
    .temp <- X[X[[periods]] <= intervention,]
    h <- optimalBandwidth(.temp, k=k, id=id, id.interest=id.interest, 
                          kernel=kernel, method=method, norm.wgt=norm.wgt, 
                          controls=controls, out = F)$opt.h
  }
  
  # Kernels
  if (kernel == 'Gaussian'){
    .K <- function(u,h){ (1/sqrt(2*pi))*exp(-(((u)/h)^2)/2) }
  } else if (kernel == 'Epanechnikov'){
    .K <- function(u,h){ ifelse(abs(u/h)<=1, 3/4*(1-abs(u/h)^2), 0) }
  } else if (kernel == 'Tri-cube'){
    .K <- function(u,h){ ifelse(abs(u/h)<=1, (1-abs(u/h)^3)^3, 0) }
  } else if (kernel == 'Triangular'){
    .K <- function(u,h){ max(1 - (u/h), 0) }
  }
  
  # Calculate weights
  .wgt <- c()
  .units <- unique(X[[id]])[-which(unique(X[[id]]) == id.interest)]
  
  for (.u in .units){
    .IN <- X[X[[id]] == id.interest & X[[periods]] <= intervention,]
    .CO <- X[X[[id]] == .u & X[[periods]] <= intervention,]

    if (!is.null(dim(.IN[,controls]))){
      .vec <- apply(.IN[,controls] - .CO[,controls], 1, 
                    function(x) norm(x ,type = '2'))
    } else {
      .vec <- .IN[,controls] - .CO[,controls]
    }
    
    if (method == 'Aggregate') .w <- .K(.vec,h) %>% mean()
    if (method == 'Norm') .w <- .K(.vec,h) %>% norm(type='2')
    .wgt <- c(.wgt, .w)
  }
  
  # Name weights
  names(.wgt) <- .units
  
  # Normalize weights
  if (norm.wgt) .wgt <- .wgt / sum(.wgt)
  
  # Separate predictors and outcome variable
  .X <- X[X[[id]] != id.interest, c(id,periods,k)] %>% 
    pivot_wider(names_from = id, values_from = k) %>% 
    select(-all_of(c(periods))) %>% 
    as.matrix()
  
  .Y <- X[X[[id]] == id.interest, k]
  
  # Perfromance metrics
  .preds <- as.vector(.X %*% .wgt); names(.preds) <- unique(X[[periods]])
  .idx <- which(X[[id]] == id.interest & X[[periods]] > intervention)
  if (intervention < max(X[X[[id]] == id.interest,periods])){
    .mse <- mean((.preds[-.idx] - .IN[[k]][-.idx])^2)
  } else {
    .mse <- mean((.preds - .IN[[k]])^2)
  }
  if (intervention < max(X[X[[id]] == id.interest,periods]) &
      !is.null(true_cf)){
    .mspe <- mean((.preds[.idx] - true_cf[.idx])^2)
  } else {
    .mspe <- NULL
  }
  
  # Return results
  .l <- list(spec = paste('Kernel:',kernel,'| Bandwidth:',h,'| Method:',method),
             predictions = .preds, wgt = .wgt, mse = .mse)
  if(!is.null(.mspe)) .l$mspe <- .mspe
  return(.l)
}


##' @description
##' This function uses kernels to estimate the weights of different control
##' groups in a non-parametric way. It uses multivariate kernels or the product
##' of different univariate kernels. 
##' @param X: The data frame used for fitting the weights.
##' @param k: The variable of interest as a string.
##' @param h: The bandwidth used for the kernel. If set to `NULL`, which is the
##' default, then the function will choose the optimal bandwidth based on 
##' cross validation.
##' @param kernel: The kernel used. The default is `kernel = "Gaussian"`, but
##' one can also choose to instead use `Epanechnikov`, `Triangular`, or 
##' `Tri-cube`.
##' @param norm.wgt: Boolean value indicating whether or not to normalize the
##' estiamted weights, such that the sum of all weights is equal to 1. The 
##' default is set to `norm.wgt = FALSE`.
##' @return The function returns a list. The first element of that list is some
##' specifications about the method. The second element are the estimated 
##' weights.

multiKernel <- function(X, k, id, id.interest, periods = NULL, 
                        intervention = NULL, controls = NULL, h = NULL, 
                        kernel = 'Gaussian', method = 'Multivariate', 
                        norm.wgt = F, true_cf = NULL, standardize = F
                        # X, k, h = NULL, kernel = 'Gaussian', method = 'Product',
                        # norm.wgt = F
                        ){
  
  # Controls
  if (is.null(controls)){ controls <- c('x_ctrl'); X <- X %>% mutate(x_ctrl=X[[k]]) }
  
  # Standardize 
  if (standardize){
    X <- X %>% mutate_at(controls, ~(scale(.) %>% as.vector))
  }
  
  # Periods and intervention
  if (is.null(periods)){
    .n <- nrow(X[X[[id]] == id.interest,]); periods = 'x_period'
    X <- X %>% group_by_(.dots = id) %>% mutate(x_period = (1:.n)) %>% 
      as.data.frame()
  }
  if (is.null(intervention)) intervention <- max(X[X[[id]] == id.interest,periods])
  
  if (method == 'Product'){
    
    # Optimal bandwidth estimation
    if (is.null(h)){
      .temp <- X[X[[periods]] <= intervention,]
      h <- optimalBandwidth(.temp, k=k, id=id, id.interest=id.interest, 
                            kernel=kernel, method=method, norm.wgt=norm.wgt, 
                            controls=controls, out = F)$opt.h
    }
    
    # Kernels
    if (kernel == 'Gaussian'){
      .K <- function(u,h){ (1/sqrt(2*pi))*exp(-(((u)/h)^2)/2) }
    } else if (kernel == 'Epanechnikov'){
      .K <- function(u,h){ ifelse(abs(u/h)<=1, 3/4*(1-abs(u/h)^2), 0) }
    } else if (kernel == 'Tri-cube'){
      .K <- function(u,h){ ifelse(abs(u/h)<=1, (1-abs(u/h)^3)^3, 0) }
    } else if (kernel == 'Triangular'){
      .K <- function(u,h){ max(1 - (u/h), 0) }
    }
    
  } else if (method == 'Multivariate') {
    .x <- X[X[[id]] == id.interest & X[[periods]] <= intervention, controls]
    .K <- ks::kde(.x)
  }
  
  # Calculate weights
  .wgt <- c()
  .units <- unique(X[[id]])[-which(unique(X[[id]]) == id.interest)]
  
  for (.u in .units){
    .IN <- X[X[[id]] == id.interest & X[[periods]] <= intervention,]
    .CO <- X[X[[id]] == .u & X[[periods]] <= intervention,]
    
    if (!is.null(dim(.IN[,controls]))){
      .vec <- apply(.IN[,controls] - .CO[,controls], 1, 
                    function(x) norm(x ,type = '2'))
    } else {
      .vec <- .IN[,controls] - .CO[,controls]
    }
    
    if (method == 'Product') .w <- .K(.vec,h) %>% prod()
    if (method == 'Multivariate') .w <- predict(.K, x = .CO[,c(4:6)]) %>% mean()
    
    .wgt <- c(.wgt, .w)
  }
  
  # Name weights
  names(.wgt) <- .units
  
  # Normalize weights
  if (norm.wgt) .wgt <- .wgt / sum(.wgt)
  
  # Separate predictors and outcome variable
  .X <- X[X[[id]] != id.interest, c(id,periods,k)] %>% 
    pivot_wider(names_from = id, values_from = k) %>% 
    select(-all_of(c(periods))) %>% 
    as.matrix()
  
  .Y <- X[X[[id]] == id.interest, k]
  
  # Perfromance metrics
  .preds <- as.vector(.X %*% .wgt); names(.preds) <- unique(X[[periods]])
  .idx <- which(X[[id]] == id.interest & X[[periods]] > intervention)
  if (intervention < max(X[X[[id]] == id.interest,periods])){
    .mse <- mean((.preds[-.idx] - .IN[[k]][-.idx])^2)
  } else {
    .mse <- mean((.preds - .IN[[k]])^2)
  }
  if (intervention < max(X[X[[id]] == id.interest,periods]) &
      !is.null(true_cf)){
    .mspe <- mean((.preds[.idx] - true_cf[.idx])^2)
  } else {
    .mspe <- NULL
  }
  
  # Return results
  .l <- list(predictions = .preds, wgt = .wgt, mse = .mse)
  if(!is.null(.mspe)) .l$mspe <- .mspe
  return(.l)
}


##' @description
##' This function allows us to easily compute the cross entropy of two given
##' vectors.
##' @param x: A vector of numerical values.
##' @param y: A vector of numerical values.
##' @param h: The bandwidth for bins in the discrete case.
##' @param discrete: A boolean value indicating whether or not to use the 
##' discrete version of the cross-entropy. By default, `discrete = T`.
##' @param standardize: A boolean value indicating whether or not to standardize
##' the input vectors such that mean is zero and standard deviation is one. By
##' default, `standardize = F`.
##' @return The function returns a single numerical value.

klDiv <- function(X, k, id, id.interest, periods = NULL, 
                  intervention = NULL, controls = NULL, h = 0.5, 
                  discrete = T, KL = T, method = 'Norm', norm.wgt = F, 
                  true_cf = NULL, standardize = F){
  
    # Controls
  if (is.null(controls)){ controls <- c('x_ctrl'); X <- X %>% mutate(x_ctrl=X[[k]]) }
  
  # Standardize 
  if (standardize){
    X <- X %>% mutate_at(controls, ~(scale(.) %>% as.vector))
  }
  
  # Periods and intervention
  if (is.null(periods)){
    .n <- nrow(X[X[[id]] == id.interest,]); periods = 'x_period'
    X <- X %>% group_by_(.dots = id) %>% mutate(x_period = (1:.n)) %>% 
      as.data.frame()
  }
  if (is.null(intervention)) intervention <- max(X[X[[id]] == id.interest,periods])
  
  
  # Calculate weights
  .units <- unique(X[[id]])[-which(unique(X[[id]]) == id.interest)]
  .mat <- matrix(NA, length(controls), length(.units))
  rownames(.mat) <- controls; colnames(.mat) <- .units
  
  for (.u in .units){
    .IN <- X[X[[id]] == id.interest,]
    .CO <- X[X[[id]] == .u,]
    
    for (.c in controls){
      .in <- .IN[[.c]]; .co <- .CO[[.c]]
      
      # Cross-entropy
      if (discrete){
        .upper <- ceiling(max(.in, .co))
        .lower <- floor(min(.in, .co))
        .x <- cut(.in, breaks = c(seq(.lower,.upper,h))) %>% table() / length(.in)
        .y <- cut(.co, breaks = c(seq(.lower,.upper,h))) %>% table() / length(.co)
        if (KL == F){
          .ce <- - sum(.x * ifelse(.y == 0, 0, log(.y))) 
        } else {
          .ce <- abs(sum(.x * ifelse(.y == 0 | .x == 0, 0, log(.x / .y))))
        }
      } else {
        .dx <- density(.in); .dy <- density(.co)
        f <- function(x){ 
          .x <- approx(.dx$x, .dx$y, xout = x)$y
          .y <- approx(.dy$x, .dy$y, xout = x)$y
          .x <- ifelse(is.na(.x),0,.x)
          .y <- ifelse(is.na(.y) | .y == 0, 0, log(.y))
          return(.x - .y)
        }
        g <- function(x){
          .x <- approx(.dx$x, .dx$y, xout = x)$y
          .y <- approx(.dy$x, .dy$y, xout = x)$y
          .x <- ifelse(is.na(.x), 0, .x)
          .z <- ifelse(is.na(.y) | .y == 0 | .x == 0, 0, log(.x / .y))
          return(.x * .z)
        }
        if (KL == F){
          .ce <- cubintegrate(f, -Inf, Inf)$integral
        } else {
          .ce <- abs(cubintegrate(g, -Inf, Inf)$integral)
        }
      }
      
      # Save results
      .mat[.c, toString(.u)] <- .ce
    }
  }
  
  # Weights
  if (method == 'Aggregate') .wgt <- apply(1-(sigmoid(.mat)-0.5)*2, 2, mean)
  if (method == 'Norm') .wgt <- apply(1-(sigmoid(.mat)-0.5)*2, 2, 
                                      function(x) norm(x, type = '2'))
  
  # Normalize weights
  if (norm.wgt) .wgt <- .wgt / sum(.wgt)
  
  # Separate predictors and outcome variable
  .X <- X[X[[id]] != id.interest, c(id,periods,k)] %>% 
    pivot_wider(names_from = id, values_from = k) %>% 
    select(-all_of(c(periods))) %>% 
    as.matrix()
  
  .Y <- X[X[[id]] == id.interest, k]
  
  # Perfromance metrics
  .preds <- as.vector(.X %*% .wgt); names(.preds) <- unique(X[[periods]])
  .idx <- which(.IN[[periods]] > intervention)
  if (intervention < max(X[X[[id]] == id.interest,periods])){
    .mse <- mean((.preds[-.idx] - .IN[[k]][-.idx])^2)
  } else {
    .mse <- mean((.preds - .IN[[k]])^2)
  }
  if (intervention < max(X[X[[id]] == id.interest,periods]) &
      !is.null(true_cf)){
    .mspe <- mean((.preds[.idx] - true_cf[.idx])^2)
  } else {
    .mspe <- NULL
  }
  
  # Return results
  .l <- list(method = method, cross_entropy = .mat, predictions = .preds, 
             wgt = .wgt, mse = .mse)
  if(!is.null(.mspe)) .l$mspe <- .mspe
  return(.l)
}


##' @description
##' This function is used to make predictions when weights are used, such as in
##' the kernel methods.
##' @param X: The covariates used for prediction.
##' @param w: The weights used for the prediction, i.e. the ones coming from 
##' the `uniKernel()` function.
##' @param intercept: Indicates the intercept value. If no value is set, the
##' function assumes no intercept. The default is set to `intercept = NULL`.
##' @param k: The variable of interest, i.e. the one that should be predicted.
##' If this is not provided, it is assumed that it is not in the data frame 
##' `X`.
##' @return The function returns a vector of the predictions. Additionally, if
##' `!is.null(k)`, it will compute the mean squared error of the predictions.
##' In that case, the function returns a list with the weights and the MSE.

predict.weighted <- function(X, w, intercept = NULL, k = NULL){
  .X <- X
  if (! is.null(k)) .X <- .X %>% select(-any_of(c(k)))
  .i <- ifelse(! is.null(intercept), intercept, 0)
  .X <- as.matrix(.X)
  .p = list(preds = as.vector(.X %*% w) + .i)
  if (! is.null(k)){ .p$mse <- mean((X[[k]]-.p$preds)^2,na.rm=T) }
  return(.p)
}


##' @description
##' This function estimates the optimal bandwidth using cross validation and a
##' range of possible bandwidths. 
##' @param x: The covariates for which we want to estimate weights.
##' @param y: The variable of interest.
##' @param kernel: The kernel to be used, same as in `uniKernel()`.
##' @param method: The metheod to be used, same as in `uniKernel()`.
##' @param intercept: Whether or not to use an intercept, boolean variable.
##' @param norm.wgt: Whether or not to normalize the weights, boolean variable.
##' @return The function returns a list containing the optimal bandwidth and the
##' different MSEs for the different bandwidths tested.

optimalBandwidth <- function(x, k, id, id.interest, controls, kernel, method, 
                             norm.wgt = NULL, out = T, plot = F){
  
  # Possible bandwidths
  .poss <- c(seq(0.01,1.99,by=0.01),seq(2,5,by=0.1),6:50)
  # Progress bar
  cat('\nSearching for best bandwidth...\n')
  pb <- progress_bar$new(total = length(.poss))
  # Value storing
  .MSES <- c()
  # Length of timeseries
  .n <- nrow(x[x[[id]] == id.interest,])
  # Indexing 
  x <- x %>% group_by_(.dots = id) %>% mutate(idx = 1:.n) %>% ungroup()
  # Plotting
  if (plot) par(mfrow = c(4,2))
  
  for (bw in .poss){
    
    .mses <- c()
    
    # Cross validation
    if (out == T){
      
      .sampledIDX <- sample(1:(.n))
      
      for (.i in 1:5){
        .IDX <- .sampledIDX[((.i-1) * floor(.n/5) + 1) : (.i * floor(.n/5))]
        .x <- x[-which(x$idx %in% .IDX),]
        # Fitting
        if (method %in% c('Aggregate','Norm')){
          .z <- uniKernel(.x, k=k, h=bw, kernel=kernel, method=method, 
                          norm.wgt=norm.wgt, controls=controls, id=id,
                          id.interest=id.interest)
        } else if (method == 'Product'){
          
          ...
          
        } else if (method == 'Multivariate'){
          
          ##' @note Actually not needed.
          ...
          
        }
        if (any(is.na(.z$wgt))){
          .mse <- NA
        } else {
          xx <- x[which(x$idx %in% .IDX),]
          yy <- xx[xx[[id]] == id.interest,][[k]]
          XX <- xx[xx[[id]] != id.interest, c(id,periods,k)] %>% 
            pivot_wider(names_from = id, values_from = k) %>% 
            select(-all_of(c(periods))) %>% 
            as.matrix()
          
          .mse <- mean((XX %*% .z$wgt - yy)^2)
        }
        .mses <- c(.mses, .mse)
      }
      
    # Based on whole sample
    } else {
      
      if (method %in% c('Aggregate','Norm')){
        .z <- uniKernel(x, k=k, h=bw, kernel=kernel, method=method, 
                        norm.wgt=norm.wgt, controls=controls, id=id,
                        id.interest=id.interest)
        
        if (is.na(.z$mse)){
          yy <- x[x[[id]] == id.interest, k]
          .mse <- mean(as.vector((.z$predictions - yy))$y^2)
        } else {
          .mse <- .z$mse
        }
        .mses <- c(.mses, .mse)
      }
      
    }
    
    if (plot){
      plot(2000:(1999+n), data %>% filter(id == 1) %>% pull(y), type = 'l', col = .c[2], 
           ylim = c(-150,150), xlab = 'Year', ylab = 'GDP')
      lines(2000:(1999+n), data %>% filter(id == 2) %>% pull(y), type = 'l', col=.c[4])
      lines(2000:(1999+n), data %>% filter(id == 3) %>% pull(y), type = 'l', col=.c[8])
      lines(2000:(1999+n), data %>% filter(id == 4) %>% pull(y), type = 'l', col=.c[10])
      abline(v = p_int, lty = 2)
      lines(2000:(1999+nrow(x[x[[id]] == id.interest,])), 
            .z$predictions, type = 'l', col = .c[12], lwd = 3)
      legend(x = 'topright', legend = c('Treatment','Control 1','Control 2','Control 3','Synthetic'), 
             col = .c[c(2,4,8,10,12)], 
             lty= rep(1, 4))
    }
    
    .MSES <- c(.MSES, mean(.mses, na.rm=T))
    pb$tick()
  }
  names(.MSES) <- .poss
  opt.h <- .poss[which.min(.MSES)]
  return(list(MSEs = .MSES, opt.h = opt.h))
}


# Simulation Study ----

##' @description
##' This function is the data generating process (DGP). It creates a time series
##' with certain desired attributes that can be chosen in the function calling.
##' @param n: The length of the series.
##' @param form: Specifies the form of the function as a string, can be either 
##' `linear`, `power`, `quadratic`, `cubic`, or `oscillating`. The default is 
##' `form = 'quadratic'`.
##' @param alpha: This can be used to specify the alpha of the AR process used
##' to add noise. If this is not explicitly chosen, a random alpha will be 
##' selected.
##' @param beta: This can be used to specify the beta of the MA process used to
##' add noise. If this is not explicitly chosen, the noise simply comes from an
##' AR process.
##' @param gamma: This parameter squeezes the noise added to the time series. 
##' If not chosen, a random factor will be selected.
##' @param plot: If true, the generated time series will be plotted when the 
##' function is called. The default is that `plot = TRUE`.
##' @param max.diff: The maximal length of the range. The default is 
##' `max.diff = 10`. See the `squeeze()` function for a better overview.
##' @return The function returns a list with the generated time series, as well
##' as a vector of the parameters used for the functional form.

dgp <- function(n, form = 'quadratic', alpha = NULL, beta = NULL, gamma = NULL,
                plot = T, max.diff = 10) {
  # Set parameters for ARMA process
  .alpha <- ifelse(is.null(alpha), runif(1, 0, 0.999), alpha)
  .beta <- ifelse(is.null(beta), 0, beta)
  # Generate data
  .noise <- getTS(.alpha, .beta, n)
  
  .a <- round(runif(1, -5, 5),2)
  .b <- round(rnorm(1,0,2),2)
  .c <- round(runif(1,-2,2),1)
  .d <- round(rnorm(1,0,2),2)
  .e <- round(rnorm(1,0,1),1)
  .f <- round(rnorm(1,0,2),2)
  .g <- round(runif(1,1,5),0)
  .h <- round(rnorm(1,1,0.2),2)
  
  if (is.null(gamma)){
    .gamma <- ifelse(round(rnorm(1,0.3,0.2),2) == 0, 0.05, round(rnorm(1,0.3,0.2),2))
  } else {
    .gamma = gamma
  }
  
  if (form == 'linear'){
    .fn <- .a + squeeze(.b * c(1:n), max.diff = max.diff) + .gamma * .noise
    .params <- c(.a,.b,.gamma)
  } else if (form == 'power'){
    .fn <- .a + squeeze(.b * c(1:n)^.c, max.diff = max.diff) + .gamma * .noise
    .params <- c(.a,.b,.c,.gamma)
  } else if (form == 'quadratic'){
    .fn <- .a + squeeze(.b * c(1:n)^2 + .d * c(1:n), 
                        max.diff = max.diff) + .gamma * .noise
    .params <- c(.a,.b,.d,.gamma)
  } else if (form == 'cubic'){
    .fn <- .a + squeeze(.b * c(1:n)^3 + (n/2) * .d * c(1:n)^2 + .f * c(1:n), 
                        max.diff = max.diff) + .gamma * .noise
    .params <- c(.a,.b,(n/2)*.d,.f,.gamma)
  } else if (form == 'oscillating'){
    .fn <- .a + squeeze(.b * sin(1:n / .g) + .d * cos(1:n / .g), 
                        max.diff = max.diff) + .gamma * .noise
    .params <- c(.a,.b,.d,.g,.h,.gamma)
  } else {
    print('ERROR: No valid form was selected.')
  }
  if (plot == T){ plot(.fn, ylab = 'Function', main = 'Generated Time Series') }
  .r <- list(); .r$ts <- .fn; .r$params <- .params
  return(.r)
}


##' @description
##' The function generates a data set that can be used for the simulation study
##' to test different estimators.
##' @param n: The length of each time series.
##' @param forms: Specifies the forms of the functions as a string vector, can 
##' contain either `linear`, `power`, `quadratic`, `cubic`, or `oscillating`. 
##' The length of the vector automatically indicates the number of groups in the 
##' data. The first string in the vector is used for the treatment group.
##' @param alpha: This can be used to specify the alpha of the AR process used
##' to add noise. If this is not explicitly chosen, a random alpha will be 
##' selected.
##' @param beta: This can be used to specify the beta of the MA process used to
##' add noise. If this is not explicitly chosen, the noise simply comes from an
##' AR process.
##' @param gamma: This parameter squeezes the noise added to the time series. 
##' If not chosen, a random factor will be selected.
##' @param plot: If true, the generated data will be plotted when the 
##' function is called. The default is that `plot = TRUE`.
##' @param max.diff: The maximal length of the range. The default is 
##' `max.diff = 10`. See the `squeeze()` function for a better overview.
##' @return The function returns a list with the generated data, as well as a 
##' vector of the parameters used for the functional form.

DGP <- function(n, forms, alpha = NULL, beta = NULL, gamma = NULL,
                plot = T, max.diff = 10, diff.total = NULL){
  # Create matrix to store the data
  .d <- matrix(NA, n, length(forms))
  colnames(.d) <- c('T', paste0('X', 1:(length(forms)-1)))
  .p = list()
  # Simulate the data
  for (i in 1:length(forms)){
    .t <- dgp(n = n, form = forms[i], 
              alpha = alpha, beta = beta, gamma = gamma,
              plot = F, max.diff = max.diff)
    .d[,i] <- .t$ts; .p[[colnames(.d)[i]]] <- .t$params
  }
  .d <- as.data.frame(.d)
  if (!is.null(diff.total)){
    .d <- stack(.d) %>% 
      mutate(y = squeeze(y, max.diff = diff.total)) %>% 
      unstack()
  }
  # Plotting
  if (plot == T){
    plot(.d$`T`, ylab = 'y', xlab = 'Periods', main = '', 
         type = 'l', ylim = c(min(.d), max(.d)))
    for (s in 2:ncol(.d)){
      .s = .d %>% pull(s)
      lines(.s, type = 'l', col = rgb(0, 0, 0, 0.2))
    }
  }
  # Return object
  return(list(data = .d, params = .p))
}


##' @description
##' The DGP from Cerulli (2020). 
##' @param N: The amount of random variables generated per observation.
##' @param n: The amount of observations in the data set.
##' @param n_covariates: The number of covariates, by default, 
##' `n_covariates = 3`.
##' @param functional.form: The functional form when generating the data. Can
##' be either the original one from Cerulli (2020), or the `New` one. By 
##' default, however, `functional.form = Cerulli`.
##' @param intervention: The intervention period. If not provided, the 
##' intervention period will be after two thirds of the observations.
##' @param plot: Whether or not to plot the data when it's generated.
##' @return The function returns a list containing the data and the 
##' true counterfactual.

DGPcerulli <- function(N, n, n_covariates = 3, functional.form = 'Cerulli', 
                       intervention = NULL, plot = F){
  
  if (is.null(intervention)){
    p_int <- 1999 + floor(2/3 * n)
  } else {
    p_int <- intervention
  }
  
  # Generate data
  .data <- as.data.frame(matrix(rnorm(N * n_x * 3), N, n_x * 3))
  names(.data) <- paste0(paste0('x',rep(1:3,each=3)),'_',rep(1:3,3))
  
  if (functional.form == 'Cerulli'){
    
    .data <- .data %>% 
      mutate(x1 = rnorm(N) + x1_1^2 + abs(x1_2)^(0.5) + x1_3^3,
             x2 = rnorm(N) + exp(x2_1) + exp(1/x2_2) + (x2_3)^(2),
             x3 = rnorm(N) + (x3_1) / x3_2 + exp(x3_2) + (x3_3)^(-5)) %>% 
      mutate(y = rnorm(N) + x1 + x2 + x3) %>% 
      filter(abs(y) < 120) %>% 
      slice(1:n) %>% 
      mutate(year = 2000:(1999+n)) %>% 
      mutate(
        y0 = case_when(
          year <= p_int ~ y,
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + rnorm(n),
          T ~ NA_real_),
        id = 1,
        y1 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 2),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 15, 30),
          T ~ NA_real_),
        y2 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 3),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 10, 15),
          T ~ NA_real_),
        y3 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 4),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 5, 10),
          T ~ NA_real_),
        id1 = 2, id2 = 3, id3 = 4) %>%
      select(year, y, y0, y1, y2, y3, x1, x2, x3, everything())
    
  } else if (functional.form == 'New'){
    
    .data <- .data %>% 
      mutate(x1 = rnorm(N) + x1_1 + x1_2^2 + abs(x1_3)^(0.5),
             x2 = rnorm(N) + exp(x2_1) + exp(1/x2_2) + (x2_3)^(2),
             x3 = rnorm(N) + (x3_1) / x3_2 + exp(x3_2) + (x3_3)^(-3)) %>% 
      mutate(y = rnorm(N) + x1 + x2 + x3) %>% 
      filter(abs(y) < 120) %>% 
      slice(1:n) %>% 
      mutate(year = 2000:(1999+n)) %>% 
      mutate(
        y0 = case_when(
          year <= p_int ~ y,
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + rnorm(n),
          T ~ NA_real_),
        id = 1,
        y1 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 1),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 10, 5),
          T ~ NA_real_),
        y2 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 2),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 5, 3),
          T ~ NA_real_),
        y3 = case_when(
          year <= p_int ~ y + rnorm(n, 0, 3),
          year > p_int ~ -100 + x1 + Re(as.complex(x2)^0.7) + log(abs(x3)) + 
            rnorm(n, 2, 3),
          T ~ NA_real_),
        id1 = 2, id2 = 3, id3 = 4) %>%
      select(year, y, y0, y1, y2, y3, x1, x2, x3, everything())
    
  }
  
  # Set aside the DGP
  .DGP <- .data %>% 
    select(year, y, y0) %>% 
    rename(factual = y, counterfactual = y0)
  
  # Re-format the data
  for (i in 1:3){
    .d <- .data %>% 
      select(all_of(c(
        paste0('id',i), 'year',
        paste0('y',i), paste0('x',i,'_',1:3))
      ))
    names(.d) <- c('id','year','y','x1','x2','x3')
    assign(paste0('.data',i),.d)
  }
  
  .data <- .data %>% 
    select(id, year, y, x1, x2, x3) %>% 
    rbind(.data1, .data2, .data3)
  
  if (plot){
    plot(2000:(1999+n), .data %>% filter(id == 1) %>% pull(y), type = 'l', col = .c2[1], 
         ylim = c(-150,150), xlab = 'Year', ylab = 'GDP')
    lines(2000:(1999+n), .data %>% filter(id == 2) %>% pull(y), type = 'l', col=.c2[2])
    lines(2000:(1999+n), .data %>% filter(id == 3) %>% pull(y), type = 'l', col=.c2[3])
    lines(2000:(1999+n), .data %>% filter(id == 4) %>% pull(y), type = 'l', col=.c2[4])
    abline(v = p_int, lty = 2)
    legend(x = 'bottomleft', legend = c('Treatment','Control 1','Control 2','Control 3'), 
           col = .c2[1:4], 
           lty= rep(1, 3))
  }
  
  .l <- list(data = .data, DGP = .DGP); return(.l)
}


##' @description
##' This function tests the performance of different models on the data that
##' was generated before / on real-world data.
##' @param x: The data to be used to test the predictive performance of the
##' different models. Generally, the output of the `DGP()` function is used
##' for that.
##' @param t: This variable indicates the time of the policy intervention in the
##' case of synthetic controls. Basically this states up to which period we use
##' the data for training. Everything after will be used for calculating the 
##' test statistics of prediction accuracy.
##' @param models: A vector of strings indicating which models are to be fitted.
##' The options here are `LASSO`, `Abadie`, `GAM`, `baggedTrees`, `rf`,
##' `boostedTrees`, `uniKernel`, `multiKernel` and `LSTM`. If nothing is chosen,
##' the default is to use all these models.
##' @param plot: A boolean variable indicating whether or not to plot the data
##' with the fitted models. The default is that `plot = TRUE`.
##' @return The function returns a matrix of performance statistics for the
##' different models fitted to the data.

testModels <- function(x, t, models = NULL, plot = T){
  # Choose models to test on data
  if (is.null(models)){ models = c('LASSO', 'GAM', 'MARS', 'baggedTrees', 
                                   'boostedTrees', 'rf', 'uniKernel',
                                   'entropy') }
  # Split the data
  .train = x[1:t,]
  .test = x[(t+1):nrow(x),]
  # Rearrange data for methods with explicit weights
  x.expl <- stack(x) %>% group_by(id) %>% mutate(year = 1:n()) %>% 
    ungroup %>% as.data.frame()
  # Fit the models
  .m <- data.frame(Statistic = c('MSE','RMSE','MAE'))
  .p <- matrix(NA, nrow(x) - t, length(models)); colnames(.p) <- models
  .l <- list()
  
  for (m in models){
    
    if (m == 'LASSO'){
      
      .cv <- cv.glmnet(as.matrix(.train[,-which(names(.train)=='T')]),
                       .train$`T`, alpha=1, intercept=F)
      .lambda <- .cv$lambda.min
      .model <- glmnet(as.matrix(.train[,-which(names(.train)=='T')]),
                      .train$`T`, alpha=1, lambda=.lambda, intercept=F)
      .preds <- predict(.model, as.matrix(.test[,-which(names(.test)=='T')])) %>% 
        as.data.frame() %>% pull(s0)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'GAM'){
      
      .form <- wrapFormula(`T` ~ ., data = .train, wrapString = 's(*)')
      .model <- try({
        gam(.form, data = .train, select=T)
      }, silent = T)
      
      if (any(class(.model) == 'try-error')){
        .preds <- rep(NA, nrow(.test))
      } else {
        .preds <- predict(.model, .test)
      }
      
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
          
    } else if (m == 'MARS'){
      
      .model <- earth(`T` ~ ., data = .train)
      .preds <- predict(.model, .test)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'baggedTrees'){
      
      .model <- baggedTree(.train, 'T', B=200)
      .preds <- predict.baggedTree(.model, .test) 
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'boostedTrees'){
      
      .model <- boostedTree(.train, 'T', B=200, info = F)
      .preds <- predict.boostedTree(.model, .test)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'rf'){
      
      .model <- randomForest(`T` ~ ., data=.train)
      .preds <- predict(.model, .test)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
    
    } else if (m == 'SVM'){
      
      .model <- svm(`T` ~ ., data=.train)
      .preds <- predict(.model, .test)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'ProjPursuit'){
      
      .form <- as.formula(paste('`T` ~ ',paste(names(.train)[-1],collapse=' + ')))
      .mses <- c()
      
      # Cross-validating number of ridge regressions
      for (nterm in 1:5){
        cv <- crossv_kfold(.train)
        cv.models <- map(cv$train, ~ppr(.form , data = ., nterms = nterm))
        get_pred  <- function(model, test_data){
          data  <- as.data.frame(test_data);pred  <- add_predictions(data, model)
          return(pred)
        }
        cv.pred <- map2_df(cv.models, cv$test, get_pred, .id = "Run")
        cv.mse <- cv.pred %>% group_by(Run) %>%
          summarise(MSE = mean( (`T` - pred)^2)) %>% pull(MSE) %>%  mean()
        .mses <- c(.mses,cv.mse)
      }
      
      .model <- ppr(`T` ~ ., data=.train, nterms = which.min(.mses))
      .preds <- predict(.model, .test)
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'uniKernel'){
      
      .model <- uniKernel(x.expl, 'y', id = 'id', id.interest = 'T', 
                          periods = 'year', intervention = t, 
                          method = 'Aggregate', norm.wgt = T,
                          standardize = F)
      .preds <- .model$predictions[(t+1):length(.model$predictions)]
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'multiKernel'){
      
      .model <- multiKernel(x.expl, 'y', id = 'id', id.interest = 'T', 
                            periods = 'year', intervention = t, 
                            method = 'Product', norm.wgt = T)
      .preds <- .model$predictions[(t+1):length(.model$predictions)]
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    } else if (m == 'entropy'){
      
      .model <- klDiv(x.expl, 'y', id = 'id', id.interest = 'T',
                      periods = 'year', intervention = t, norm.wgt = T,
                      standardize = F)
      .preds <- .model$predictions[(t+1):length(.model$predictions)]
      .l[[m]] <- c(x$`T`[1:t], .preds) %>% unname()
      
    }
    
    # Stats
    .mse <- mean((x$`T`[(t+1):nrow(x)] - .preds)^2, na.rm = T)
    .rmse <- sqrt(.mse)
    .mae <- mean(abs(x$`T`[(t+1):nrow(x)] - .preds), na.rm = T)
    # Save estimates
    .p[,m] <- .preds
    .m[[m]] <- c(.mse,.rmse,.mae)
    
  }
  
  # Change to data frame
  .p <- as.data.frame(.p)
  # Plot
  if (plot == T){
    # Variable of interest
    plot(x$`T`, ylab = 'Function', xlab = 'Periods', main = 'Model Fits', 
         type = 'l', ylim = c(min(x), max(x)))
    # Control groups
    for (s in 2:ncol(x)){
      .s = x %>% pull(s)
      lines(.s, type = 'l', col = rgb(0, 0, 0, 0.2))
    }
    # Estimates
    # .c <- c('red','blue','green','purple','orange','pink','aquamarine4',
    #         'forestgreen','coral3')
    .i <- 1
    for (m in models){
      lines((t+1):nrow(x), .p[[m]], type = 'l', col = .c[.i])
      .i <- .i + 1
    }
    # Add legend
    legend(x = 'topright', legend = models, col = .c[1:length(models)], 
           lty= rep(1, length(models)))
  }
  # Return object
  .l$stats <- .m
  return(.l)
}


##' @description
##' This function is used to create the desired plots for simulation studies,
##' namely boxplots and line plots for the different performance metrics used.
##' @param x: The data frame with the performance metrics of the different 
##' methods on simulated data.
##' @param simulation: The type of simulation that is performed as a string. 
##' This will be used for the axis labeling of the plots. E.g., when we simulate
##' data with different time series lengths, we can put 
##' `simulation = "Series Length"`.
##' @param metrics: The performance metrics to be calculated as a string in a
##' vector. The possible options are `MSE`, `RMSE`, and `MAE`. The default is 
##' to use all three of them. 
##' @param file: If this is not specified, the plots will simply be showed, but
##' not saved anywhere. However, one can give the plots a name here (i.e., a 
##' string value) and then the plots will automatically be saved in the Figures
##' folder.
##' @param aspect.ratio: Gives the aspect ratio of the plot, by default, height
##' is 4 and width is 12.

simPlot <- function(x, simulation, metrics = c('MSE','RMSE','MAE'), file = NULL,
                    aspect.ratio = c(4,12), line.plot = T){
  
  if (length(metrics)<=3){
    par(mfrow = c(1,length(metrics)))
  } else if (length(metrics)==4){
    par(mfrow = c(2,2))
  }
  
  if (line.plot){ types <- c('boxplot','lines') } else { types <- c('boxplot') }
  
  for (type in types){
    
    if (!is.null(file)){
      pdf(paste0(here(),'/../Paper/Figures/Simulation/',file,'_',type,'.pdf'),
          height = aspect.ratio[1],
          width = aspect.ratio[2]) 
    }
    
    for (metric in metrics) {
      # Prep
      .Stats <- x %>% filter(Statistic == metric)
      for (i in 3:ncol(.Stats)){
        if (i == 3){ 
          .v <- .Stats[[i]]; .g <- rep(names(.Stats)[i], nrow(.Stats))
        } else { 
          .v <- c(.v, .Stats[[i]]); .g <- c(.g, rep(names(.Stats)[i], nrow(.Stats)))
        }
      }
      # Boxplot
      if (type == 'boxplot'){
        # Fix the order so that the coloring is the same in both plot types
        .data <- data.frame(value = .v, group = .g)
        .data$group <- factor(.data$group , levels=names(.Stats[-c(1:2)]))
        # Plot
        boxplot(value ~ group, data = .data,
                xlab = '', ylab = metric, col = .c2[1:(ncol(.Stats)-2)])
        # Line plots
      } else if (type == 'lines'){
        if (length(metrics) > 1){
          .xaxis <- unique(.Stats$Value); .order <- 1:length(.xaxis)
        } else {
          .xaxis <- sort(.Stats$Value); .order <- order(.Stats$Value)
        }
        plot(.xaxis, .Stats[[3]][.order], type = 'l', col = .c2[1],
             ylab = metric, xlab = simulation, 
             ylim = c(0, max(.Stats[-c(1:2)])))
        if (ncol(.Stats) > 3){
          .i <- 2
          for (c in 4:ncol(.Stats)){
            lines(.xaxis, .Stats[[c]][.order], type = 'l', col = .c2[.i])
            .i <- .i + 1
          }
        }
        # Add Legend
        legend(x = 'topright', legend = names(.Stats)[3:ncol(.Stats)], 
               col = .c2[1:(ncol(.Stats) - 2)], 
               lty= rep(1, (ncol(.Stats) - 2)),
               bg = rgb(1, 1, 1, 0.7))
      }
    }
    if (!is.null(file)){ dev.off() }
  }
}


# Other ----

##' @description
##' This function generates a time series.
##' @param alpha: The alpha parameter of the AR-part of the model.
##' @param beta: The beta parameter of the MA-part of the model, if needed.
##' @param n: The desired length of the time series.
##' @return The function returns the time series as a vector.

getTS <- function(alpha, beta, n) {
  if (beta == 0){
    .ts <- arima.sim(list(order = c(1,0,0), ar = alpha), n = n)
  } else {
    .ts <- arima.sim(list(order = c(1,0,1), ar = alpha, ma = beta), n = n)
  }
  return(.ts)
}


##' @description
##' This function squeezes a time series to a desired range.
##' @param ts: The time series to be squeezed.
##' @param max.diff: The maximal length of the range. The default is 
##' `max.diff = 10`.
##' @return The function returns the squeezed time series.

squeeze <- function(ts, max.diff = 10){
  if (max(ts) - min(ts) > max.diff){
    .ts <- ts - min(ts)
    .ts <- .ts * (max.diff / max(.ts))
    ts <- .ts + (ts[1] - .ts[1])
    
    # .factor <- abs(ifelse(abs(max(ts))>abs(min(ts)),max(ts),min(ts)))
    # ts <- ts * max.diff / .factor
  }
  return(ts)
}


##' @description
##' This function takes a data frame and stacks all the columns on top of each
##' other.
##' @param x: The data frame to stack.
##' @return A data frame with one column indicating the name of the previous 
##' column and one with the output variables.

stack <- function(x){
  for (i in 1:ncol(x)){
    if (i == 1){
      .d <- x[i] %>% mutate(id = names(x)[i]); names(.d) <- c('y','id')
    } else {
      .t <- x[i] %>% mutate(id = names(x)[i]); names(.t) <- c('y','id') 
      .d <- rbind(.d, .t)
    }
  }
  return(.d %>% select(id, y))
}


##' @description
##' This function takes a data frame and unstacks it into different columns.
##' @param x: The data frame to unstack.
##' @return A data frame with multiple columns.

unstack <- function(x){
  .ls <- list()
  for (id in unique(x[[1]])){
    .ls[[id]] <- x[x[[1]] == id, 2]
  }
  return(data.frame(.ls))
}


##' @description
##' The sigmoid function.
##' @param x: A scalar or numerical vector.
##' @return A scalar or numerical vector of the same length as the input.

sigmoid <- function(x){ y <- 1/(1+exp(-x)); return(y) }


##' @description
##' This function is there to completely mute all possible outputs in the 
##' console, e.g., warnings, messages, etc.

inSilence = function(...){
  .mc = match.call()[-1]
  .a = capture.output(
    tryCatch(suppressMessages(suppressWarnings(eval(as.list(.mc)[[1]]))), 
             error = function(e) ""))
}

