##' @title  Synthetic Control Groups -- Dissimilarity Matrix
##' @author Pascal Amiet
##' @date   13.11.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('here','Synth','stringr','readr','xtable') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')


# Performance -------------------------------------------------------------

# Set parameters
SEED <- 1337
set.seed(SEED)
reps <- 200
plot = F

# Matrix for results
results <- matrix(NA_real_,reps,17)
methods <- c('Synth','LASSO','Manhattan','Euclidean','Correlation')
colnames(results) <- c(
  paste0(rep(c('mse_','mspe_'),each=5),methods),
  paste0(paste0(rep(c('mse_','mspe_'),each=5),methods)[c(3:5,8:10)],'_norm'),
  'length'
)

functional.forms <- c('New','Cerulli')
progress <- progress_bar$new(total = reps)


# Simulation
for (rep in 1:reps){
  
  # Parameters
  N <- 200
  n <- ceiling(runif(1,10,100))
  n_x <- 3
  p_int <- 1999 + floor(2/3 * n)
  
  # Simulate data
  temp.obj <- DGPcerulli(N, n, functional.form = 'New')
  data <- temp.obj$data; DGP <- temp.obj$DGP
  
  # Compute synthetic control group
  dataprep.out <- dataprep(foo = data, predictors = c('x1', 'x2', 'x3'), 
                           predictors.op = "mean",
                           time.predictors.prior = 2000:2009, dependent = "y",
                           unit.variable = "id", time.variable = "year", 
                           treatment.identifier = 1, controls.identifier = c(2:4),
                           time.optimize.ssr = 2000:p_int, time.plot = 2000:(1999+n))
  
  inSilence({
    Synth.out <<- synth(dataprep.out)
  })
  
  # LASSO
  for (.i in unique(data$id)){
    .d <- data %>% filter(id == .i) %>% select(-id)
    names(.d) <- c('year',paste0(names(data)[3:6],'_',.i))
    .n <- paste0('.d',.i); assign(.n, .d)
  }
  
  data.lasso <- .d1 %>% 
    left_join(.d2, by = 'year') %>% 
    left_join(.d3, by = 'year') %>% 
    left_join(.d4, by = 'year')
  
  x.lasso <- data.lasso %>% filter(year <= p_int) %>% select(-y_1) %>% as.matrix()
  y.lasso <- data.lasso %>% filter(year <= p_int) %>% pull(y_1)
  cv.model <- cv.glmnet(x.lasso, y.lasso, alpha=1,intercept=F)
  best.lambda <- cv.model$lambda.min
  lasso.out <- glmnet(x.lasso, y.lasso, alpha=1, lambda=best.lambda, intercept=F)
  
  # Kernel synthetic control groups
  for (bool.norm in c(T, F)){
    for (dist in methods[-1]){
      inSilence({
        x <<- dissimilarityMatrix(X = data, k = 'y', id, id.interest = 1, 
                                  periods = 'year', intervention = p_int, 
                                  controls = c('x1','x2','x3'), 
                                  distance = dist, method = 'Aggregate', 
                                  norm.wgt = bool.norm, 
                                  true_cf = DGP$counterfactual,
                                  standardize = T)
      })
      .n <- ifelse(bool.norm, paste0(dist,'_norm.out'), paste0(dist,'.out'))
      assign(.n, x)
    }
  }
  
  # Save results
  preds_synth <- as.vector(dataprep.out$Y0plot %*% synth.out$solution.w)
  preds_lasso <- predict(lasso.out, data.lasso %>% select(-y_1) %>% as.matrix()) %>% 
    as.data.frame() %>% pull(s0)
  idx <- (2000:p_int)-1999
  mse <- function(e){ mean((e)^2) }
  
  results[rep,'mse_Synth'] <- mse(preds_synth[idx] - 
                                    DGP$counterfactual[idx])
  results[rep,'mse_LASSO'] <- mse(preds_lasso[idx] - 
                                    DGP$counterfactual[idx])
  results[rep,'mse_Euclidean'] <- mse(Euclidean.out$predictions[idx] - 
                                       DGP$counterfactual[idx])
  results[rep,'mse_Manhattan'] <- mse(Manhattan.out$predictions[idx] - 
                                           DGP$counterfactual[idx])
  results[rep,'mse_Correlation'] <- mse(Correlation.out$predictions[idx] - 
                                       DGP$counterfactual[idx])
  
  results[rep,'mspe_Synth'] <- mse(preds_synth[-idx] - 
                                     DGP$counterfactual[-idx])
  results[rep,'mspe_LASSO'] <- mse(preds_lasso[-idx] - 
                                     DGP$counterfactual[-idx])
  results[rep,'mspe_Euclidean'] <- mse(Euclidean.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  results[rep,'mspe_Manhattan'] <- mse(Manhattan.out$predictions[-idx] - 
                                            DGP$counterfactual[-idx])
  results[rep,'mspe_Correlation'] <- mse(Correlation.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])

  results[rep,'mse_Euclidean_norm'] <- mse(Euclidean_norm.out$predictions[idx] - 
                                            DGP$counterfactual[idx])
  results[rep,'mse_Manhattan_norm'] <- mse(Manhattan_norm.out$predictions[idx] - 
                                                DGP$counterfactual[idx])
  results[rep,'mse_Correlation_norm'] <- mse(Correlation_norm.out$predictions[idx] - 
                                            DGP$counterfactual[idx])
  
  results[rep,'mspe_Euclidean_norm'] <- mse(Euclidean_norm.out$predictions[-idx] - 
                                             DGP$counterfactual[-idx])
  results[rep,'mspe_Manhattan_norm'] <- mse(Manhattan_norm.out$predictions[-idx] - 
                                                 DGP$counterfactual[-idx])
  results[rep,'mspe_Correlation_norm'] <- mse(Correlation_norm.out$predictions[-idx] - 
                                             DGP$counterfactual[-idx])
  results[rep,'length'] <- n
  
  # Track progress
  progress$tick()
  
}

# Save results
saveRDS(results, paste0(here(),'/data/simulations/dissimilarity_standard[',SEED,'].rds'))


# Outputs -----------------------------------------------------------------

# Load data
data  <- readRDS(paste0(here(),'/data/simulations/dissimilarity_standard[',SEED,'].rds'))

# Re-arrange data
d.plots <- data %>% 
  as.data.frame() %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         L1 = mspe_Manhattan, L2 = mspe_Euclidean, Cor = mspe_Correlation,
         L1_norm = mspe_Manhattan_norm, L2_norm = mspe_Euclidean_norm,
         Cor_norm = mspe_Correlation_norm) %>% 
  select(Statistic, Value, everything())

# Plots
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'), 
        file = paste0('dissimilarity_comp'), aspect.ratio = c(6,6))

# Table
tbl <- matrix(NA, 8, 4); colnames(tbl) <- c('Method','Normalized','Mean RMSPE','Variance')
tbl[,1] <- c(names(d.plots)[3:4], rep(names(d.plots)[5:7], 2))
tbl[,2] <- c('Yes', 'No', rep(c('No','Yes'), each = 3))
tbl[1:5,3] <- d.plots %>% select(3:7) %>% apply(2, mean) %>% round(3)
tbl[1:5,4] <- d.plots %>% select(3:7) %>% apply(2, var) %>% round(3)
tbl[6:8,3] <- d.plots %>% select(8:10) %>% apply(2, mean) %>% round(3)
tbl[6:8,4] <- d.plots %>% select(8:10) %>% apply(2, var) %>% round(3)

saveRDS(tbl, paste0(here(),'/data/simulations/tbl_dissimilarity_standard[',SEED,'].rds'))

print(xtable(tbl,type='latex',caption='Statistics of Dissimilarity Comparison',
             label='tab:dissimilarity_comp',align=c('l','l',rep('r',3)),
             table.placement='b',
             hline.after = c(3)),
      include.rownames=FALSE,
      file=paste0(here(),'/../Paper/Tables/dissimilarity_comp.tex'))


pdf(paste0(here(),'/../Paper/Figures/Simulation/dissimilarity_series_length.pdf'),
    height = 5, width = 8) 

plot(rep(d.plots$Value, 5), stack(d.plots %>% select(c(3:7)))$y,
     xlab = 'Time Series Length', ylab = 'RMSPE', cex = 0.5, pch = 16,
     col = rgb(0, 0, 0, 0.2))
.m <- names(d.plots)[3:7]
for (i in 1:length(.m)){
  .val <- paste0('`',.m[i],'`')
  .frm <- as.formula(paste(.val,'~ Value'))
  .reg <- lm(.frm, data = d.plots)
  lines(11:100, .reg$coefficients[1] + .reg$coefficients[2] * c(11:100),
        col = .c2[i])
}
legend(x = 'topright', legend = .m, 
       col = .c2[1:length(.m)], 
       lty= rep(1, length(.m)))

dev.off()
