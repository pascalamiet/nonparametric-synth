##' @title  Synthetic Control Groups -- Kernel Methods
##' @author Pascal Amiet
##' @date   13.11.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('here','Synth','stringr','readr','xtable','gplm') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')


# Kernel Plots ------------------------------------------------------------

# Kernels
.GK <- function(u,h){ (1/sqrt(2*pi))*exp(-(((u)/h)^2)/2) }
.EK <- function(u,h){ ifelse(abs(u/h)<=1, 3/4*(1-abs(u/h)^2), 0) }
.TK <- function(u,h){ ifelse(abs(u/h)<=1, (1-abs(u/h)^3)^3, 0) }

# Save plot
pdf(paste0(here(),'/../Paper/Figures/kernel_overview.pdf'),
    height = 5,
    width = 10) 

plot(seq(-3,3,0.01), .TK(seq(-3,3,0.01),1), type = 'l', main = '', 
     ylim = c(0,1), xlab = '', ylab = '', col = .c2[1], lwd = 3)

lines(seq(-3,3,0.01), .EK(seq(-3,3,0.01),1), col = .c2[2], lwd = 3)
lines(seq(-3,3,0.01), .GK(seq(-3,3,0.01),1), col = .c2[3], lwd = 3)

legend(x = 'topright', legend = c('Gaussian','Epanechnikov','Tri-Cube'), 
       col = .c2[c(3,2,1)], lty= rep(1, 3), lwd = 3)

dev.off()


# Multivariate Kernels ----------------------------------------------------

# ks::amise.mixt()
# ks::Hamise.mixt()
# ks::kde()
# 
# 
# # Simulated data from a bivariate normal
# n <- 200
# set.seed(35233)
# x <- mvtnorm::rmvnorm(n = n, mean = c(0, 0),
#                       sigma = rbind(c(1.5, 0.25), c(0.25, 0.5)))
# 
# null <- cbind(rnorm(100),rnorm(100),rnorm(100),rnorm(100),rnorm(100),rnorm(100))
# 
# # Compute kde for a diagonal bandwidth matrix (trivially positive definite)
# H <- diag(c(1.25, 0.75))
# kde <- ks::kde(x = x, H = H)
# 
# # The eval.points slot contains the grids on x and y
# str(kde$eval.points)
# ## List of 2
# ##  $ : num [1:151] -8.58 -8.47 -8.37 -8.26 -8.15 ...
# ##  $ : num [1:151] -5.1 -5.03 -4.96 -4.89 -4.82 ...
# 
# # The grids in kde$eval.points are crossed in order to compute a grid matrix
# # where to evaluate the estimate
# dim(kde$estimate)
# ## [1] 151 151
# 
# # Manual plotting using the kde object structure
# image(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate,
#       col = viridis::viridis(20), xlab = "x", ylab = "y")
# points(kde$x) # The data is returned in $x
# 
# approxfun(kde$x, kde$estimate)
# 
# H <- Hscv(..., amise = TRUE)
# 
# 
# .mkd <- gplm::kde(.IN[,c(4:6)])
# .mkd <- ks::kde(.IN[,c(4:6)])

x <- multiKernel(X = data, k = 'y', id = 'id', id.interest = 1,
                periods = 'year', intervention = p_int, 
                controls = c('x1','x2','x3'), kernel = 'Gaussian',
                method = 'Product', true_cf = DGP$counterfactual,
                norm.wgt = T)


# Performance -------------------------------------------------------------

# Set parameters
SEED <- 1337
set.seed(SEED)
reps <- 100
plot = F

# Matrix for results
results <- matrix(NA_real_,reps,43)
methods <- c('Synth','LASSO','Gaussian','Epanechnikov','Tri-cube')
colnames(results) <- c(
  c(paste0(rep(c('mse_','mspe_'),each=5),methods),
    paste0('h_',methods[c(3:5)])),
  paste0(c(paste0(rep(c('mse_','mspe_'),each=5),methods),
           paste0('h_',methods[c(3:5)]))[c(3:5,8:13)],'_norm'),
  'mse_multi','mse_prod','mspe_multi','mspe_prod',
  'mse_multi_norm','mse_prod_norm','mspe_multi_norm','mspe_prod_norm',
  'rho_Synth','rho_LASSO',
  'rho_Gaussian','rho_Epanechnikov','rho_Tri-cube','rho_multi','rho_prod',
  'rho_Gaussian_norm','rho_Epanechnikov_norm','rho_Tri-cube_norm',
  'rho_multi_norm','rho_prod_norm',
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
  
  # Rho
  .rho <- treatEff(DGP$factual, DGP$counterfactual, p_int - 1999)
  
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
    for (kernel in methods[c(3:5)]){
      inSilence({
        x <<- uniKernel(X = data, k = 'y', id = 'id', id.interest = 1,
                        periods = 'year', intervention = p_int, 
                        controls = c('x1','x2','x3'), kernel = kernel,
                        method = 'Aggregate', true_cf = DGP$counterfactual,
                        norm.wgt = bool.norm, standardize = T)
      })
      .n <- ifelse(bool.norm, 
                   paste0(paste0(str_split(kernel,'-')[[1]],collapse='_'),'_norm.out'),
                   paste0(paste0(str_split(kernel,'-')[[1]],collapse='_'),'.out'))
      assign(.n, x)
    }
    
    for (method in c('Multivariate','Product')){
      inSilence({
        x <<- multiKernel(X = data, k = 'y', id = 'id', id.interest = 1,
                          periods = 'year', intervention = p_int, 
                          controls = c('x1','x2','x3'), kernel = 'Gaussian',
                          method = method, true_cf = DGP$counterfactual,
                          norm.wgt = bool.norm, standardize = T)
      })
      .n <- ifelse(bool.norm, paste0(method,'_norm.out'), paste0(method,'.out'))
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
  results[rep,'mse_Gaussian'] <- mse(Gaussian.out$predictions[idx] - 
                                       DGP$counterfactual[idx])
  results[rep,'mse_Epanechnikov'] <- mse(Epanechnikov.out$predictions[idx] - 
                                           DGP$counterfactual[idx])
  results[rep,'mse_Tri-cube'] <- mse(Tri_cube.out$predictions[idx] - 
                                       DGP$counterfactual[idx])
  
  results[rep,'mspe_Synth'] <- mse(preds_synth[-idx] - 
                                     DGP$counterfactual[-idx])
  results[rep,'mspe_LASSO'] <- mse(preds_lasso[-idx] - 
                                     DGP$counterfactual[-idx])
  results[rep,'mspe_Gaussian'] <- mse(Gaussian.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  results[rep,'mspe_Epanechnikov'] <- mse(Epanechnikov.out$predictions[-idx] - 
                                            DGP$counterfactual[-idx])
  results[rep,'mspe_Tri-cube'] <- mse(Tri_cube.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  
  results[rep,'h_Gaussian'] <- as.numeric(str_split(
    str_split(Gaussian.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  results[rep,'h_Epanechnikov'] <- as.numeric(str_split(
    str_split(Epanechnikov.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  results[rep,'h_Tri-cube'] <- as.numeric(str_split(
    str_split(Tri_cube.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  
  
  results[rep,'mse_Gaussian_norm'] <- mse(Gaussian_norm.out$predictions[idx] - 
                                            DGP$counterfactual[idx])
  results[rep,'mse_Epanechnikov_norm'] <- mse(Epanechnikov_norm.out$predictions[idx] - 
                                                DGP$counterfactual[idx])
  results[rep,'mse_Tri-cube_norm'] <- mse(Tri_cube_norm.out$predictions[idx] - 
                                            DGP$counterfactual[idx])
  
  results[rep,'mspe_Gaussian_norm'] <- mse(Gaussian_norm.out$predictions[-idx] - 
                                             DGP$counterfactual[-idx])
  results[rep,'mspe_Epanechnikov_norm'] <- mse(Epanechnikov_norm.out$predictions[-idx] - 
                                                 DGP$counterfactual[-idx])
  results[rep,'mspe_Tri-cube_norm'] <- mse(Tri_cube_norm.out$predictions[-idx] - 
                                             DGP$counterfactual[-idx])
  
  results[rep,'h_Gaussian_norm'] <- as.numeric(str_split(
    str_split(Gaussian_norm.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  results[rep,'h_Epanechnikov_norm'] <- as.numeric(str_split(
    str_split(Epanechnikov_norm.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  results[rep,'h_Tri-cube_norm'] <- as.numeric(str_split(
    str_split(Tri_cube_norm.out$spec,pattern='Bandwidth: ')[[1]][2],' | ')[[1]][1])
  
  results[rep,'mse_multi'] <- mse(Multivariate.out$predictions[idx] - 
                                         DGP$counterfactual[idx])
  results[rep,'mse_prod'] <- mse(Product.out$predictions[idx] - 
                                        DGP$counterfactual[idx])
  results[rep,'mse_multi_norm'] <- mse(Multivariate_norm.out$predictions[idx] - 
                                            DGP$counterfactual[idx])
  results[rep,'mse_prod_norm'] <- mse(Product_norm.out$predictions[idx] - 
                                                DGP$counterfactual[idx])
  results[rep,'mspe_multi'] <- mse(Multivariate.out$predictions[-idx] - 
                                    DGP$counterfactual[-idx])
  results[rep,'mspe_prod'] <- mse(Product.out$predictions[-idx] - 
                                   DGP$counterfactual[-idx])
  results[rep,'mspe_multi_norm'] <- mse(Multivariate_norm.out$predictions[-idx] - 
                                         DGP$counterfactual[-idx])
  results[rep,'mspe_prod_norm'] <- mse(Product_norm.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  
  results[rep,'rho_Synth'] <- abs(.rho - treatEff(DGP$factual,preds_synth, 
                                                  p_int - 1999))
  results[rep,'rho_LASSO'] <- abs(.rho - treatEff(DGP$factual,preds_lasso, 
                                                  p_int - 1999))
  results[rep,'rho_Gaussian'] <- abs(.rho - treatEff(DGP$factual,
                                                 Gaussian.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_Epanechnikov'] <- abs(.rho - treatEff(DGP$factual,
                                                 Epanechnikov.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_Tri-cube'] <- abs(.rho - treatEff(DGP$factual,
                                                      Tri_cube.out$predictions %>% unname(), 
                                                      p_int - 1999))
  results[rep,'rho_Gaussian_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                     Gaussian_norm.out$predictions %>% unname(), 
                                                     p_int - 1999))
  results[rep,'rho_Epanechnikov_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                         Epanechnikov_norm.out$predictions %>% unname(), 
                                                         p_int - 1999))
  results[rep,'rho_Tri-cube_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                     Tri_cube_norm.out$predictions %>% unname(), 
                                                     p_int - 1999))
  results[rep,'rho_multi'] <- abs(.rho - treatEff(DGP$factual,
                                                     Multivariate.out$predictions %>% unname(), 
                                                     p_int - 1999))
  results[rep,'rho_prod'] <- abs(.rho - treatEff(DGP$factual,
                                                 Product.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_multi_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                     Multivariate_norm.out$predictions %>% unname(), 
                                                     p_int - 1999))
  results[rep,'rho_prod_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                       Product_norm.out$predictions %>% unname(), 
                                                       p_int - 1999))
  
  results[rep,'length'] <- n
  
  # Track progress
  progress$tick()
  
}

# Save results
if ('rho_Synth' %in% colnames(results)){
  saveRDS(results, paste0(here(),'/data/simulations/kernel_rho[',SEED,'].rds'))
} else {
  saveRDS(results, paste0(here(),'/data/simulations/kernel_standard[',SEED,'].rds'))
}


# Outputs -----------------------------------------------------------------

# Load data
data  <- readRDS(paste0(here(),'/data/simulations/kernel_standard[',SEED,'].rds'))

# Re-arrange data
d.plots <- data %>% 
  as.data.frame() %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         Gaussian = mspe_Gaussian, Epanechnikov = mspe_Epanechnikov, 
         `Tri-Cube` = `mspe_Tri-cube`, Gaussian_norm = mspe_Gaussian_norm,
         Epanechnikov_norm = mspe_Epanechnikov_norm, `Tri-Cube_norm` = 
           `mspe_Tri-cube_norm`, Multivariate = mspe_multi,
         Multivariate_norm = mspe_multi_norm, Product = mspe_prod,
         Product_norm = mspe_prod_norm) %>% 
  select(Statistic, Value, everything())

# Plots
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'), 
        file = paste0('kernel_comp'), aspect.ratio = c(6,6))


tbl <- matrix(NA, 12, 4); colnames(tbl) <- c('Method','Normalized','Mean RMSPE','Variance')
tbl[,1] <- c(names(d.plots)[3:4], rep(names(d.plots)[5:7], 2), rep(names(d.plots)[11:12],2))
tbl[,2] <- c('Yes', 'No', rep(c('No','Yes'), each = 3), rep(c('No','Yes'),each=2))
tbl[1:12,3] <- d.plots %>% select(3:14) %>% apply(2, mean) %>% round(3)
tbl[1:12,4] <- d.plots %>% select(3:14) %>% apply(2, var) %>% round(3)

saveRDS(tbl, paste0(here(),'/data/simulations/tbl_kernel_standard[',SEED,'].rds'))

print(xtable(tbl,type='latex',caption='Statistics of Kernel Comparison',
             label='tab:kernel_comp',align=c('l','l',rep('r',3)),
             table.placement='b',
             hline.after = c(3)),
      include.rownames=FALSE,
      file=paste0(here(),'/../Paper/Tables/kernel_comp.tex'))


pdf(paste0(here(),'/../Paper/Figures/Simulation/kernel_series_length.pdf'),
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