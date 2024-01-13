##' @title  Synthetic Control Groups -- Simulation w/out Covariates
##' @author Pascal Amiet
##' @date   04.11.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('tidyverse','here','progress','Synth','xtable') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')

SEED <- 1337
set.seed(SEED)


# With Covariates ---------------------------------------------------------

# Load data
d.dissimilarity <- readRDS(paste0(here(),'/data/simulations/dissimilarity[',SEED,'].rds'))
d.kernel <- readRDS(paste0(here(),'/data/simulations/kernel_standard[',SEED,'].rds'))
d.entropy <- readRDS(paste0(here(),'/data/simulations/entropy_standard[',SEED,'].rds'))

data <- d.dissimilarity[,1:16] %>% 
  cbind(d.kernel[,c(3:5,8:42)]) %>% 
  cbind(d.entropy[,c(3:6,9:12,15:19)]) %>% 
  as.data.frame()

# MSE all
d.plots <- data %>% 
  select(starts_with('mse'), length) %>% 
  mutate(Statistic = 'RMSE',
         across(starts_with('mse'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mse_Synth, LASSO = mse_LASSO, 
         L1 = mse_Manhattan, L2 = mse_Euclidean, Correlation = mse_Correlation,
         Gaussian = mse_Gaussian, Epanechnikov = mse_Epanechnikov,
         `Tri-Cube` = `mse_Tri-cube`, `MV Kernel` = mse_multi, 
         `Product Kernel` = mse_prod,
         `KL Discrete` = mse_disc, `KL Continuous` = mse_cont) %>% 
  select(Statistic, Value, Synth, LASSO, Gaussian, Epanechnikov, `Tri-Cube`,
         `MV Kernel`, `Product Kernel`,
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('RMSE'))
simPlot(d.plots, 'Number of Periods', metrics = c('RMSE'),
        file = 'mse_all', aspect.ratio = c(5,15), line.plot = F)


# MSE all normalized
d.plots <- data %>% 
  select(starts_with('mse'), length) %>% 
  mutate(Statistic = 'RMSE',
         across(starts_with('mse'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mse_Synth, LASSO = mse_LASSO, 
         L1 = mse_Manhattan_norm, L2 = mse_Euclidean_norm, 
         Correlation = mse_Correlation_norm,
         Gaussian = mse_Gaussian_norm, Epanechnikov = mse_Epanechnikov_norm,
         `Tri-Cube` = `mse_Tri-cube_norm`,`MV Kernel` = mse_multi_norm, 
         `Product Kernel` = mse_prod_norm, `KL Discrete` = mse_disc_norm, 
         `KL Continuous` = mse_cont_norm) %>% 
  select(Statistic, Value, Synth, LASSO, Gaussian, 
         Epanechnikov, `Tri-Cube`,`MV Kernel`, `Product Kernel`, 
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('RMSE'))
simPlot(d.plots, 'Number of Periods', metrics = c('RMSE'),
        file = 'mse_all_norm', aspect.ratio = c(5,15), line.plot = F)

# MSPE all
d.plots <- data %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         L1 = mspe_Manhattan, L2 = mspe_Euclidean, Correlation = mspe_Correlation,
         Gaussian = mspe_Gaussian, Epanechnikov = mspe_Epanechnikov,
         `Tri-Cube` = `mspe_Tri-cube`, `MV Kernel` = mspe_multi, 
         `Product Kernel` = mspe_prod, `KL Discrete` = mspe_disc, 
         `KL Continuous` = mspe_cont) %>% 
  select(Statistic, Value, Synth, LASSO, Gaussian, 
         Epanechnikov, `Tri-Cube`,`MV Kernel`, `Product Kernel`,  
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'),
        file = 'mspe_all', aspect.ratio = c(5,15), line.plot = F)


# MSPE all normalized
d.plots <- data %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         L1 = mspe_Manhattan_norm, L2 = mspe_Euclidean_norm, 
         Correlation = mspe_Correlation_norm,
         Gaussian = mspe_Gaussian_norm, Epanechnikov = mspe_Epanechnikov_norm,
         `Tri-Cube` = `mspe_Tri-cube_norm`, `MV Kernel` = mspe_multi_norm, 
         `Product Kernel` = mspe_prod_norm, `KL Discrete` = mspe_disc_norm, 
         `KL Continuous` = mspe_cont_norm) %>% 
  select(Statistic, Value, Synth, LASSO, Gaussian, 
         Epanechnikov, `Tri-Cube`, `MV Kernel`, `Product Kernel`, 
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'),
        file = 'mspe_all_norm', aspect.ratio = c(5,15), line.plot = F)

# Table
t.dissimilarity <- readRDS(paste0(here(),'/data/simulations/tbl_dissimilarity[',SEED,'].rds'))
t.kernel <- readRDS(paste0(here(),'/data/simulations/tbl_kernel[',SEED,'].rds'))
t.entropy <- readRDS(paste0(here(),'/data/simulations/tbl_entropy[',SEED,'].rds'))

t.dissimilarity_std <- readRDS(paste0(here(),'/data/simulations/tbl_dissimilarity_standard[',SEED,'].rds'))
t.kernel_std <- readRDS(paste0(here(),'/data/simulations/tbl_kernel_standard[',SEED,'].rds'))
t.entropy_std <- readRDS(paste0(here(),'/data/simulations/tbl_entropy_standard[',SEED,'].rds'))

tbl <- rbind(t.dissimilarity, t.kernel[3:12,], t.entropy[3:6,],
             t.dissimilarity_std[3:8,], t.kernel_std[3:12,], t.entropy_std[3:6,])
tbl[c(3:8,23:28),'Method'] <- rep(c('L-1 Norm','L-2 Norm','Correlation'),4)
tbl[c(9:18,29:38),'Method'] <- paste0(tbl[c(9:18,29:38),'Method'],' Kernel')
tbl[c(19:22,39:42),'Method'] <- paste0(tbl[c(19:22,39:42),'Method'],' KL')

tbl <- tbl %>% as.data.frame() %>% 
  mutate(Standardized = c(rep('No',22),rep('Yes',20))) %>% 
  select(Method, Standardized, everything()) %>% as.matrix()

print(xtable(tbl,type='latex',caption='Statistics of the First Simulation',
             label='tab:sim_comp',align=c('l','l',rep('r',4)),
             table.placement='b',
             hline.after = c(3)),
      include.rownames=FALSE,
      file=paste0(here(),'/../Paper/Tables/sim_comp.tex'))

# Best methods
d.plots <- data %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         Gaussian = mspe_Gaussian_norm, Epanechnikov = mspe_Epanechnikov_norm,
         `Tri-Cube` = `mspe_Tri-cube_norm`,`MV Kernel` = mspe_multi_norm, 
         `Product Kernel` = mspe_prod_norm, `KL Discrete` = mspe_disc_norm, 
         `KL Continuous` = mspe_cont_norm) %>% 
  select(Statistic, Value, Synth, Gaussian, 
         Epanechnikov, `Tri-Cube`, `MV Kernel`, `Product Kernel`, 
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'),
        file = 'best_norm', aspect.ratio = c(7,12), line.plot = F)

# Time series length
pdf(paste0(here(),'/../Paper/Figures/Simulation/sim_length.pdf'),
    height = 5, width = 8) 

plot(rep(d.plots$Value, 6), stack(d.plots %>% select(c(3:8)))$y,
     xlab = 'Time Series Length', ylab = 'RMSPE', cex = 0.5, pch = 16,
     col = rgb(0, 0, 0, 0.2))
.m <- names(d.plots)[3:10]
for (i in 1:length(.m)){
  .val <- paste0('`',.m[i],'`')
  .frm <- as.formula(paste(.val,'~ Value'))
  .reg <- lm(.frm, data = d.plots)
  lines(11:100, .reg$coefficients[1] + .reg$coefficients[2] * c(11:100),
        col = .c2[i])
}
legend(x = 'topright', legend = .m, 
       col = .c2[1:length(.m)], 
       lty= rep(1, length(.m)),
       bg = rgb(1, 1, 1, 0.7))

dev.off()

# Plots of rho deviance
d.plots <- data %>% 
  select(starts_with('rho'), length) %>% 
  mutate(Statistic = 'True - Estimated CTE') %>% 
  rename(Value = length, Synth = rho_Synth, LASSO = rho_LASSO, 
         Gaussian = rho_Gaussian_norm, Epanechnikov = rho_Epanechnikov_norm,
         `Tri-Cube` = `rho_Tri-cube_norm`,`MV Kernel` = rho_multi_norm, 
         `Product Kernel` = rho_prod_norm, 
         `KL Discrete` = rho_disc_norm, `KL Continuous` = rho_cont_norm) %>% 
  select(Statistic, Value, Synth, 
         Gaussian, Epanechnikov, `Tri-Cube`, `MV Kernel`, `Product Kernel`, 
         `KL Discrete`, `KL Continuous`)

simPlot(d.plots, 'Number of Periods', metrics = c('True - Estimated CTE'), line.plot = F)
simPlot(d.plots, 'Number of Periods', metrics = c('True - Estimated CTE'),
        file = 'rho', aspect.ratio = c(7,12), line.plot = F)


# DGP Plots w/out Covariates ----------------------------------------------

set.seed(19)
DGP(100, c('cubic',rep('linear',4)), gamma = 4, max.diff = 200)

pdf(paste0(here(),'/../Paper/Figures/Simulation/sim2_ex_different.pdf'),
    height = 5, width = 5) 
DGP(100, c('cubic',rep('linear',4)), gamma = 1, max.diff = 200)
dev.off()

pdf(paste0(here(),'/../Paper/Figures/Simulation/sim2_ex_same.pdf'),
    height = 5, width = 5) 
DGP(100, rep('cubic',10), gamma = 3, max.diff = 200)
dev.off()


# For Appendix
set.seed(19)
for (i in 1:6){
  
  poss.forms <- c('linear','quadratic','cubic','oscillating')
  forms <- sample(poss.forms, 5, replace = T)
  
  pdf(paste0(here(),'/../Paper/Figures/Simulation/sim2_ex_app_',i,'.pdf'),
      height = 5, width = 5) 
  DGP(100, forms, gamma = 3, max.diff = runif(1, 20, 200))
  dev.off() 
  
}


# Comparison w/out Covariates ---------------------------------------------

# General Comparison
for (diff.forms in c(F, T)){
    
  rep <- 200
  set.seed(SEED)
  pb <- progress_bar$new(total = rep)
  
  for (i in 1:rep){
    
    groups <- rpois(1,5); groups <- ifelse(groups > 0, groups, 1)
    len <- round(runif(1,50,200),0)
    # Specify the forms
    poss.forms <- c('linear','quadratic','cubic','oscillating')
    if (diff.forms) poss.forms <- c('quadratic','cubic','oscillating')
    treat.form <- poss.forms[round(runif(1,1,length(poss.forms)),0)]
    if (diff.forms) poss.forms <- c('linear')
    forms <- c(treat.form, poss.forms[round(runif(groups,1,length(poss.forms)),0)])
    # Generate data
    gamma <- runif(1, 0.2, 2)
    max.diff <- runif(1, 50, 200); diff.total <- runif(1, 5, 15)
    X <- DGP(len, forms = forms, plot = F, gamma = gamma, max.diff = max.diff,
             diff.total = diff.total)
    # Test models
    inSilence({
      .temp <- testModels(X$data, t = floor(len * 2/3), 
                           models = c('LASSO','rf','GAM','MARS','SVM',
                                      'ProjPursuit','uniKernel','multiKernel',
                                      'entropy'), 
                           plot = F)
      stats <<- .temp$stats
      preds <<- .temp[-length(.temp)]
    })
    stats$Value <- i; stats <- stats %>% select(Statistic, Value, everything())
    if (i == 1){ Stats <- stats } else { Stats <- rbind(Stats, stats) }
    
    # Show progress
    pb$tick()
    
  }
  
  # Save files
  if (diff.forms){
    write_rds(Stats, paste0(here(),'/data/simulations/stats_diff[',SEED,'].rds'))
  } else {
    write_rds(Stats, paste0(here(),'/data/simulations/stats_same[',SEED,'].rds'))
  }
}


stats_same <- readRDS(paste0(here(),'/data/simulations/stats_same[',SEED,'].rds'))
stats_diff <- readRDS(paste0(here(),'/data/simulations/stats_diff[',SEED,'].rds'))

# Plot
d.plots <- stats_same %>% 
  filter(Statistic == 'RMSE') %>% 
  mutate(Statistic = 'RMSPE')

simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'), line.plot = F)



# Time Series Length

##' @note This simulation tests the performance of different models given 
##' different time series lengths. For each length of the time series the 
##' process is repeated multiple times (can be set via the `rep` parameter)
##' such that the variance can be decreased. The maximum length of the time 
##' series can be chosen differently, although after a certain amount the 
##' changes will be negligible. The functional forms can be adjusted as well. 

# Simulation
min <- 8
max <- 100
rep <- 10
forms <- c('cubic',rep('linear',5),rep('quadratic',3),'cubic')

set.seed(SEED)
pb <- progress_bar$new(total = (max-min))

for (l in min:max){
  # Generate data multiple times to decrease variance
  for (i in 1:rep){
    gamma <- rnorm(1, 0.3, 0.2)
    max.diff <- runif(1, 6, 20)
    X <- DGP(l, forms = forms, plot = F, gamma = gamma, max.diff = max.diff)
    # Test models
    s <- testModels(X$data, t = floor(l * 2/3), 
                    models = c('LASSO','rf','baggedTrees','boostedTrees','MARS',
                               'GAM','uniKernel','entropy'), 
                    plot = F)
    if (i == 1){ stats <- s } else { stats <- rbind(stats, s)}
  }
  # Save error stats
  stats <- stats %>% group_by(Statistic) %>% 
    summarise_at(names(stats)[2:ncol(stats)], mean, na.rm=T)
  stats$Value <- l; stats <- stats %>% select(Statistic, Value, everything())
  if (l == min){ Stats <- stats } else { Stats <- rbind(Stats, stats) }
  # Show progress
  pb$tick()
}

# Export data
write_rds(Stats, paste0(here(),'/data/simulation/length[',SEED,'].rds'))

# Plots
simPlot(Stats, file = paste0('length[',SEED,']'))







# Quick Results

# Load previously computed statistics
stats_same <- readRDS(paste0(here(),'/data/simulations/stats_same[',SEED,'].rds'))
stats_diff <- readRDS(paste0(here(),'/data/simulations/stats_diff[',SEED,'].rds'))
  
# Plot
d.plots <- stats_same %>% 
  filter(Statistic == 'RMSE') %>% 
  mutate(Statistic = 'RMSPE') %>% 
  rename(RF = rf, Kernel = uniKernel, `MV Kernel` = multiKernel,
         KL = entropy, PPR = ProjPursuit) %>% 
  select(Statistic, Value, LASSO, RF, GAM, MARS, PPR, SVM, Kernel, KL)

# simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'), line.plot = F)
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'),
        file = 'stats_same', aspect.ratio = c(7,12), line.plot = F)

# d.plots <- stats_diff %>% 
#   filter(Statistic == 'RMSE') %>% 
#   mutate(Statistic = 'RMSPE') %>% 
#   rename(RF = rf, Kernel = uniKernel, KL = entropy) %>% 
#   select(Statistic, Value, LASSO, RF, GAM, MARS, Kernel, KL)
# simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'), line.plot = F)
# simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'),
#         file = 'stats_diff', aspect.ratio = c(6,6), line.plot = F)

stats_same <- stats_same %>% select(-multiKernel)
stats_diff <- stats_diff %>% select(-multiKernel)

.x <- apply(stats_same %>% filter(Statistic == 'RMSE') %>% 
        select(-all_of(c('Statistic','Value'))), 
      2, function(x) mean(x, na.rm=T)) 

.a <- apply(stats_same %>% filter(Statistic == 'RMSE') %>% 
              select(-all_of(c('Statistic','Value'))), 
            2, function(x) var(x, na.rm=T)) 

.y <- apply(stats_diff %>% filter(Statistic == 'RMSE') %>% 
              select(-all_of(c('Statistic','Value'))), 
      2, function(x) mean(x, na.rm=T))

.b <- apply(stats_diff %>% filter(Statistic == 'RMSE') %>% 
              select(-all_of(c('Statistic','Value'))), 
            2, function(x) var(x, na.rm=T)) 

.c <- rep('-', 2 * length(.x))
.c[2*(1:length(.x))] <- paste(round(100 * (.y - .x) / .x, 3),'%')

.X <- rep(NA, 2*length(.x))
.X[2*(1:length(.x))-1] <- .x; .X[2*(1:length(.x))] <- .y

.Y <- rep(NA, 2*length(.a))
.Y[2*(1:length(.a))-1] <- .a; .Y[2*(1:length(.a))] <- .b

tbl_sim2 <- data.frame(
  Method = rep(c('LASSO','RF','GAM','MARS','SVM','PPR','Kernel','KL'), each=2),
  fform = rep(c('Same','Different'), 8),
  RMSPE = .X,
  Variance = .Y,
  Change = .c) %>% 
  rename(`Functional Form` = fform,
         `RMSPE Change` = Change)

print(xtable(tbl_sim2,type='latex',caption='Statistics of the Second Simulation',
             label='tab:sim2_comp',align=c('l','l',rep('r',4)),
             table.placement='T',
             hline.after = c(3)),
      include.rownames=FALSE,
      file=paste0(here(),'/../Paper/Tables/sim2_comp.tex'))








# Load previously computed statistics
stats.length   <- readRDS(paste0(here(),'/data/simulations/length[',SEED,'].rds'))
stats.nrGroups <- readRDS(paste0(here(),'/data/simulations/nr-groups[',SEED,'].rds'))
stats.fnForms  <- readRDS(paste0(here(),'/data/simulations/non-similar-fn-forms[',SEED,'].rds'))
stats.linCtrl  <- readRDS(paste0(here(),'/data/simulations/linear-ctrl-gr[',SEED,'].rds'))

# Plot
simPlot(stats.length, 'Time Series Length')
simPlot(stats.nrGroups, 'Number of Control Groups')
simPlot(stats.fnForms, 'Iterations')
simPlot(stats.linCtrl, 'Iterations')

