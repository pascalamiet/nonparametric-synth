##' @title  Synthetic Control Groups -- Cross-Entropy
##' @author Pascal Amiet
##' @date   12.12.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('here','dplyr','Synth') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')


# Cross-Entropy Function --------------------------------------------------

klDiv(data, 'y', 1, 'years', 2009, c('x1','x2','x3'))

x <- klDiv(X = data, k = 'y', id = 'id', id.interest = 1,
           periods = 'year', intervention = 2009, 
           controls = c('x1','x2','x3'), h = 0.1,
           discrete = F,
           method = 'Aggregate', true_cf = DGP$counterfactual,
           norm.wgt = T)


# Entropy Plot ------------------------------------------------------------

data <- readRDS(paste0(here(), '/npsynth/sim_data.rds'))

.var <- 'x2'
d.plot <- data.frame(x1 = data[data$id == 1, .var],
                     x2 = data[data$id == 4, .var],
                     x3 = data[data$id == 3, .var])

# Save plot
pdf(paste0(here(),'/../Paper/Figures/entropy_ill.pdf'),
    height = 6,
    width = 10) 

ggplot(d.plot) +
  geom_density(aes(x=x1), color=.c2[5], fill=.c2[9], alpha = .5) +
  geom_density(aes(x=x3), color=.c2[10], fill=.c2[6], alpha = .5) +
  geom_density(aes(x=x2), color=.c2[11], fill=.c2[7], alpha = .5) +
  xlim(-7,24) +
  theme_void() +
  theme(legend.position = "none")

dev.off()


# Performance -------------------------------------------------------------

# Setup
SEED <- 1337
set.seed(SEED)
reps <- 200
results <- matrix(NA_real_,reps,19)
colnames(results) <- c(paste0(rep(c('mse_','mspe_'),each=6),
                              rep(c('Synth','LASSO','disc','cont','disc_norm','cont_norm'),2)),
                       'rho_Synth','rho_LASSO',
                       'rho_disc','rho_cont','rho_disc_norm','rho_cont_norm',
                       'length')
functional.form <- 'New' 
progress <- progress_bar$new(total = reps)
plot = F

# Simulation
for (rep in 1:reps){
  
  # Parameters
  N <- 200
  n <- ceiling(runif(1,10,100))
  n_x <- 3
  p_int <- 1999 + floor(2/3 * n)
  
  # Simulate data
  temp.obj <- DGPcerulli(N, n, functional.form = functional.form)
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
    for (bool.disc in c(T, F)){
      inSilence({
        x <<- klDiv(X = data, k = 'y', id = 'id', id.interest = 1,
                    periods = 'year', intervention = p_int, 
                    controls = c('x1','x2','x3'), h = 0.1,
                    discrete = bool.disc, method = 'Aggregate', 
                    true_cf = DGP$counterfactual, norm.wgt = bool.norm,
                    standardize = T)
      })
      .n <- ifelse(bool.disc,
                   paste0('KL_disc',ifelse(bool.norm, '_norm.out','.out')),
                   paste0('KL_cont',ifelse(bool.norm, '_norm.out','.out')))
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
  results[rep,'mse_disc'] <- mse(KL_disc.out$predictions[idx] - 
                                   DGP$counterfactual[idx])
  results[rep,'mse_cont'] <- mse(KL_cont.out$predictions[idx] - 
                                   DGP$counterfactual[idx])
  results[rep,'mse_disc_norm'] <- mse(KL_disc_norm.out$predictions[idx] - 
                                        DGP$counterfactual[idx])
  results[rep,'mse_cont_norm'] <- mse(KL_cont_norm.out$predictions[idx] - 
                                        DGP$counterfactual[idx])

  results[rep,'mspe_Synth'] <- mse(preds_synth[-idx] - 
                                     DGP$counterfactual[-idx])
  results[rep,'mspe_LASSO'] <- mse(preds_lasso[-idx] - 
                                    DGP$counterfactual[-idx])
  results[rep,'mspe_disc'] <- mse(KL_disc.out$predictions[-idx] - 
                                   DGP$counterfactual[-idx])
  results[rep,'mspe_cont'] <- mse(KL_cont.out$predictions[-idx] - 
                                   DGP$counterfactual[-idx])
  results[rep,'mspe_disc_norm'] <- mse(KL_disc_norm.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  results[rep,'mspe_cont_norm'] <- mse(KL_cont_norm.out$predictions[-idx] - 
                                        DGP$counterfactual[-idx])
  
  results[rep,'rho_Synth'] <- abs(.rho - treatEff(DGP$factual,preds_synth, 
                                                  p_int - 1999))
  results[rep,'rho_LASSO'] <- abs(.rho - treatEff(DGP$factual,preds_lasso, 
                                                  p_int - 1999))
  results[rep,'rho_disc'] <- abs(.rho - treatEff(DGP$factual,
                                                 KL_disc.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_cont'] <- abs(.rho - treatEff(DGP$factual,
                                                 KL_cont.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_disc_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                 KL_disc_norm.out$predictions %>% unname(), 
                                                 p_int - 1999))
  results[rep,'rho_cont_norm'] <- abs(.rho - treatEff(DGP$factual,
                                                 KL_cont_norm.out$predictions %>% unname(), 
                                                 p_int - 1999))

  results[rep,'length'] <- n
  
  # Track progress
  progress$tick()
  
}

# Save results
saveRDS(results, paste0(here(),'/data/simulations/entropy_standard[',SEED,'].rds'))


# Outputs -----------------------------------------------------------------

# Load data
data  <- readRDS(paste0(here(),'/data/simulations/entropy_standard[',SEED,'].rds'))

# Re-arrange data
d.plots <- data %>% 
  as.data.frame() %>% 
  select(starts_with('mspe'), length) %>% 
  mutate(Statistic = 'RMSPE',
         across(starts_with('mspe'), ~ sqrt(.x))) %>% 
  rename(Value = length, Synth = mspe_Synth, LASSO = mspe_LASSO, 
         Discrete = mspe_disc, Continuous = mspe_cont, 
         Discrete_norm = mspe_disc_norm, Continuous_norm = mspe_cont_norm) %>% 
  select(Statistic, Value, Synth, LASSO, Discrete, Continuous, Discrete_norm,
         Continuous_norm)

# Plots
simPlot(d.plots, 'Number of Periods', metrics = c('RMSPE'))
simPlot(d.plots %>% select(-all_of(c('Discrete','Continuous'))), 
        'Number of Periods', metrics = c('RMSPE'))

simPlot(d.plots %>% select(-all_of(c('Discrete','Continuous'))), 
        'Number of Periods', metrics = c('RMSPE'), 
        file = 'entropy_comp', aspect.ratio = c(6,6))

# Table
tbl <- matrix(NA, 6, 4); colnames(tbl) <- c('Method','Normalized','Mean RMSPE','Variance')
tbl[,1] <- c(names(d.plots)[3:4], rep(names(d.plots)[5:6], 2))
tbl[,2] <- c('Yes', 'No', rep(c('No','Yes'), each = 2))
tbl[1:4,3] <- d.plots %>% select(3:6) %>% apply(2, mean) %>% round(3)
tbl[1:4,4] <- d.plots %>% select(3:6) %>% apply(2, var) %>% round(3)
tbl[5:6,3] <- d.plots %>% select(7:8) %>% apply(2, mean) %>% round(3)
tbl[5:6,4] <- d.plots %>% select(7:8) %>% apply(2, var) %>% round(3)

saveRDS(tbl, paste0(here(),'/data/simulations/tbl_entropy_standard[',SEED,'].rds'))

print(xtable(tbl,type='latex',caption='Statistics of Entropy Comparison',
             label='tab:entropy_comp',align=c('l','l',rep('r',3)),
             table.placement='b',
             hline.after = c(3)),
      include.rownames=FALSE,
      file=paste0(here(),'/../Paper/Tables/entropy_comp.tex'))


pdf(paste0(here(),'/../Paper/Figures/Simulation/entropy_series_length.pdf'),
    height = 5, width = 8) 

plot(rep(d.plots$Value, 6), stack(d.plots %>% select(c(3:8)))$y,
     xlab = 'Time Series Length', ylab = 'RMSPE', cex = 0.5, pch = 16,
     col = rgb(0, 0, 0, 0.2))
.m <- names(d.plots)[3:8]
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


