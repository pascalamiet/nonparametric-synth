##' @title  Synthetic Control Groups -- Unemployment Data
##' @author Pascal Amiet
##' @date   28.09.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('tidyverse','here','glmnet','rpart','rpart.plot','randomForest','mgcv',
         'sfsmisc','earth','xtable') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')


# Unemployment Data -------------------------------------------------------

unemployment_raw <- read_csv(paste0(here(),'/data/unemployment.csv'), skip = 1)

# Rearranging the data
.n <- unemployment_raw$`Series ID`
unemployment <- as.data.frame(t(unemployment_raw[c(4:ncol(unemployment_raw))]))
names(unemployment) <- .n

# Splits
.data = unemployment[1:349,]
.train = .data[1:228,]
.test = .data[229:349,]

# LASSO

cv.model <- cv.glmnet(as.matrix(.train[,-which(names(.train)=='NHUR')]),
                      .train$NHUR,
                      alpha=1,
                      intercept=F)

best.lambda <- cv.model$lambda.min
model <- glmnet(as.matrix(.train[,-which(names(.train)=='NHUR')]),
                .train$NHUR,
                alpha=1,
                lambda=best.lambda,
                intercept=F)

predict(model, as.matrix(.data[,-which(names(.data)=='NHUR')])) %>% 
  as.data.frame() %>% pull(s0) -> synth.lasso


# Trees

# Single Tree
tree <- rpart(NHUR ~ ., data=.train,
              control = rpart.control(cp = 0.0, minsplit = 2))
predict(tree, .data[,-which(names(.data)=='NHUR')]) -> synth.tree
# lines(dates, synth.tree, type = 'l', col = 'green')

# Random Forest
rf <- randomForest(NHUR ~ ., data=.train)
predict(rf, .data) -> synth.rf

# Bagged Tree
bag <- baggedTree(.train, 'NHUR', B=200)
synth.bag <- predict.baggedTree(bag, .data)

# Boosted Trees
boost <- boostedTree(.train, 'NHUR', B=200)
synth.boost <- predict.boostedTree(boost, .data)


# GAM 

##' @note There are easily too many predictors for the amount of observations
##' when working with GAMs. Hence, we would have to implement a clever form
##' of model selection.

# Select States closest to NH
closeNH <- paste0(c('ME','NH','VT','MA','CT','RI','NY','PA','NJ','MD','DE','VA',
                    'WV','NC','SC','GA','FL','OH','DC','MI','IN','KY','TN','IL',
                    'WI','MN'),'UR')

# Select variables based on model performance, forward selection

# .x <- names(.data)[-which(names(.data)=='NHUR')]
# .m <- c('NHUR')
# 
# for (i in 1:25){
#   .i <- floor(0.8*nrow(.train))
#   .t <- .train[1:.i,]
#   .v <- .train[(.i+1):nrow(.train),]
#   .lmse <- c()
#   for (j in .x){
#     .t.temp <- .t[c(.m,j)]
#     .v.temp <- .v[c(.m,j)]
#     .f <- wrapFormula(NHUR ~ ., data = .t.temp, wrapString = 's(*)')
#     .g <- gam(.f, data = .t.temp, select = T)
#     .mse <- mean((.v$NHUR - predict(.g,.v.temp))^2)
#     .lmse <- c(.lmse,.mse)
#   }
#   .m <- c(.m,.x[which.min(.lmse)])
#   .x <- .x[-which.min(.lmse)]
#   print(paste0(i,'/25'))
# }
# 
# saveRDS(.m,paste0(here(),'/data/gam_model_selection_unemployment.RDS'))

modelVars <- readRDS(paste0(here(),'/data/gam_model_selection_unemployment.RDS'))

# Fit model
form <- wrapFormula(NHUR ~ ., data = .train[,closeNH], wrapString = 's(*)')
gam <- gam(form, data = .train[,closeNH], select=T)
synth.gam <- predict(gam, .data)


# MARS

mars <- earth(NHUR ~ ., data = .train)
synth.mars <- predict(mars, .data)


# Plot --------------------------------------------------------------------

# Saving settings
pdf(paste0(here(),'/../Paper/Figures/basic_np_without.pdf')) 

# Plot
dates <- as.Date(rownames(.data), format="%Y-%m-%d")
plot(dates, .data$NHUR,
     type='l',
     ylim = c(0,13),
     ylab = 'Unemployment Rate (%)',
     xlab = 'Years')
abline(v = as.Date('1995-01-01', format = "%Y-%m-%d"), lty = 'dashed')
for (i in names(.data)){
  .s = .data %>% pull(i)
  lines(dates, .s, 
        type = 'l', col = rgb(0, 0, 0, 0.1))
}

legend(x = 'topright', 
       legend = c('New Hampshire'), 
       col = c('black'), 
       lty=c(1))

dev.off()

pdf(paste0(here(),'/../Paper/Figures/basic_np.pdf')) 

plot(dates, .data$NHUR,
     type='l',
     ylim = c(0,13),
     ylab = 'Unemployment Rate (%)',
     xlab = 'Years')
abline(v = as.Date('1995-01-01', format = "%Y-%m-%d"), lty = 'dashed')
for (i in names(.data)){
  .s = .data %>% pull(i)
  lines(dates, .s, 
        type = 'l', col = rgb(0, 0, 0, 0.1))
}

# Predictions
lines(dates, synth.lasso, type = 'l', col = .c2[1])
lines(dates, synth.rf, type = 'l', col = .c2[2])
lines(dates, synth.bag, type = 'l', col = .c2[3])
lines(dates, synth.boost, type = 'l', col = .c2[4])
lines(dates, synth.gam, type = 'l', co = .c2[5])


# Add legend
legend(x = 'topright', 
       legend = c('LASSO','Random Forest','Bagged Trees', 'Boosted Trees', 'GAM'), 
       col = .c2[1:5], 
       lty=c(1,1,1,1,1))

dev.off()


pdf(paste0(here(),'/../Paper/Figures/basic_np_wide.pdf'),
    height = 6,
    width = 12) 

plot(dates, .data$NHUR,
     type='l',
     ylim = c(0,13),
     ylab = 'Unemployment Rate (%)',
     xlab = 'Years')
abline(v = as.Date('1995-01-01', format = "%Y-%m-%d"), lty = 'dashed')
for (i in names(.data)){
  .s = .data %>% pull(i)
  lines(dates, .s, 
        type = 'l', col = rgb(0, 0, 0, 0.1))
}

# Predictions
lines(dates[229:349], synth.lasso[229:349], type = 'l', col = .c2[1], lwd = 3)
lines(dates[229:349], synth.rf[229:349], type = 'l', col = .c2[2], lwd = 3)
lines(dates[229:349], synth.bag[229:349], type = 'l', col = .c2[3], lwd = 3)
# lines(dates[229:349], synth.boost[229:349], type = 'l', col = .c2[4], lwd = 3)
lines(dates[229:349], synth.gam[229:349], type = 'l', co = .c2[4], lwd = 3)
lines(dates[229:349], synth.mars[229:349], type = 'l', co = .c2[5], lwd = 3)


# Add legend
legend(x = 'topright', 
       legend = c('New Hampshire','LASSO','Random Forest','Bagged Trees', 
                  'GAM','MARS'), 
       col = c('black',.c2[1:5]), 
       lty=c(1,1,1,1,1),
       lwd = c(1,3,3,3,3,3))

dev.off()


# Statistics --------------------------------------------------------------

y <- .data$NHUR[229:349]

p.lasso <- synth.lasso[229:349] %>% unname()
p.rf <- synth.rf[229:349] %>% unname()
p.bag <- synth.bag[229:349] %>% unname()
p.boost <- synth.boost[229:349] %>% unname()
p.gam <- synth.gam[229:349] %>% unname()
p.mars <- synth.mars[229:349] %>% unname()

stats <- matrix(NA, 3, 6)
stats[ ,1] <- c('MSE','RMSE','MAE')
colnames(stats) <- c('Metric','LASSO','Random Forest','Bagged Trees','GAM','MARS')

i <- 1

for (p in list(p.lasso,p.rf,p.bag,p.gam,p.mars)){
  .mse <- mean((y-p)^2) %>% round(4)
  .rmse <- sqrt(.mse) %>% round(4)
  .mae <- mean(abs(y-p)) %>% round(4)
  stats[,i+1] <- c(.mse,.rmse,.mae)
  i <- i+1
}

stats <- as.data.frame(stats)
rownames(stats) <- stats$Metric
stats <- select(stats, -Metric)

print(xtable(stats,type='latex',caption='Unemployment Prediction Errors',
             label='tab:unemp',align=c('l',rep('r',ncol(stats))),
             table.placement='t',),
      file=paste0(here(),'/../Paper/Tables/unemp.tex'))



