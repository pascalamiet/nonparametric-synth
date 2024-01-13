##' @title  Synthetic Control Groups -- Playground
##' @author Pascal Amiet
##' @date   24.07.2023


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('Synth','tidyverse','here','glmnet','rpart','rpart.plot','randomForest') 
for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}


# Synth Package -----------------------------------------------------------

##' Get the Basque Country data
data('basque')

dataprep.out <- dataprep(
  foo = basque,
  predictors = c("school.illit", "school.prim", "school.med",
                 "school.high", "school.post.high", "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1960:1969, "mean"),
    list("sec.agriculture",seq(1961, 1969, 2), "mean"),
    list("sec.energy",seq(1961, 1969, 2), "mean"),
    list("sec.industry",seq(1961, 1969, 2), "mean"),
    list("sec.construction",seq(1961, 1969, 2), "mean"),
    list("sec.services.venta",seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

##' One-dimensional
dataprep.out <- dataprep(
  foo = basque,
  predictors = c('gdpcap'),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

synth.out = synth(dataprep.out,'BFGS')

path.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(0, 12), Legend = c("Basque country",
                                      "synthetic Basque country"), 
          Legend.position = "bottomright")

gaps.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "gap in real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(-1.5, 1.5), Main = NA)

dataprep.out$Y0plot %*% synth.out$solution.w



# Own 1-d approach --------------------------------------------------------

##' Re-formatting
basque %>% 
  select(regionname, year, gdpcap) %>% 
  pivot_wider(names_from = regionname, values_from = gdpcap) %>% 
  as.data.frame() %>% 
  arrange(year) %>% 
  select(-year) -> data
  
rownames(data) <- basque %>% pull(year) %>% unique() %>% sort()

##' Train/test split
train = data[1:which(rownames(data) == 1969),]
test = data[(which(rownames(data) == 1969)+1):nrow(data),]

y_train = train %>% pull('Basque Country (Pais Vasco)')
X_train = train %>% select(-'Basque Country (Pais Vasco)') %>% data.matrix()

y_test = test %>% pull('Basque Country (Pais Vasco)')
X_test = test %>% select(-'Basque Country (Pais Vasco)') %>% data.matrix

##' Fitting
cv.model <- cv.glmnet(X_train,y_train,alpha=1,intercept=F)
best.lambda <- cv.model$lambda.min
model <- glmnet(X_train,y_train,alpha=1,lambda=best.lambda,intercept=F)
wgt <- glmnet(X_train,y_train,alpha=1,lambda=best.lambda,intercept=F) %>% coef()

data.matrix(wgt) %>% as.data.frame() %>%  rename(coef = s0) %>% 
  slice(2:n()) %>% t() -> wgt

predict(model,rbind(X_train,X_test)) %>% as.data.frame() %>% 
  pull(s0) -> synth.basque

basePlot <- function(){
  plot(basque %>% pull(year) %>% unique() %>% sort(),
       c(y_train,y_test),
       type='l',
       ylim=c(0,12),
       ylab = 'gdpcap',
       xlab = 'years')
  abline(v=1969,lty = 'dashed')
  for (i in 1:ncol(X_train)){
    .s = rbind(X_train, X_test)[,i] %>% unname()
    lines(basque %>% pull(year) %>% unique() %>% sort(), .s, 
          type = 'l', col = rgb(0, 0, 0, 0.1))
  }
}

plot(basque %>% pull(year) %>% unique() %>% sort(),
     c(y_train,y_test),
     type='l',
     ylim=c(0,12),
     ylab = 'gdpcap',
     xlab = 'years')
lines(basque %>% pull(year) %>% unique() %>% sort(),synth.basque,type='l',
      col='red')
abline(v=1969,lty = 'dashed')

##' Add predictor lines
for (i in 1:ncol(X_train)){
  .s = rbind(X_train, X_test)[,i] %>% unname()
  lines(basque %>% pull(year) %>% unique() %>% sort(), .s, 
        type = 'l', col = rgb(0, 0, 0, 0.1))
}

##' Trees
.train = as.data.frame(cbind(y_train, X_train))
names(.train) = c('y',paste0('x',1:ncol(X_train)))
.test = as.data.frame(X_test)
names(.test) = paste0('x',1:ncol(X_test))
tree <- rpart(y ~ ., data=.train,
              control = rpart.control(cp = 0.0, minsplit = 2))

predict(tree, rbind(.train[,2:ncol(.train)],.test)) -> synth.tree

lines(basque %>% pull(year) %>% unique() %>% sort(), synth.tree, type = 'l',
      col = 'green')

##' @note Problem with extrapolation of trees if the domain is not the same.
##' Thought of demeaning, but then that takes away a lot of the variance. How
##' could this problem be solved?
##' However, as long as the domain stays the same, it could work pretty well.

##' Attempt to deduct trend
d.trend <- as.data.frame(cbind(y_train, c(1:length(y_train))))
names(d.trend) <- c('y','p')
m.trend <- lm(y ~ p, data = d.trend)
trend <- predict(m.trend, data.frame(p = c(1:nrow(data))))
lines(basque %>% pull(year) %>% unique() %>% sort(),
      trend,
      col = 'black', type='l')

trend <- predict(m.trend, data.frame(p = c(1:length(y_train))))
trend.complete <- predict(m.trend, data.frame(p = c(1:nrow(data))))
y_train.new = y_train - trend
X_train.new = X_train - trend
X_test.new = X_test - trend.complete[16:43]

.train = as.data.frame(cbind(y_train.new, X_train.new))
names(.train) = c('y',paste0('x',1:ncol(X_train.new)))
.test = as.data.frame(X_test.new)
names(.test) = paste0('x',1:ncol(X_test.new))
tree <- rpart(y ~ ., data=.train,
              control = rpart.control(cp = 0.0, minsplit = 2))

predict(tree, rbind(.train[,2:ncol(.train)],.test)) + trend.complete -> synth.tree

lines(basque %>% pull(year) %>% unique() %>% sort(), synth.tree, type = 'l',
      col = 'purple')

##' @note Still far from good...
##' The problem is non-constant variance/deviance from de-trended time series.



# Own n-d approach --------------------------------------------------------



.c1 = (c(-50:50)^(2) + 50) / 250
.c2 = (c(-50:50)^(2) - 50) / 500
.c3 = c(50:-50) / 5
.c = cbind(.c1,.c2,.c3) %>% data.matrix()
.t = c(-50:50) / 5
plot(.t,type = 'l')
lines(.c1, col='red')
lines(.c2, col='red')
lines(.c3, col='red')

fit.ols <- lm(.t ~ -1 + .c)
lines(fit.ols$fitted.values, col='green')

cv.model <- cv.glmnet(.c,.t,alpha=1,intercept=F)
best.lambda <- cv.model$lambda.min
fit.lasso <- glmnet(.c,.t,alpha=1,intercept=F,lambda=best.lambda)
lines(predict(fit.lasso, .c), col='purple')










