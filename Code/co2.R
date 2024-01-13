##' @title  Synthetic Control Groups -- CO2 Data
##' @author Pascal Amiet
##' @date   07.10.2023


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


# CO2 Data ----------------------------------------------------------------

# emissions_raw <- read_csv(paste0(here(),'/data/emissions.csv'))
emissions_raw <- read_csv(paste0(here(),'/data/emissions_total.csv'))
names(emissions_raw)[4] <- 'CO2'

.l <- list()
.codes <- unique(emissions_raw$Code)[-which(is.na(unique(emissions_raw$Code)))]

# Get distinct time series
for (c in .codes){
  .e <- emissions_raw %>% 
    filter(Code == c) %>% 
    mutate(CO2 = CO2 / 1000000000) %>% 
    select(Year,CO2)
  names(.e)[2] <- c
  .l[[c]] <- .e
}

suppressMessages({
  # Join time series
  emissions <- .l %>%
    reduce(full_join) %>% 
    arrange(Year) %>% 
    filter(Year > 1949)
  # Removing NAs
  emissions <- emissions[colSums(is.na(emissions)) == 0]
})

# Splits
.data = emissions[2:ncol(emissions)]
.train = .data[1:60,]
.test = .data[61:72,]

.cty <- 'IND'

# LASSO

cv.model <- cv.glmnet(as.matrix(.train[,-which(names(.train)==.cty)]),
                      .train[[.cty]],
                      alpha=1,
                      intercept=F)

best.lambda <- cv.model$lambda.min
model <- glmnet(as.matrix(.train[,-which(names(.train)==.cty)]),
                .train[[.cty]],
                alpha=1,
                lambda=best.lambda,
                intercept=F)

predict(model, as.matrix(.data[,-which(names(.data)==.cty)])) %>% 
  as.data.frame() %>% pull(s0) -> synth.lasso


# Trees

# Single Tree
form <- as.formula(paste(.cty,'~ .'))
tree <- rpart(form, data=.train,
              control = rpart.control(cp = 0.0, minsplit = 2))
predict(tree, .data[,-which(names(.data)==.cty)]) -> synth.tree
# lines(dates, synth.tree, type = 'l', col = 'green')

# Random Forest
rf <- randomForest(form, data=.train)
predict(rf, .data) -> synth.rf

# Bagged Tree
bag <- baggedTree(.train, .cty, B=200)
synth.bag <- predict.baggedTree(bag, .data)

# Boosted Trees
boost <- boostedTree(.train, .cty, B=200, info=F)
synth.boost <- predict.boostedTree(boost, .data)


# GAM 

##' @note There are easily too many predictors for the amount of observations
##' when working with GAMs. Hence, we would have to implement a clever form
##' of model selection.

# Select States closest to CHN
closeCHN <- c('CHN','IND','PAK','RUS','KOR','VNM','THA')
closeCHN <- c('CHN','IND','RUS','OWID_WRL','KOR','VNM','BRA')

# Fit model
form2 <- wrapFormula(as.formula(paste(.cty,'~ .')), data = .train[closeCHN], 
                    wrapString = 's(*)')
gam <- gam(form2, data = .train[closeCHN], select=T)
synth.gam <- predict(gam, .data)



# MARS

mars <- earth(form, data = .train)
synth.mars <- predict(mars, .data)




# Plot --------------------------------------------------------------------

pdf(paste0(here(),'/../Paper/Figures/co2_wide.pdf'),
    height = 6,
    width = 12)

plot(emissions$Year, emissions[[.cty]],
     type='l',
     ylim = c(0,5),
     xlim = c(1990,2021),
     ylab = expression(paste('Tons of CO'[2],' Emission (Billions)')),
     xlab = 'Years')
abline(v = 2010, lty = 'dashed')
for (i in names(emissions[2:ncol(emissions)])){
  .s = emissions %>% pull(i)
  lines(emissions$Year, .s, 
        type = 'l', col = rgb(0, 0, 0, 0.1))
}

# Predictions
lines(emissions$Year[61:72], synth.lasso[61:72], type = 'l', col = .c2[1],lwd=3)
lines(emissions$Year[61:72], synth.rf[61:72],    type = 'l', col = .c2[2],lwd=3)
lines(emissions$Year[61:72], synth.bag[61:72],   type = 'l', col = .c2[3],lwd=3)
lines(emissions$Year[61:72], synth.gam[61:72],   type = 'l', co = .c2[4],lwd=3)
lines(emissions$Year[61:72], synth.mars[61:72],   type = 'l', co = .c2[5],lwd=3)
# lines(emissions$Year[61:72], synth.boost[61:72], type = 'l', col = .c2[6],lwd=3)
# lines(emissions$Year[61:72], synth.tree[61:72], type = 'l', col = .c2[7],lwd=3)

# Add legend
legend(x = 'topright', 
       legend = c('India','LASSO','Random Forest','Bagged Trees', 
                  'GAM', 'MARS'), 
       col = c('black',.c2[1:5]), 
       lty = c(1,1,1,1,1),
       lwd = c(1,rep(3,5)))

dev.off()


pdf(paste0(here(),'/../Paper/Figures/co2_gam.pdf'),
    height = 6,
    width = 9)

par(mfrow=c(2,3))
plot(gam,rug=T)
par(mfrow=c(1,1))

dev.off()



