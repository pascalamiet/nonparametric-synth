##' @title  Synthetic Control Groups -- Npsynth Simulation
##' @author Pascal Amiet
##' @date   22.11.2023


##' @note Replication exercise of the simulation of the paper by Cerulli (2020)
##' in R instead of Stata. The translation has been made by ChatGPT, with slight
##' adjustments made by me. In the second part, I then re-write the code in a
##' more R-friendly way to run the simulation many times to get a more 
##' consistent performance measure.

##' @url http://fmwww.bc.edu/repec/bocode/s/Simulation_npsynth.do


# Setup -------------------------------------------------------------------

##' Load required packages
.pkg = c('here','Synth') 

for (.p in .pkg){
  if(!require(.p,character.only=T)){
    install.packages(.p)
  }
}

source('functions.R')


# Stata Code --------------------------------------------------------------

##' @note In this section I replicate the Stata code, as it is given by Cerulli 
##' (2020). In the section below, I then re-write it in a more appealing
##' way for R. The first translation to R code was done using ChatGPT, with
##' slight adaptions by me.


# clear
# set obs 20
# set seed 101
# drawnorm e

set.seed(101)
n <- 50
e <- rnorm(n)

# generate x1_1=invnorm(uniform())
# generate x2_1=invnorm(uniform())
# generate x3_1=invnorm(uniform())

x1_1 <- qnorm(runif(n))
x2_1 <- qnorm(runif(n))
x3_1 <- qnorm(runif(n))

# generate x1_2=invnorm(uniform())
# generate x2_2=invnorm(uniform())
# generate x3_2=invnorm(uniform())

x1_2 <- qnorm(runif(n))
x2_2 <- qnorm(runif(n))
x3_2 <- qnorm(runif(n))

# generate x1_3=invnorm(uniform())
# generate x2_3=invnorm(uniform())
# generate x3_3=invnorm(uniform())

x1_3 <- qnorm(runif(n))
x2_3 <- qnorm(runif(n))
x3_3 <- qnorm(runif(n))

# generate x1=invnorm(uniform())+x1_1^2+abs(x2_1)^(0.5)+x3_1^3
# generate x2=invnorm(uniform())+exp(x1_2)+exp(1/x2_2)+(x3_2)^(2)
# generate x3=invnorm(uniform())+(x1_3)/x2_3+exp(x2_3)+(x3_3)^(-5)

x1 <- qnorm(runif(n)) + x1_1^2 + abs(x2_1)^(0.5) + x3_1^3
x2 <- qnorm(runif(n)) + exp(x1_2) + exp(1/x2_2) + (x3_2)^(2)
x3 <- qnorm(runif(n)) + (x1_3)/x2_3 + exp(x2_3) + (x3_3)^(-5)

# drop in 9
# drop in 13
# drop in 5
# generate y1=x1+x2+x3+e
# gen year=2000+_n-1
# tw (connected y1 year , xlabel(2000(5)2016))

year <- 2000:(2000 + n - 1)
data <- data.frame(
  year = year,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x1_1 = x1_1,
  x2_1 = x2_1,
  x3_1 = x3_1,
  x1_2 = x1_2,
  x2_2 = x2_2,
  x3_2 = x3_2,
  x1_3 = x1_3,
  x2_3 = x2_3,
  x3_3 = x3_3
  )
data$y1 <- data$x1 + data$x2 + data$x3 + e
data <- data %>% select(year, y1, everything())

##' @note It is not clear, why certain columns are dropped in the end. I assume
##' it is because of certain very high values. I thus do a similar thing.

data <- data %>% 
  filter(abs(y1) < 120) %>% 
  slice(1:16) %>% 
  mutate(year = 2000:2015)

##' @note I assume the tw command is the same as the twoway command which is 
##' used for plotting.

plot(data$year, data$y1, xlab = 'Year', ylab = 'GDP', type = 'l')

# generate y0=y1 if year<=2010
# replace y0=-100+x1+x2^0.7+ln(abs(x3))+e if year>2010
# gen y=y1
# gen id=1

n <- nrow(data)
data$y0 <- data$y1
data$y0[data$year > 2009] <- -100 + data$x1[data$year > 2009] + 
  pmax(0,data$x2[data$year > 2009])^0.7 +
  log(abs(data$x3[data$year > 2009])) + 
  rnorm(n-10)
data$y <- data$y1
data$id <- 1

# gen y_1=y1+rnormal(0,2) if year<=2010
# replace y_1=-100+x1+x2^0.7+ln(abs(x3))+rnormal(15,30) if year>2010
# gen y_2=y1+rnormal(0,3) if year<=2010
# replace y_2=-100+x1+x2^0.7+ln(abs(x3))+rnormal(10,15) if year>2010
# gen y_3=y1+rnormal(0,4) if year<=2010
# replace y_3=-100+x1+x2^0.7+ln(abs(x3))+rnormal(5,10) if year>2010

data$y_1 <- data$y1 + rnorm(n, 0, 2)
data$y_2 <- data$y1 + rnorm(n, 0, 3)
data$y_3 <- data$y1 + rnorm(n, 0, 4)

data$y_1[data$year > 2009] <- -100 + data$x1[data$year > 2009] + 
  Re(as.complex(data$x2[data$year > 2009])^0.7) +
  log(abs(data$x3[data$year > 2009])) + 
  rnorm(sum(data$year > 2009), 15, 30)

data$y_2[data$year > 2009] <- -100 + data$x1[data$year > 2009] + 
  Re(as.complex(data$x2[data$year > 2009])^0.7) +
  log(abs(data$x3[data$year > 2009])) + 
  rnorm(sum(data$year > 2009), 10, 15)

data$y_3[data$year > 2009] <- -100 + data$x1[data$year > 2009] + 
  Re(as.complex(data$x2[data$year > 2009])^0.7) +
  log(abs(data$x3[data$year > 2009])) + 
  rnorm(sum(data$year > 2009), 5, 10)

# drop if year==2007
# replace year=2000+_n-1

##' @note It is not clear, why certain columns are dropped in the end. I decided
##' not to do it, since it doesn't seem to make sense.


# gen TE=y1-y0
# order id year y y1 y0 TE

data$TE <- data$y1 - data$y0
data <- data[order(data$id, data$year, data$y, data$y1, data$y0, data$TE), ]

# sort year
# tw (connected y1 year , sort(year)) (connected y0 year , sort(year) lp(dash)) , ///
# xline(2009 , lp(solid)) scheme(s2mono) ///
# legend(label(1 "Factual") label(2 "Counterfactual"))

plot(x = data$year, y = data$y1, type = 'l', col = .c[3], 
     xlab = 'Year', ylab = 'GDP', ylim = c(-150,150))
lines(x = data$year, y = data$y0, type = 'l', col = .c[9])
abline(v = 2009, lty = 'dashed')
legend(x = 'topleft', legend = c('Factual','Counterfactual'), 
       col = .c[c(3,9)], 
       lty= rep(1, 2))

# line TE year , sort(year) xline(2009 , lp(solid)) ylabel(-150(50)200) scheme(s1mono) 

plot(x = data$year, y = data$TE, type = 'l', lty = 'dashed', col = .c[10],
     xlab = 'Year', ylab = 'GDP')

# tw (line y year) (line y_1 year) (line y_2 year) (line y_3 year) , xline(2009 , lp(solid)) scheme(s2mono) ///
# legend(label(1 "Treated") label(2 "Donor 1") label(3 "Donor 2") label(4 "Donor 3"))

plot(x = data$year, y = data$y, type = 'l', col = .c[3], 
     xlab = 'Year', ylab = 'GDP', ylim = c(-150,150))
lines(x = data$year, y = data$y_1, type = 'l', col = .c[9])
lines(x = data$year, y = data$y_2, type = 'l', col = .c[10])
lines(x = data$year, y = data$y_3, type = 'l', col = .c[5])
abline(v = 2009, lty = 'dashed')
legend(x = 'topleft', legend = c('Treated','Donor 1','Donor 2','Donor 3'), 
       col = .c[c(3,9,10,5)], 
       lty= rep(1, 4))

# preserve
# keep year y0 
# rename y0 _y_0_dgp 
# save DGP , replace
# restore

##' @note In Stata, the preserve command is used to save the current state of 
##' the dataset in memory, including the data, variable labels, value labels, 
##' and other dataset attributes. This is useful when you want to make changes 
##' to the dataset and then later return to the original state, via the restore
##' command.

DGP <- data %>% 
  select(year, y0) %>% 
  rename(`_y_0_gdp` = y0)

# gen id_1=2
# gen id_2=3
# gen id_3=4

data <- data %>% 
  mutate(id_1 = 2,
         id_2 = 3,
         id_3 = 4)

# order id year y x1 x2 x3 ///
# id_1 y_1 x1_1 x2_1 x3_1  ///
# id_2 y_2 x1_2 x2_2 x3_2  ///
# id_3 y_3 x1_3 x2_3 x3_3  
#
# keep id year y x1 x2 x3 ///
# id_1 y_1 x1_1 x2_1 x3_1  ///
# id_2 y_2 x1_2 x2_2 x3_2  ///
# id_3 y_3 x1_3 x2_3 x3_3  

data <- data %>% 
  arrange(id, year, y, x1, x2, x3, id_1, y_1, x1_1, x2_1, x3_1, id_2, y_2, x1_2, 
          x2_2, x3_2, id_3, y_3, x1_3, x2_3, x3_3) %>% 
  select(id, year, y, x1, x2, x3, id_1, y_1, x1_1, x2_1, x3_1, id_2, y_2, x1_2, 
         x2_2, x3_2, id_3, y_3, x1_3, x2_3, x3_3)

# save data , replace
# preserve
# keep id year y x1 x2 x3
# save data , replace
# restore

data.backup <- data
data <- data %>% select(id, year, y, x1, x2, x3)

# preserve
# keep id_1 year y_1 x1_1 x2_1 x3_1
# rename id_1 id
# rename y_1 y
# rename x1_1 x1
# rename x2_1 x2
# rename x3_1 x3
# save data1, replace
# restore

data1 <- data.backup %>% 
  select(id_1, year, y_1, x1_1, x2_1, x3_1) %>% 
  rename(y = y_1,
         id = id_1,
         x1 = x1_1,
         x2 = x2_1,
         x3 = x3_1)

# preserve
# keep id_2 year y_2 x1_2 x2_2 x3_2
# rename id_2 id
# rename y_2 y
# rename x1_2 x1
# rename x2_2 x2
# rename x3_2 x3
# save data2 , replace
# restore

data2 <- data.backup %>% 
  select(id_2, year, y_2, x1_2, x2_2, x3_2) %>% 
  rename(y = y_2,
         id = id_2,
         x1 = x1_2,
         x2 = x2_2,
         x3 = x3_2)

# preserve
# keep id_3 year y_3 x1_3 x2_3 x3_3
# rename id_3 id
# rename y_3 y
# rename x1_3 x1
# rename x2_3 x2
# rename x3_3 x3
# save data3 , replace
# restore

data3 <- data.backup %>% 
  select(id_3, year, y_3, x1_3, x2_3, x3_3) %>% 
  rename(y = y_3,
         id = id_3,
         x1 = x1_3,
         x2 = x2_3,
         x3 = x3_3)

# use data , clear
# append using data1 data2 data3

data <- rbind(data, data1, data2, data3)

# tsset id year
# global xvars "x1 x2 x3"

##' @note In Stata, the `tsset` command is used to declare the dataset as a 
##' time-series dataset. This command sets the dataset's time variable and, if 
##' applicable, any panel variable, which is crucial for Stata to understand the 
##' time structure of the data. Time-series datasets are essential for many 
##' time-series analyses, including panel data analyses where observations are 
##' repeated for different units over time.

##' @note In Stata, the `global` command is used to create or modify global 
##' macros, which are variables that store strings of text. The purpose of using 
##' global macros is to define variables that can be accessed and used across 
##' different Stata commands or in subsequent analyses.

xvars <- c('x1', 'x2', 'x3')

# * PARAMETRIC 
# synth y  $xvars , trunit(1) trperiod(2009) figure keep(SYNTH_data , replace)

dataprep.out <- dataprep(foo = data, predictors = xvars, predictors.op = "mean",
                         time.predictors.prior = 2000:2009, dependent = "y",
                         unit.variable = "id", time.variable = "year", 
                         treatment.identifier = 1, controls.identifier = c(2:4),
                         time.optimize.ssr = 2000:2009, time.plot = 2000:2015)

synth.out <- synth(dataprep.out)

# * NON-PARAMETRIC 
# label define LAB 1 "ITALY" 2 "GERMANY" 3 "FRANCE" 4 "UK"  , replace
# label val id LAB
# npsynth y $xvars ,  npscv panel_var(id) time_var(year) trperiod(2009)  trunit(1) bandw(1.55) kern(triangular) gr1 gr2 gr3  ///
# save_gr1(gr1) save_gr2(gr2) save_gr3(gr3) gr_y_name("graph") gr_tick(5) save_res(NPSYNTH_data) 
# return list

data <- data %>% 
  mutate(id = case_when(
    id == 1 ~ 'ITALY',
    id == 2 ~ 'GERMANY',
    id == 3 ~ 'FRANCE',
    id == 4 ~ 'UK'
  ))

##' From `npsynth` documentation for stata: 
##' https://www.stata.com/meeting/uk17/slides/uk17_Cerulli.pdf

possible.bw <- c(seq(0.01,1.99,by=0.01),seq(2,5,by=0.1),6:50)
mses <- matrix(NA, length(possible.bw), 3)
colnames(mses) <- c('bw','mse','mspe')

plot = F
if (plot == T) par(mfrow=c(3,2))

# Get a feeling for the function and find optimal bandwidth
for (bw in possible.bw){
  
  x <- npsynth(X = data,
               k = 'y',
               id = 'id',
               id.interest = 'ITALY',
               periods = 'year',
               intervention = 2009,
               controls = c('x1','x2','x3'),
               h = bw)
  
  if (bw == possible.bw[1]) i = 1
  mses[i, ] <- c(bw, x$mse, x$mspe); i <- i + 1
  
  # Plotting
  if (plot == T){
    plot(2000:2015, preds$synthetic, type = 'l', col = .c[2], ylim = c(-150,150),
         xlab = 'Year', ylab = 'GDP', main = paste('BW:',bw))
    lines(2000:2015, preds$dgp, type = 'l', col=.c[4])
    lines(2000:2015, x$predictions, type = 'l', col=.c[8])
    abline(v = 2009, lty = 2)
    legend(x = 'topleft', legend = c('Parametric','DGP','Non-Parametric'), 
           col = .c[c(2,4,8)], 
           lty= rep(1, 3))
  }
  
  if (bw == tail(possible.bw, 1)){
    cat(paste('h.opt MSPE:', mses[which.min(mses[,'mspe'])]))
    cat(paste('\nh.opt MSE: ', mses[which.min(mses[,'mse'])]))
  } 
}

# Use optimal bandwidth
npsynth.out <- npsynth(X = data,
               k = 'y',
               id = 'id',
               id.interest = 'ITALY',
               periods = 'year',
               intervention = 2009,
               controls = c('x1','x2','x3'),
               h = 14,
               dgp = DGP$`_y_0_gdp`)

# preserve
# use NPSYNTH_data , clear
# keep year _Y0_
# rename _Y0_ _y_0_npsynth
# save NPSYNTH_data , replace
# restore
#
# preserve
# use SYNTH_data , clear
# keep _Y_synthetic _time
# rename _time year
# rename _Y_synthetic _y_0_synth
# save SYNTH_data , replace
# restore
#
# use DGP , clear
# merge 1:1 year using SYNTH_data
# cap drop _merge
# merge 1:1 year using NPSYNTH_data

DGP <- DGP %>% 
  mutate(`_y_0_synth` = dataprep.out$Y0plot %*% synth.out$solution.w,
         `_y_0_npsynth` = npsynth.out$predictions)
  

# tw (line _y_0_dgp year) (line _y_0_synth year) (line _y_0_npsynth year) if year>=2009 , scheme(s2mono) ///
# legend(label(1 "DGP") label(2 "SYNTH") label(3 "NPSYNTH")) xlabel(2009(1)2015)

plot(2000:2015, DGP$`_y_0_gdp`, type = 'l', col = .c[2], ylim = c(-150,150),
     xlab = 'Year', ylab = 'GDP')
lines(2000:2015, DGP$`_y_0_synth`, type = 'l', col=.c[4])
lines(2000:2015, DGP$`_y_0_npsynth`, type = 'l', col=.c[8])
abline(v = 2009, lty = 2)
legend(x = 'topright', legend = c('DGP','Parametric','Non-Parametric'), 
       col = .c[c(2,4,8)], 
       lty= rep(1, 3))

# gen DEV_gdp_synth=(_y_0_dgp - _y_0_synth)^2
# qui sum DEV_gdp_synth
# global RMSPE_gdp_synth=sqrt(r(mean))
#
# gen DEV_gdp_npsynth=(_y_0_dgp - _y_0_npsynth)^2
# qui sum DEV_gdp_npsynth
# global RMSPE_gdp_npsynth=sqrt(r(mean))
#
# di  $RMSPE_gdp_synth
# di  $RMSPE_gdp_npsynth

##' @note In Stata, the `qui` command is a shorthand abbreviation for quietly. 
##' The purpose of the quietly command is to suppress the output of a command 
##' while still allowing the command to execute. It's often used to perform 
##' calculations or operations without cluttering the Stata results window with 
##' unnecessary information.

##' @note In Stata, the `di` command is short for display. It is used to display 
##' output, messages, or the values of variables in the Results window.

RMSPE_gdp_synt <- sqrt(mean((DGP$`_y_0_gdp` - DGP$`_y_0_synth`)^2))
RMSPE_gdp_npsynt <- sqrt(mean((DGP$`_y_0_gdp` - DGP$`_y_0_npsynth`)^2))

cat(paste('RMSE COMPARISON:\nSynth:  ',RMSPE_gdp_synt,
          '\nNpsynth:',RMSPE_gdp_npsynt,'\n'))


# R Code ------------------------------------------------------------------

##' @note This code replicates the code from the previous section, just it is
##' written in a cleaner way for R.

set.seed(101)

# Parameters
N <- 50
n <- 16
n_x <- 3
p_int <- 2009

# Generate data
data <- as.data.frame(matrix(rnorm(N * n_x * 3), N, n_x * 3))
names(data) <- paste0(paste0('x',rep(1:3,each=3)),'_',rep(1:3,3))
data <- data %>% 
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

# Set aside the DGP
DGP <- data %>% 
  select(year, y, y0) %>% 
  rename(factual = y, counterfactual = y0)

# Re-format the data
for (i in 1:3){
  .d <- data %>% 
    select(all_of(c(
      paste0('id',i), 'year',
      paste0('y',i), paste0('x',i,'_',1:3))
    ))
  names(.d) <- c('id','year','y','x1','x2','x3')
  assign(paste0('data',i),.d)
}

data <- data %>% 
  select(id, year, y, x1, x2, x3) %>% 
  rbind(data1, data2, data3)

# See time series plots
pdf(paste0(here(),'/../Paper/Figures/Simulation/dgp_realization.pdf'),
    height = 5, width = 8) 

plot(2000:(1999+n), data %>% filter(id == 1) %>% pull(y), type = 'l', col = .c2[1], 
     ylim = c(-150,150), xlab = 'Year', ylab = 'y')
lines(2000:(1999+n), DGP$counterfactual, type = 'l', col=.c2[2])
lines(2000:(1999+n), data %>% filter(id == 2) %>% pull(y), type = 'l', col=.c2[3])
lines(2000:(1999+n), data %>% filter(id == 3) %>% pull(y), type = 'l', col=.c2[4])
lines(2000:(1999+n), data %>% filter(id == 4) %>% pull(y), type = 'l', col=.c2[5])
abline(v = p_int, lty = 2)
legend(x = 'topright', legend = c('Treatment','Counterfactual','Control 1',
                                  'Control 2','Control 3'), 
       col = .c2[1:5], 
       lty= rep(1, 5))

dev.off()

# Compute synthetic control groups
dataprep.out <- dataprep(foo = data, predictors = c('x1', 'x2', 'x3'), 
                         predictors.op = "mean",
                         time.predictors.prior = 2000:2009, dependent = "y",
                         unit.variable = "id", time.variable = "year", 
                         treatment.identifier = 1, controls.identifier = c(2:4),
                         time.optimize.ssr = 2000:2009, time.plot = 2000:2015)

synth.out <- synth(dataprep.out)

possible.bw <- c(seq(0.01,1.99,by=0.01),seq(2,5,by=0.1),6:50)
mses <- matrix(NA, length(possible.bw), 3)
colnames(mses) <- c('bw','mse','mspe')

plot = F
if (plot == T) par(mfrow=c(3,2))

# Get a feeling for the function and find optimal bandwidth
for (bw in possible.bw){
  
  x <- uniKernel(X = data,
                 k = 'y',
                 id = 'id',
                 id.interest = 1,
                 periods = 'year',
                 intervention = 2009,
                 controls = c('x1','x2','x3'),
                 kernel = 'Triangular',
                 method = 'Aggregate',
                 h = bw,
                 true_cf = DGP$counterfactual)
  
  if (bw == possible.bw[1]) i = 1
  mses[i, ] <- c(bw, x$mse, x$mspe); i <- i + 1
  
  # Plotting
  if (plot == T){
    plot(2000:2015, preds$synthetic, type = 'l', col = .c[2], ylim = c(-150,150),
         xlab = 'Year', ylab = 'GDP', main = paste('BW:',bw))
    lines(2000:2015, preds$dgp, type = 'l', col=.c[4])
    lines(2000:2015, x$predictions, type = 'l', col=.c[8])
    abline(v = 2009, lty = 2)
    legend(x = 'topleft', legend = c('Parametric','DGP','Non-Parametric'), 
           col = .c[c(2,4,8)], 
           lty= rep(1, 3))
  }
  
  if (bw == tail(possible.bw, 1)){
    cat(paste('h.opt MSPE:', mses[which.min(mses[,'mspe'])]))
    cat(paste('\nh.opt MSE: ', mses[which.min(mses[,'mse'])]))
  } 
}

# Use optimal bandwidth
npsynth.out <- uniKernel(X = data,
                         k = 'y',
                         id = 'id',
                         id.interest = 1,
                         periods = 'year',
                         intervention = 2009,
                         controls = c('x1','x2','x3'),
                         kernel = 'Triangular',
                         method = 'Aggregate',
                         h = 3.7,
                         true_cf = DGP$counterfactual)

# Let the function select the optimal bandwidth by not defining it
npsynth.new <- uniKernel(X = data,
                         k = 'y',
                         id = 'id',
                         id.interest = 1,
                         periods = 'year',
                         intervention = 2009,
                         controls = c('x1','x2','x3'),
                         kernel = 'Triangular',
                         method = 'Aggregate',
                         true_cf = DGP$counterfactual)

# Append DGP data frame
DGP <- DGP %>% 
  mutate(synth = as.vector(dataprep.out$Y0plot %*% synth.out$solution.w),
         npsynth = npsynth.out$predictions)

# Plots
plot(2000:2015, DGP$counterfactual, type = 'l', col = .c[2], ylim = c(-150,150),
     xlab = 'Year', ylab = 'GDP')
lines(2000:2015, DGP$synth, type = 'l', col=.c[4])
lines(2000:2015, DGP$npsynth, type = 'l', col=.c[8])
abline(v = 2009, lty = 2)
legend(x = 'topright', legend = c('DGP','Parametric','Non-Parametric'), 
       col = .c[c(2,4,8)], 
       lty= rep(1, 3))

# Performance measures
sqrt(mean((DGP$counterfactual - DGP$synth)^2))
sqrt(mean((DGP$counterfactual - DGP$npsynth)^2))

# MSPE
idx <- 2010:2015 - 1999
sqrt(mean((DGP$counterfactual[idx] - DGP$synth[idx])^2))
sqrt(mean((DGP$counterfactual[idx] - DGP$npsynth[idx])^2))


# Save data ---------------------------------------------------------------

saveRDS(data, paste0(here(), '/npsynth/sim_data.rds'))

