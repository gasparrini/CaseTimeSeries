################################################################################
# Updated version of the R code for the simulation study in:
#
#  "The case time series analysis"
#  Antonio Gasparrini
#  http://www.ag-myresearch.com/2021_gasparrini_epidemiol.html
#
# * an updated version of this code is available at:
#   https://github.com/gasparrini/CaseTimeSeries
################################################################################

################################################################################
# SIMULATED SCENARIO: BINARY OUTCOME
################################################################################

# PACKAGES
library(splines)
library(gnm); library(survival)
library(data.table) ; library(scales)
library(foreach) ; library(doParallel)

# SET SEED
set.seed(13041975)

# SET SIMULATION SETTINGS
nsim <- 50000#

# SET PARAMETERS: NUMBER OF SUBJECTS, FOLLOW-UP, TIME VARIABLES
n <- 500
dstart <- as.Date("2019-01-01")
dend <- as.Date("2019-12-31")
date <- seq(dstart, dend, by=1)

################################################################################
# PREPARE THE DATA

# GENERATE DATASET
data <- data.table(id = factor(rep(seq(n), each=length(date))), date = date)

# TRUE EFFECT OF CONTINUOUS EXPOSURE 
beta <- log(1.0025)

################################################################################
# EXAMPLE

# SIMULATE SEASONAL CONTINUOUS EXPOSURE
fexp <- function(date) exp(cos(yday(date)*2*pi/365) + 1 +
    as.numeric(arima.sim(list(ar=0.5), length(date), sd=0.5)))
data[, x := fexp(date), by=id]
plot(fexp(date)~date, data[id%in%1:10], type="p", ylab="Exposure", xlab="Date",
  main="Exposure distribution")
mtext("First 10 subjects")

# DEFINE EXPOSURE HISTORY AND THE CUMULATED EXPOSURE OVER LAG TIMES 0-10
exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
data$cumexp <- rowSums(exphist)

# COMPUTE BASELINE PROBABILITY AND RISK ASSOCIATED TO EXPOSURE
pbase <- 0.10
data$orexp <- exp(data$cumexp*beta)

# SIMULATE THE OUTCOME (BINARY)
data[, y := rbinom(length(date), 1, plogis(qlogis(pbase) + log(orexp))), by=id]

# RUN THE MODEL
mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="binomial")

# ESTIMATED EFFECT
beta ; coef(mod)[[1]]

################################################################################
# SIMULATIONS

# SET THE PARALLELIZATION
ncores <- detectCores()
cl <- makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# ITERATIONS USING FOREACH
pkg <- c("splines","gnm","survival","data.table")
eff <- foreach(i=seq(nsim), .packages=pkg, .combine=rbind) %dopar% {
  
  # SET SEED
  set.seed(13041975+i)
  
  # SIMULATE EXPOSURE, RISKS, OUTCOME (AS ABOVE)
  data[, x := fexp(date), by=id]
  exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
  data$cumexp <- rowSums(exphist)
  data$orexp <- exp(data$cumexp*beta)
  data[, y := rbinom(length(date), 1, plogis(qlogis(pbase) + log(orexp))), by=id]
  
  # RUN THE MODEL AND EXTRACT COEF-VCOV
  mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="binomial")
  coef <- coef(mod)[[1]]
  vcov <- as.matrix(vcov(mod))[1,1]
  cilow <- (coef-qnorm(0.975)*sqrt(vcov)) 
  cihigh <- (coef+qnorm(0.975)*sqrt(vcov))
  
  # STORE THE RESULTS
  cbind(est = coef, cov = cilow < beta & cihigh > beta)
}

# STOP CLUSTER
stopCluster(cl)

################################################################################
# RESULTS

# RELATIVE BIAS (%) 
(bias <- (mean(eff[,1])-beta)/beta * 100)

# COVERAGE
(cov <- mean(eff[,2]))

# RELATIVE RMSE (%) 
(rmse <- sqrt(mean((eff[,1]-beta)^2))/beta * 100)
