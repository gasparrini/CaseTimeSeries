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
# SIMULATED SCENARIO: RARE OUTCOME/EXPOSURE
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

# TRUE EFFECT OF EXPOSURE EPISODE
beta <- log(1.15)

################################################################################
# EXAMPLE

# SIMULATE EXPOSURE EPISODES (RARE: 1-5 EPISODES)
nexp <- sample(1:5, 1)
data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]

# DEFINE EXPOSURE HISTORY AND THE CUMULATED EXPOSURE OVER LAG TIMES 0-10
exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
data$cumexp <- rowSums(exphist)

# COMPUTE RISK ASSOCIATED TO EXPOSURE
data$rrexp <- exp(data$cumexp*beta)

# SIMULATE THE OUTCOME (RARE: 1-3 EVENTS)
data[, y := rmultinom(1, sample(1:3, 1), prob=rrexp), by=id]

# RUN THE MODEL
mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson")

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
  nexp <- sample(1:5, 1)
  data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]
  exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
  data$cumexp <- rowSums(exphist)
  data$rrexp <- exp(data$cumexp*beta)
  data[, y := rmultinom(1, sample(1:3, 1), prob=rrexp), by=id]
  
  # RUN THE MODEL AND EXTRACT COEF-VCOV-CI
  mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson")
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
