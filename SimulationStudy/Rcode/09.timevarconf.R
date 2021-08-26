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
# SIMULATED SCENARIO: TIME-VARYING CONFOUNDER
################################################################################

# PACKAGES
library(splines) ; library(pbs)
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

# SIMULATE THE TREND (COMMON BASELINE + SUBJECT-SPECIFIC DEVIATIONS)
tpar <- c(-0.15,0.06,-3.5,-1.2)
frrtrend <- function(date, tpar) 
  exp(tpar[1]*cos(yday(date)*tpar[3]*pi/365) + 
      tpar[2]*cos(yday(date)*tpar[4]*pi/365))
plot(date, frrtrend(date,tpar), type="l", ylim=c(0.5,1.5), col=1, lwd=2, 
  ylab="RR", main="Seasonal trend")
for(i in 1:7) lines(date, frrtrend(date, tpar+rnorm(4,0,rep(0.05,4))),
  type="l", lty=2, col=i)
abline(h=1)

# SIMULATE A TIME-VARYING CONFOUNDER WITH SEASONAL DISTRIBUTION
data$zvar <- data$x + rnorm(nrow(data),0,3)
cor(data[,.(x,zvar)])
plot(zvar~date, data[id%in%1:10], ylab="Confounder",
  main="Confounder distribution")
mtext("First 10 subjects")

# DEFINE EXPOSURE HISTORY AND THE CUMULATED EXPOSURE OVER LAG TIMES 0-10
exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
data$cumexp <- rowSums(exphist)

# COMPUTE RISK ASSOCIATED TO EXPOSURE, TREND, AND TIME-VAR CONFOUNDER
data$rrexp <- exp(data$cumexp*beta)
data[, rrtrend := frrtrend(date, tpar+rnorm(4,0,rep(0.05,4))), by=id]
data$rrzvar <- exp(0.01*data$zvar)

# SIMULATE THE OUTCOME (5-20 EVENTS)
data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp*rrtrend*rrzvar), by=id]

# DEFINE CYCLIC SPLINES FOR SEASONALITY
splseas <- pbs(yday(data$date), df=6)

# DEFINE SUBJECT-MONTH STRATA
data$stratum <- with(data, factor(paste(id, year(date), month(date), sep="-")))

# RUN THE MODEL
mod <- gnm(y ~ cumexp + date + splseas + zvar, eliminate=stratum, data=data,
  family="poisson")

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
  data$rrexp <- exp(data$cumexp*beta)
  tpar <- c(runif(2,-0.2,0.2),runif(2,-4,4))
  data[, rrtrend := frrtrend(date, tpar+rnorm(4,0,rep(0.05,4))), by=id]
  data$zvar <- data$x + rnorm(nrow(data),0,3)
  data$rrzvar <- exp(0.01*data$zvar)
  data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp*rrtrend*rrzvar), by=id]
  
  # RUN THE MODEL AND EXTRACT COEF-VCOV
  mod <- gnm(y ~ cumexp + date + splseas + zvar, eliminate=stratum, data=data,
    family="poisson")
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
