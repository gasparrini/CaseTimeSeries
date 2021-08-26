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
# SIMULATED SCENARIO: COMPLEX LAG STRUCTURE
################################################################################

# PACKAGES
library(splines) ; library(pbs)
library(gnm); library(survival)
library(data.table) ; library(scales)
library(foreach) ; library(doParallel)
library(dlnm)

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

# TRUE EFFECT OF CONTINUOUS EXPOSURE CUMULATED OVER LAG 0-10
# NB: SUM OF AN AVERAGE EFFECT ALONG THE 11 LAG TIMES
beta <- log(1.0025) * 11

################################################################################
# EXAMPLE

# SIMULATE AN (UNOBSERVED) BASELINE CONFOUNDER
data[, zfix := runif(1,0,100), by=id]

# SIMULATE SEASONAL CONTINUOUS EXPOSURE CORRELATED WITH BASELINE CONFOUNDER
fexp <- function(date) exp(cos(yday(date)*2*pi/365) + 1 +
    as.numeric(arima.sim(list(ar=0.5), length(date), sd=0.5)))
data[, x := fexp(date) * zfix/50, by=id]
cor(data[,.(x,zfix)])
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

# DEFINE EXPOSURE HISTORY
exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])

# DEFINE THE LAG WEIGHT FUNCTION AND LAGGED EFFECT MATRIX
fwlag <- function(lag) 5*dnorm(lag,2,2)
plot(0:100/10, fwlag(0:100/10), type="l", ylab="Weight", main="Lag structure")
betalag <- beta * fwlag(0:10) / sum(fwlag(0:10))
matbetalag <- matrix(betalag, nrow(data), 11, byrow=T)

# COMPUTE RISK ASSOCIATED TO EXPOSURE, TREND, AND TIME-VAR CONFOUNDER
data$rrexp <- exp(rowSums(exphist*matbetalag))
data[, rrtrend := frrtrend(date, tpar+rnorm(4,0,rep(0.05,4))), by=id]
data$rrzvar <- exp(0.01*data$zvar)

# SIMULATE THE OUTCOME (SUBJECT-SPECIFIC RISK DEPENDENT ON BASELINE CONFOUNDER)
data[, y := rmultinom(1, round(first(zfix)/5+1), prob=rrexp*rrtrend*rrzvar),
  by=id]

# DEFINE CYCLIC SPLINES FOR SEASONALITY
splseas <- pbs(yday(data$date), df=6)

# DEFINE SUBJECT-MONTH STRATA
data$stratum <- with(data, factor(paste(id, year(date), month(date), sep="-")))

# DEFINE THE CROSSBASIS FOR EXPOSURE
cb <- crossbasis(exphist, lag=10, argvar=list(fun="lin"),
  arglag=list(fun="ns", knots=1:3*2.5))
#summary(cb)

# RUN THE MODEL (WITHOUT BASELINE CONFOUNDER)
mod <- gnm(y ~ cb + date + splseas + zvar, eliminate=stratum, data=data,
  family="poisson")

# PREDICTIONS AND PLOTTED LAG-RESPONSE
cp <- crosspred(cb, mod, at=1, cen=0, bylag=0.1)
plot(cp, var=1, ylab="RR", main="Lag-response association")
#lines(0:100/10, exp(beta*fwlag(0:100/10)/sum(fwlag(0:10))), lty=2)
#legend("top", c("True", "Estimated"), lty=2:1, bty="n")

################################################################################
# SIMULATIONS

# SET THE PARALLELIZATION
ncores <- detectCores() 
cl <- makeCluster(max(1,ncores-2))
registerDoParallel(cl)

# ITERATIONS USING FOREACH
pkg <- c("dlnm","splines","gnm","survival","data.table")
eff <- foreach(i=seq(nsim), .packages=pkg, .combine=rbind) %dopar% {
  
  # SET SEED
  set.seed(13041975+i)
  
  # SIMULATE EXPOSURE, RISKS, OUTCOME (AS ABOVE)
  data[, zfix := runif(1,0,100), by=id]
  data[, x := fexp(date) * zfix/50, by=id]
  exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
  data$rrexp <- exp(rowSums(exphist*matbetalag))
  tpar <- c(runif(2,-0.2,0.2),runif(2,-4,4))
  data[, rrtrend := frrtrend(date, tpar+rnorm(4,0,rep(0.05,4))), by=id]
  data$zvar <- data$x + rnorm(nrow(data),0,3)
  data$rrzvar <- exp(0.01*data$zvar)
  data[, y := rmultinom(1, round(first(zfix)/5+1), prob=rrexp*rrtrend*rrzvar),
    by=id]
  
  # DEFINE THE CROSSBASIS
  cb <- crossbasis(exphist, lag=10, argvar=list(fun="lin"),
    arglag=list(fun="ns", knots=1:3*2.5))
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  mod <- gnm(y ~ cb + date + splseas + zvar, eliminate=stratum, data=data,
    family="poisson")
  cp <- crosspred(cb, mod, at=1, cen=0, bylag=0.1)
  
  # STORE THE RESULTS (INCLUDING LAG-SPECIFIC)
  cbind(est = cp$allfit, cov = log(cp$allRRlow)<beta & log(cp$allRRhigh)>beta,
    cp$matfit)
}

# STOP CLUSTER
stopCluster(cl)

# EXTRACT THE LAG ESTIMATES (MEAN AND SAMPLE) AND THEN ERASE
efflagmean <- colMeans(eff[,-c(1:2)])
efflagsample <- eff[seq(min(nsim,100)),-c(1:2)]
eff <- eff[,1:2]

################################################################################
# RESULTS

# RELATIVE BIAS (%) 
(bias <- (mean(eff[,1])-beta)/beta * 100)

# COVERAGE
(cov <- mean(eff[,2]))

# RELATIVE RMSE (%) 
(rmse <- sqrt(mean((eff[,1]-beta)^2))/beta * 100)
