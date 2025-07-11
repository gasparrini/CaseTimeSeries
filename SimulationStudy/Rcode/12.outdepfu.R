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
# SIMULATED SCENARIO: OUTCOME-DEPENDENT FOLLOW-UP PERIOD
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

# SIMULATE EXPOSURE EPISODES (10% OF TIMES)
nexp <- round(length(date)*0.1)
data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]

# DEFINE EXPOSURE HISTORY AND THE CUMULATED EXPOSURE OVER LAG TIMES 0-10
exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
data$cumexp <- rowSums(exphist)

# COMPUTE RISK ASSOCIATED TO EXPOSURE
data$rrexp <- exp(data$cumexp*beta)

# SIMULATE THE OUTCOME (5-20 EVENTS)
data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp), by=id]

# SIMULATE CENSORING AFTER OUTCOME
data[, cens := rbinom(length(y), 1, 0.2) & y>0, by=id]
data[, keep := seq(length(y)) < min(c(length(y), which(cens))) + 1, by=id]

# PLOT
plot(date, date, ylim=c(0.5,10+0.5), yaxt="n", ylab="", xlab="Follow-up",
  frame.plot=F, main="Outcome-dependent censoring")
axis(1, at=range(date)+c(-500,500), labels=F, lwd.ticks=0)
axis(2, at=rev(seq(10)), labels=paste("Sub",seq(10)), cex.axis=0.8, lwd=0, las=1)
for(i in seq(10)) {
  sub <- subset(data, id==rev(seq(10))[i])
  lines(date, rep(i, length(date)), col=grey(0.9))
  lines(sub$date[sub$keep], rep(i, sum(sub$keep)))
  for(t in seq(length(sub$date))[sub$y>0]) {
    points(rep(sub$date[t], sub$y[t]), i+(seq(sub$y[t])-mean(seq(sub$y[t])))/4,
      pch=21, cex=1.5, bg=c(grey(0.9),2)[sub$keep[t]+1],
      col=c(grey(0.9),1)[sub$keep[t]+1])
  }
}

# RUN THE MODEL
mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson", subset=keep)

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
  data[, x := seq(length(date)) %in% sample(length(date), nexp) + 0, by=id]
  exphist <- as.matrix(data[, shift(x, 0:10, fill=0), by=id][,2:12])
  data$cumexp <- rowSums(exphist)
  data$rrexp <- exp(data$cumexp*beta)
  data[, y := rmultinom(1, sample(5:20, 1), prob=rrexp), by=id]
  
  # SIMULATE CENSORING AFTER OUTCOME
  data[, cens := rbinom(length(y), 1, 0.2) & y>0, by=id]
  data[, keep := seq(length(y)) < min(c(length(y), which(cens))) + 1, by=id]
  
  # RUN THE MODEL AND EXTRACT COEF-VCOV-CI
  mod <- gnm(y ~ cumexp, eliminate=id, data=data, family="poisson", subset=keep)
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
