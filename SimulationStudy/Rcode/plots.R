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
# PRODUCE THE PLOTS
################################################################################

# SAVE AS 5X10 OR 5X5 PDF

################################################################################
# eFIGURE 1

par(las=1, mar=c(5,1,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

load("results/01.basic.RData")

plot(1:10, xlim=c(0.05,0.25), ylim=c(-0.5,0.5), type="n", yaxt="n", ylab="",
  xlab=expression(beta))
points(eff[,1], rnorm(nsim, 0, 0.07), col=grey(0.7))
abline(v=mean(eff[,1]), col=1, lwd=2)
abline(v=beta, col=2, lty=2, lwd=2)
title("Scenario 1: Basic")
mtext("Simulation estimates")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.7)), 
  lty=c(1,2,NA), lwd=2, pch=c(NA,NA,1), inset=0.02, y.intersp=1.3, cex=0.8)

par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))

load("results/10.complexlag.RData")

matplot(0:100, t(exp(efflagsample)), type="l", lty=1,
  col=grey(0.9), ylab="RR", xlab="Lag (days)")
lines(0:100, exp(beta*fwlag(0:100/10)/sum(fwlag(0:10))), col=1, lwd=2)
lines(0:100, exp(efflagmean), col=2, lty=2, lwd=2)
abline(h=1)
title("Scenario 10: Complex lag structure")
mtext("Lag-response association")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.9)), 
  lty=c(1,2,1), lwd=2, inset=0.02, y.intersp=1.3, cex=0.8)

################################################################################
# FIGURE S1

par(las=1, mar=c(5,1,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

load("results/01.basic.RData")

plot(1:10, xlim=c(0.05,0.25), ylim=c(-0.5,0.5), type="n", yaxt="n", ylab="",
  xlab=expression(beta))
points(eff[,1], rnorm(nsim, 0, 0.07), col=grey(0.7))
abline(v=mean(eff[,1]), col=1, lwd=2)
abline(v=beta, col=2, lty=2, lwd=2)
title("Scenario 1: Basic")
mtext("Simulation estimates")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.7)), 
  lty=c(1,2,NA), lwd=2, pch=c(NA,NA,1), inset=0.02, y.intersp=1.3, cex=0.8)

load("results/02.rare.RData")

plot(1:10, xlim=c(0.05,0.25), ylim=c(-0.5,0.5), type="n", yaxt="n", ylab="",
  xlab=expression(beta))
points(eff[,1], rnorm(nsim, 0, 0.07), col=grey(0.7))
abline(v=mean(eff[,1]), col=1, lwd=2)
abline(v=beta, col=2, lty=2, lwd=2)
title("Scenario 2: Rare outcome/exposure")
mtext("Simulation estimates")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.7)), 
  lty=c(1,2,NA), lwd=2, pch=c(NA,NA,1), inset=0.02, y.intersp=1.3, cex=0.8)

################################################################################
# FIGURE S2

# RUN FIRST PART OF 03.expcont.R, BFORE THE MODEL FITTING
lines <- scan("03.expcont.R", what="character", sep="\n")
linescoll <- paste(lines[seq(which(lines == "# RUN THE MODEL"))], collapse='\n')
source(textConnection(linescoll))

par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

plot(fexp(date)~date, data[id%in%1:10], type="p", ylab="Exposure", xlab="Date",
  main="Scenario 3: Continuous exposure")
mtext("Exposure distribution: first 10 subjects")

load("results/03.expcont.RData")

matplot(0:50, exp(sapply(eff[1:100,1]*11,"*",0:50)), type="l", lty=1,
  col=grey(0.9), ylab="RR", xlab="Lag (days)")
lines(0:50, exp(log(1.0025)*11*0:50), col=1, lwd=2)
lines(0:50, exp(mean(eff[,1])*11*0:50), col=2, lty=2, lwd=2)
abline(h=1)
title("Scenario 3: Continuous exposure")
mtext("Cumulative exposure-response")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.9)), 
  lty=c(1,2,1), lwd=2, inset=0.02, y.intersp=1.3, cex=0.8)

################################################################################
# FIGURE S3

# RUN FIRST PART OF 05.outcont.R, BFORE THE MODEL FITTING
lines <- scan("05.outcont.R",what="character",sep="\n")
linescoll <- paste(lines[seq(which(lines == "# RUN THE MODEL"))], collapse='\n')
source(textConnection(linescoll))
  
par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

plot(y~date, data[id%in%1:5], ylab="Outcome", xlab="Date",
  col=alpha(data[id%in%1:5]$id, 0.4))
abline(h=data[id%in%1:5,.(y=mean(y)),by=id]$y, col=1:5, lwd=2)
title("Scenario 5: Continuous outcome")
mtext("Outcome distribution: first 5 subjects")

load("results/05.outcont.RData")

matplot(0:50, sapply(eff[1:100,1]*11,"*",0:50), type="l", lty=1,
  col=grey(0.9), ylab="Effect", xlab="Exposure")
lines(0:50, beta*11*0:50, col=1, lwd=2)
lines(0:50, mean(eff[,1])*11*0:50, col=2, lty=2, lwd=2)
abline(h=0)
title("Scenario 5: Continuous outcome")
mtext("Cumulative exposure-response")
legend("topleft",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.9)), 
  lty=c(1,2,1), lwd=2, inset=0.02, y.intersp=1.3, cex=0.8)

################################################################################
# FIGURE S4

# RUN FIRST PART OF 06.trendcomm.R, BFORE THE MODEL FITTING
lines <- scan("06.trendcomm.R",what="character",sep="\n")
linescoll <- paste(lines[seq(which(lines == "# RUN THE MODEL"))], collapse='\n')
source(textConnection(linescoll))

par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

set.seed(13041975)
x <- sapply(seq(1:7), function(x) frrtrend(date,c(runif(2,-0.2,0.2),runif(2,-4,4))))
plot(date, rep(1, length(date)), type="n", ylim=range(x)+c(-0.05,0.05),
  ylab="RR", xlab="Date")
matlines(date, x, type="l", lty=1, col=1:7)
abline(h=1)
title("Scenario 6: Common trend")
mtext("Simulated seasonal trend")

set.seed(13041975)
tpar <- c(-0.15,0.06,-3.5,-1.2)
x <- sapply(seq(1:7), function(x) frrtrend(date,tpar+rnorm(4,0,rep(0.05,4))))
plot(date, frrtrend(date,tpar), type="l", ylim=range(x)+c(-0.05,0.05), col=1,
  lwd=2, ylab="RR", xlab="Date")
matlines(date, x, type="l", lty=2, col=1:7)
abline(h=1)
title("Scenario 7: Subject-specific trends")
mtext("Simulated seasonal trend")

################################################################################
# FIGURE S5

par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))
layout(1)

load("results/10.complexlag.RData")

matplot(0:100, t(exp(efflagsample)), type="l", lty=1,
  col=grey(0.9), ylab="RR", xlab="Lag (days)")
lines(0:100, exp(beta*fwlag(0:100/10)/sum(fwlag(0:10))), col=1, lwd=2)
lines(0:100, exp(efflagmean), col=2, lty=2, lwd=2)
abline(h=1)
title("Scenario 10: Complex lag structure")
mtext("Lag-response association")
legend("topright",c("True","Avg simul", "Indiv simul"), col=c(1,2,grey(0.9)), 
  lty=c(1,2,1), lwd=2, inset=0.02, y.intersp=1.3, cex=0.8)

################################################################################
# FIGURE S6

# RUN FIRST PART OF 11.outdeprisk.R, BEFORE THE PLOTTING
lines <- scan("11.outdeprisk.R",what="character",sep="\n")
linescoll <- paste(lines[seq(which(lines == "# PLOT"))], collapse='\n')
source(textConnection(linescoll))

par(las=1, mar=c(5,4,4,1), mgp=c(2.5,1,0))
layout(t(1:2))

plot(date, date, ylim=c(0.5,10+0.5), yaxt="n", ylab="", xlab="Follow-up",
  frame.plot=F)
axis(1, at=range(date)+c(-500,500), labels=F, lwd.ticks=0)
axis(2, at=rev(seq(10)), labels=paste("Sub",seq(10)), cex.axis=0.8, lwd=0, las=1)
title("Scenarios 11: Outcome-dependent risk")
mtext("Change of risk status dependent on outcome")
for(i in seq(10)) {
  sub <- subset(data, id==rev(seq(10))[i])
  lines(date, rep(i, length(date)), lty=2)
  lines(date[!sub$change], rep(i, sum(!sub$change)), lty=1)
  for(t in seq(length(sub$date))[sub$y>0]) {
    points(rep(sub$date[t], sub$y[t]), i+(seq(sub$y[t])-mean(seq(sub$y[t])))/4,
      pch=21, cex=1.5, bg=2)
  }
}

# RUN FIRST PART OF 12.outdepfu.R, BEFORE THE PLOTTING
lines <- scan("12.outdepfu.R",what="character",sep="\n")
linescoll <- paste(lines[seq(which(lines == "# PLOT"))], collapse='\n')
source(textConnection(linescoll))

plot(date, date, ylim=c(0.5,10+0.5), yaxt="n", ylab="", xlab="Follow-up",
  frame.plot=F)
axis(1, at=range(date)+c(-500,500), labels=F, lwd.ticks=0)
axis(2, at=rev(seq(10)), labels=paste("Sub",seq(10)), cex.axis=0.8, lwd=0, las=1)
title("Scenario 12: Outcome-dependent follow-up")
mtext("Censoring process linked with outcome")
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
