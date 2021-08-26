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
# RUN ALL THE SCRIPTS
################################################################################

# LIST OF SCENARIOS
scen <- c("01.basic","02.rare","03.expcont","04.outbin","05.outcont",
  "06.trendcomm", "07.trendsub","08.baseconf","09.timevarconf","10.complexlag",
  "11.outdeprisk","12.outdepfu","13.outdepexp","14.varbaserisk")

# CHANGE THE NUMBER OF SIMULATIONS IN THE SCRIPTS
# NB: TO REPRODUCE THE RESULTS EXACTLY, UNCOMMENT AND SET THE NUMBER TO 50000
# for(j in seq(scen)){
#   x <- readLines(paste(scen[j],"R",sep="."))
#   y <- gsub("nsim <- 10000#","nsim <- 50000#", x, fixed=T)
#   cat(y, file=paste(scen[j],"R",sep="."), sep="\n")
# }
# rm(j,x,y)

################################################################################
# RUN ALL THE SCENARIOS

for(j in seq(scen)) {
  cat(scen[j],"\n")
  source(paste(scen[j],"R",sep="."))
  obj <- c("nsim", "fwlag", "beta", "trueeff", "eff", "efflagmean", "efflagsample",
    "bias", "cov", "rmse")
  save(list=ls()[ls() %in% obj], file=paste0("results/",scen[j],".RData"))
}

################################################################################
# RESULTS

# PRODUCE THE TABLE
tab <- t(sapply(seq(scen), function(j) {
  load(paste0("results/",scen[j],".RData"))
  c(paste0(formatC(bias, format="f", digits=1),"%"),
    formatC(cov, format="f", digits=3),
    paste0(formatC(rmse, format="f", digits=1),"%")
  )
}))

# NAMES AND PRINT
dimnames(tab) <- list(scen, c("Bias", "Coverage", "RMSE"))
tab
