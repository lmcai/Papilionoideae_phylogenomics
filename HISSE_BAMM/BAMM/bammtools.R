install.packages("BAMMtools")
library(BAMMtools)

#assess convergence
tree <- read.tree("sp214.pt_ITS.treePL.tre")
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation,type='l',ylim=c(-3300,-3150))

#determine bur-in
burnstart <- floor(0.5 * nrow(mcmcout))
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)



#analyzing rate shifts

post_probs <- table(postburn$N_shifts) / nrow(postburn)
edata <- getEventData(tree, eventdata = "event_data.txt",burnin=0.5)
shift_probs <- summary(edata)

#bayes factor of shifts
bf <- computeBayesFactors(postburn, expectedNumberOfShifts=1)

#rate-through-time
plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="netdiv")

#phylorate plot
plot.bammdata(edata, lwd=2)
#then mask high values
plot.bammdata(edata, lwd=2,legend=T,color.interval = c(0.067,0.3))
x=getTipRates(edata,returnNetDiv = TRUE,statistic = 'median') 
xx=as.data.frame(x$netdiv.avg)