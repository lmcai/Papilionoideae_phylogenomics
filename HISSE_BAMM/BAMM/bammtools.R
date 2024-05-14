install.packages("BAMMtools")
library(BAMMtools)

#assess convergence
tree <- read.tree("legume_matK.BEAST_treePL.tre")
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation,type='l')

#determine bur-in
burnstart <- floor(0.4 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
#218.4
effectiveSize(postburn$logLik)
#408.0
#make sure both values >200


#analyzing rate shifts

post_probs <- table(postburn$N_shifts) / nrow(postburn)
edata <- getEventData(tree, eventdata = "event_data.txt",burnin=0.4)
shift_probs <- summary(edata)

#bayes factor of shifts
bf <- computeBayesFactors(postburn, expectedNumberOfShifts=1)

#rate-through-time
plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="netdiv")

#phylorate plot
pdf('phylorate.pdf',width=10,height=10)
plot.bammdata(edata, lwd=2,  tau=0.001,breaksmethod='jenks', method='polar')
addBAMMshifts(edata, par.reset=FALSE, cex=1)
addBAMMlegend(q, location=c(0, 1, 140, 220))
dev.off()

#then mask high values
#plot.bammdata(edata, lwd=2,legend=T,color.interval = c(0.067,0.3))
#x=getTipRates(edata,returnNetDiv = TRUE,statistic = 'median') 
#xx=as.data.frame(x$netdiv.avg)