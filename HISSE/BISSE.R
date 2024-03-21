#First navigate to the folder where you put the tree file 'BisseRtree.nex' and the data file 'BisseRdata.txt'
setwd('Dropbox/legume_phylogenomics/Diversification/')
#load the 'diversitree' package into R's working memory

library(diversitree)
library(phytools)

#Load the tree file
tree<-read.tree("round2.spname.rooted.treePL.tre")

#The tree is not ultrametric due to rounding issues, force it to be an ultrametric
new.tree=force.ultrametric(tree, method="nnls")

#Load the data file
data<-read.csv("BISSE_trait.csv",row.names=1)

#Convert the data into a vector that is a component of the tree file, so that diversitree can work with it.

data.v<-data[,2]
names(data.v)<-row.names(data)
new.tree$tip.state<-data.v


#Define the bisse equation we will analyze initially with maximum likelihood and then with a Bayesian MCMC.
#The sampling.f tells the equation the proportion of species with each state (0 and 1) that we have sampled from the total number of extant species.  I have roughly estimated those in this example.

equation<-make.bisse(tree=new.tree,states=new.tree$tip.state,sampling.f=c(0.2,0.7))



##################################################
#Fitting the equation above via maximum likelihood requires a heuristic starting point for parameter values.  
#This function creates these heuristic starting values, here defined as object 'p'.
p<-starting.point.bisse(new.tree)


#Now we estimate the maximum likelihood solution with all 6 parameters unconstrained
mle.fit<-find.mle(equation,p,method="subplex")


#Take a look at the maximum likelihood parameter estimates
mle.fit$par
lambda0      lambda1          mu0          mu1 
4.113151e-02 4.864682e-02 5.589475e-06 1.946241e-06 
         q01          q10 
5.560346e-04 4.129330e-03 

mle.fit$lnLik
[1] -1471.226

#Model test: does two states have different rates?
#constrain equal rates
equation.l <- constrain(equation, lambda1 ~ lambda0)
equation.l <- constrain(equation.l, mu1 ~ mu0)

#then start the ML search again:
#Note that the statement “p[argnames(lik.l)]” drops the λ1 element from the starting parameter vector). This fit has quite different parameters to the full model (compare μ0)
mle.fit.l <- find.mle(equation.l, p[argnames(equation.l)])
mle.fit.l$lnLik
[1] -1472.635

#compare
round(rbind(full=coef(mle.fit), equal.l=coef(mle.fit.l, TRUE)), 3)

      lambda0 lambda1   mu0 mu1   q01   q10
full      0.041   0.049 0.000   0 0.001 0.004
equal.l   0.046   0.046 0.003   0 0.001 0.003

#ANOVA test: the TRUE argument forces coef to return values for constrained parameters). However, the difference in fits is not statistically supported.
anova(mle.fit, equal.l=mle.fit.l)
        Df   lnLik    AIC ChiSq Pr(>|Chi|)
full     6 -1471.2 2954.4                 
equal.l  4 -1472.6 2953.3 2.819     0.2443


##################################################
#To run the Bayesian MCMC, we also need a distribution for the prior probabilities of parameter values.  Here we make the prior exponentially distributed

prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))


#Now we can run our preliminary MCMC with a 'tuning parameter' (w) set somewhat arbitrarily to 0.1.  The tuning parameter defines how much the MCMC process varies the parameter values in each step.  We will run a second and final MCMC after this one with the tuning parameter based on the preliminary posterior parameter distributions.

mcmc.bisse<-mcmc(equation,mle.fit$par,nsteps=1000,prior=prior,w=0.1)


#Take a look at the output

mcmc.bisse


#The tuning parameter for the final MCMC should be equal to approximately the width of the middle 90% of the posterior samples for each parameter value from the preliminary MCMC.  We can calculate such a tuning parameter this way:

w=diff(sapply(mcmc.bisse[2:7],quantile,c(0.05,0.95)))


#Take a look at the new tuning parameters.  Now there are separate ones for each of the 6 model parameters to reflect the differences among them in their distributions

w


#Run the final MCMC

mcmc.bisse2<-mcmc(equation,mle.fit$par,nsteps=10000,w=w,prior=prior)

#Now we can examine the 95% credible intervals of the posterior samples for each parameter.  If the intervals do not overlap, then we have posterior Bayesian support for a difference in rates
sapply(mcmc.bisse2[,2:7],quantile,c(0.025,0.975))

#plot
col <- c("blue", "red")
profiles.plot(mcmc.bisse2[, c("lambda0", "lambda1")], col.line = col, las = 1, legend = "topright")
profiles.plot(mcmc.bisse2[, c("mu0", "mu1")], col.line = col, las = 1, legend = "topright")

