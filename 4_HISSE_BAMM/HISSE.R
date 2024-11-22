library(hisse)
library(phytools)
library(diversitree)
setwd('Dropbox/legume_phylogenomics/Diversification/matK')


tree<-read.tree("legume_matK.BEAST_treePL.tre")
#The tree is not ultrametric due to rounding issues, force it to be an ultrametric
phy=force.ultrametric(tree, method="nnls")

trait.data<-read.csv("HISSE.csv")
trait.dat=trait.data[,c("Tip_name","Floral_architecture_code")]


#check if name matches
name.check(phy, trait.dat)

################################################
##1. BISSE in hisse

#sampling proportion for the two states
f <- c(0.1768,0.123)

trans.rates.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)

#Did not run BISSE null model
#BISSE NULL
#turnover <- c(1,1)
#extinction.fraction <- c(1,1)

#dull.null <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
#                     eps=extinction.fraction, hidden.states=FALSE, 
#                     trans.rate=trans.rates.bisse)

#dull.null

#Fit 
#            lnL             AIC            AICc 
#      -1479.709        2967.418        2967.542 
#         n.taxa n.hidden.states 
#        326.000           1.000 

#Model parameters: 

# turnover0A  turnover1A       eps0A       eps1A 
#1.030659083 1.030659083 0.871233151 0.871233151 
#      q0A1A       q1A0A 
#0.002002025 0.001315571 

#BISSE
turnover <- c(1,2)
extinction.fraction <- c(1,1)
BiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=FALSE, 
                     trans.rate=trans.rates.bisse)

Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -9859.611       19729.222       19729.240        3336.000           1.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A       q0A1A       q1A0A 
2.875057362 3.259826568 0.949865814 0.949865814 0.002498468 0.001177123 



################################################

#2. HISSE
#decouple turover rate for 0A, 1A, 0B, 1B
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)

HiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)


Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -9400.047       18820.095       18820.161        3336.000           2.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A        q0A0B        q1A1B 
0.0291997872 0.0932612343 0.9629306771 0.9629306771 0.0011344477 0.0013525362 0.2372463830 0.2372463830 
  turnover0B   turnover1B        eps0B        eps1B        q0B1B        q1B0B        q0B0A        q1B1A 
6.6652247045 7.1864873173 0.9629306771 0.9629306771 0.0031944077 0.0009583496 0.2372463830 0.2372463830 




#plotting

#Ancestral State Estimation based on Marginal Reconstruction for the HiSSE model.
hisse_figure <- MarginReconHiSSE(phy=phy, data=trait.dat, f=f,par=HiSSE$solution, hidden.states=1,AIC=HiSSE$AIC)
plot_hisse <- plot.hisse.states(hisse_figure, rate.param = "net.div", show.tip.label = FALSE)

#Ancestral State Estimation for BiSSE
bisse_figure <- MarginReconHiSSE(phy=phy, data=trait.dat, f=f,par=BiSSE$solution, hidden.states=0,AIC=BiSSE$AIC)
plot_hisse <- plot.hisse.states(bisse_figure, rate.param = "net.div", show.tip.label = FALSE)

############################
#3. CID-2 (null HISSE)
turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate = ParDrop(trans.rate, c(3,5,8,10))

CID2 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)

CID2

Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -9863.327       19736.654       19736.672        3336.000           2.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A       q0A1A       q1A0A 
0.006208926 0.006208926 0.949894298 0.949894298 0.002406056 0.001231213 
 turnover0B  turnover1B       eps0B       eps1B       q0B1B       q1B0B 
3.125897725 3.125897725 0.949894298 0.949894298 0.002406056 0.001231213 

############################
#4. CID-4
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
CID4 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)

Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -9317.661       18651.323       18651.366        3336.000           4.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A       q0A1A       q1A0A       q0A0B       q0A0C 
8.058328270 8.058328270 0.945715423 0.945715423 0.002163516 0.001171302 0.106109050 0.106109050 
      q0A0D       q1A1B       q1A1C       q1A1D  turnover0B  turnover1B       eps0B       eps1B 
0.106109050 0.106109050 0.106109050 0.106109050 0.037168531 0.037168531 0.945715423 0.945715423 
      q0B1B       q1B0B       q0B0A       q0B0C       q0B0D       q1B1A       q1B1C       q1B1D 
0.002163516 0.001171302 0.106109050 0.106109050 0.106109050 0.106109050 0.106109050 0.106109050 
 turnover0C  turnover1C       eps0C       eps1C       q0C1C       q1C0C       q0C0A       q0C0B 
1.043019320 1.043019320 0.945715423 0.945715423 0.002163516 0.001171302 0.106109050 0.106109050 
      q0C0D       q1C1A       q1C1B       q1C1D  turnover0D  turnover1D       eps0D       eps1D 
0.106109050 0.106109050 0.106109050 0.106109050 1.011376297 1.011376297 0.945715423 0.945715423 
      q0D1D       q1D0D       q0D0A       q0D0B       q0D0C       q1D1A       q1D1B       q1D1C 
0.002163516 0.001171302 0.106109050 0.106109050 0.106109050 0.106109050 0.106109050 0.106109050 


############################
#5. MISSE (trait-free, tip rates)
turnover <- c(1)
eps <- c(1)
one.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=19302

turnover <- c(1,2)
eps <- c(1,1)
two.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18400

#rate classes A:C
turnover <- c(1,2,3)
eps <- c(1,1,1)
three.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18279

#rate classes A:D
turnover <- c(1,2,3,4)
eps <- c(1,1,1,1)
four.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=

#rate classes A:E
turnover <- c(1,2,3,4,5)
eps <- c(1,1,1,1,1)
five.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=

################################
#6. Papilionoideae-specific analysis
tree<-read.tree("papilionoideae_matK.BEAST_treePL.tre")
f <- c(0.1768,0.123)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)

HiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)


turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
CID4 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)
