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
##1. Fabales BISSE

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

#2. Fabales HISSE
#2 hidden state (standard)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)

HiSSE_two <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
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

#3 hidden states
turnover <- c(1, 2, 3, 4, 5, 6)
extinction.fraction <- rep(1, 6) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=2)
HiSSE_three <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)

HiSSE_three

Fit
            lnL             AIC            AICc          n.taxa n.hidden.states
      -9310.949       18643.898       18643.977        3326.000           3.000

Model parameters:

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A        q0A0B        q0A0C        q1A1B
7.426899e-01 7.426899e-01 9.541735e-01 9.541735e-01 3.877899e-03 3.921268e-03 1.474672e-01 1.474672e-01 1.474672e-01
       q1A1C   turnover0B   turnover1B        eps0B        eps1B        q0B1B        q1B0B        q0B0A        q0B0C
1.474672e-01 2.433228e-01 2.433228e-01 9.541735e-01 9.541735e-01 3.590887e-09 2.061154e-09 1.474672e-01 1.474672e-01
       q1B1A        q1B1C   turnover0C   turnover1C        eps0C        eps1C        q0C1C        q1C0C        q0C0A
1.474672e-01 1.474672e-01 7.874236e+00 7.874236e+00 9.541735e-01 9.541735e-01 2.536519e-03 2.061154e-09 1.474672e-01
       q0C0B        q1C1A        q1C1B
1.474672e-01 1.474672e-01 1.474672e-01 

#4 hidden states
turnover <- c(1, 2, 3, 4, 5, 6, 7, 8)
extinction.fraction <- rep(1, 8)
trans.rate <- TransMatMakerHiSSE(hidden.traits=3)
HiSSE_four <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)

HiSSE_four

Fit
            lnL             AIC            AICc          n.taxa n.hidden.states
      -9281.947       18599.893       18600.100        3326.000           4.000

Model parameters:

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A        q0A0B        q0A0C        q0A0D
8.836121e-01 8.077196e+00 9.459972e-01 9.459972e-01 9.448870e-03 2.061154e-09 1.044834e-01 1.044834e-01 1.044834e-01
       q1A1B        q1A1C        q1A1D   turnover0B   turnover1B        eps0B        eps1B        q0B1B        q1B0B
1.044834e-01 1.044834e-01 1.044834e-01 7.030263e-01 1.063395e+00 9.459972e-01 9.459972e-01 7.989273e-04 2.937228e-03
       q0B0A        q0B0C        q0B0D        q1B1A        q1B1C        q1B1D   turnover0C   turnover1C        eps0C
1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01 6.641535e-09 2.101791e-01 9.459972e-01
       eps1C        q0C1C        q1C0C        q0C0A        q0C0B        q0C0D        q1C1A        q1C1B        q1C1D
9.459972e-01 2.070454e-09 2.070469e-09 1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01 1.044834e-01
  turnover0D   turnover1D        eps0D        eps1D        q0D1D        q1D0D        q0D0A        q0D0B        q0D0C
7.880544e+00 1.189505e+00 9.459972e-01 9.459972e-01 2.064055e-09 2.756798e-03 1.044834e-01 1.044834e-01 1.044834e-01
       q1D1A        q1D1B        q1D1C
1.044834e-01 1.044834e-01 1.044834e-01 


############################
#3. Fabales CID-2 (null HISSE)
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
#4. Fabales CID-4
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


############################################################
#5. Fabales MISSE (trait-free, tip rates)
################################

#Greedy test for the best number of hidden rates
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
#AICc=18228

#rate classes A:E
turnover <- c(1,2,3,4,5)
eps <- c(1,1,1,1,1)
five.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18201

turnover <- c(1,2,3,4,5,6,7,8,9)
eps <- c(1,1,1,1,1,1,1,1,1)
nine.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18156


turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
eps <- rep(1,19)
nineteen.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AIC= 18137

turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
eps <- rep(1,20)
twenty.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18134.2


turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
eps <- rep(1,21)
twentyone.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18132.78

turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
eps <- rep(1,22)
twentytwo.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc= 18138.92


turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
eps <- rep(1,23)
twentythree.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc= 18136

turnover <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
eps <- rep(1,26)
twentysix.rate <- MiSSE(phy, f=1, turnover=turnover, eps=eps)
#AICc=18137


#misse_recon <- MarginReconMiSSE(phy = model.set[[model_index]]$phy, f = 1, 
#                                    hidden.states = nturnover, 
#                                    pars = model.set[[model_index]]$solution, 
#                                    AIC = model.set[[model_index]]$AIC)

#reconstruct tip rate with number of hidden.states=21
twentyone.rate.recon <- MarginReconMiSSE(phy=phy, f=1,  hidden.states=21, pars=twentyone.rate$solution, n.cores=4, AIC=twentyone.rate$AIC)
plot.misse.states(twentyone.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",fsize=.25, legend="none")
tip.rates <- GetModelAveRates(twentyone.rate.recon, type = c("tips"))


################################################################
#6. Papilionoideae-specific analysis
tree<-read.tree("papilionoideae_matK.BEAST_treePL.tre")
f <- c(0.1768,0.123)
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)
#HiSSE
HiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)

Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -6163.053       12346.107       12346.209        2157.000           2.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A 
6.442302e-01 2.061915e-09 9.499371e-01 9.499371e-01 5.348714e-03 1.734043e-03 
       q0A0B        q1A1B   turnover0B   turnover1B        eps0B        eps1B 
2.694958e-01 2.694958e-01 4.130285e+00 6.656230e+00 9.499371e-01 9.499371e-01 
       q0B1B        q1B0B        q0B0A        q1B1A 
3.156031e-02 2.096788e-09 2.694958e-01 2.694958e-01 


#BiSSE
turnover <- c(1,2)
extinction.fraction <- c(1,1)
BiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=FALSE, 
                     trans.rate=trans.rates.bisse)
Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -6397.342       12804.685       12804.712        2157.000           1.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A 
2.0989113285 3.1860417579 0.9366492254 0.9366492254 0.0188668477 0.0007871378 


#CID2
turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate = ParDrop(trans.rate, c(3,5,8,10))

CID2 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)

Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -6409.027       12828.054       12828.082        2157.000           2.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A 
3.1205264937 3.1205264937 0.9401849433 0.9401849433 0.0175096333 0.0008643827 
  turnover0B   turnover1B        eps0B        eps1B        q0B1B        q1B0B 
0.3168891373 0.3168891373 0.9401849433 0.9401849433 0.0175096333 0.0008643827 

#CID4
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
CID4 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)


Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -6127.718       12271.437       12271.504        2157.000           4.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A 
7.9563656913 7.9563656913 0.9382626326 0.9382626326 0.0180036837 0.0008305166 
       q0A0B        q0A0C        q0A0D        q1A1B        q1A1C        q1A1D 
0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 
  turnover0B   turnover1B        eps0B        eps1B        q0B1B        q1B0B 
0.7478200518 0.7478200518 0.9382626326 0.9382626326 0.0180036837 0.0008305166 
       q0B0A        q0B0C        q0B0D        q1B1A        q1B1C        q1B1D 
0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 
  turnover0C   turnover1C        eps0C        eps1C        q0C1C        q1C0C 
0.7444243819 0.7444243819 0.9382626326 0.9382626326 0.0180036837 0.0008305166 
       q0C0A        q0C0B        q0C0D        q1C1A        q1C1B        q1C1D 
0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 
  turnover0D   turnover1D        eps0D        eps1D        q0D1D        q1D0D 
0.8309344692 0.8309344692 0.9382626326 0.9382626326 0.0180036837 0.0008305166 
       q0D0A        q0D0B        q0D0C        q1D1A        q1D1B        q1D1C 
0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 0.1259020905 


