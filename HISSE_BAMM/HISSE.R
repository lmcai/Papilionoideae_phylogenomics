library(hisse)
library(phytools)
library(diversitree)
setwd('Dropbox/legume_phylogenomics/Diversification/')


tree<-read.tree("round2.spname.rooted.treePL.tre")
#The tree is not ultrametric due to rounding issues, force it to be an ultrametric
phy=force.ultrametric(tree, method="nnls")

trait.data<-read.csv("BISSE_trait.csv")
trait.dat=trait.data[,c("Species","Architecture_CODE")]

#check if name matches
name.check(phy, trait.dat)

################################################
##BISSE in hisse

#sampling proportion for the two states
f <- c(0.0251,0.0183)

trans.rates.bisse <-  TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)

#BISSE NULL
turnover <- c(1,1)
extinction.fraction <- c(1,1)

dull.null <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=FALSE, 
                     trans.rate=trans.rates.bisse)

dull.null

Fit 
            lnL             AIC            AICc 
      -1479.709        2967.418        2967.542 
         n.taxa n.hidden.states 
        326.000           1.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A 
1.030659083 1.030659083 0.871233151 0.871233151 
      q0A1A       q1A0A 
0.002002025 0.001315571 

#BISSE
turnover <- c(1,2)
extinction.fraction <- c(1,1)
BiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=FALSE, 
                     trans.rate=trans.rates.bisse)

BiSSE

Fit 
            lnL             AIC            AICc 
      -1431.114        2872.229        2872.416 
         n.taxa n.hidden.states 
        326.000           1.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A 
0.413406650 0.857617567 0.783589437 0.783589437 
      q0A1A       q1A0A 
0.001313956 0.001894974 


################################################

#HISSE
#decouple turover rate for 0A, 1A, 0B, 1B
turnover <- c(1,2,3,4)
extinction.fraction <- rep(1, 4) 
trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
print(trans.rate.hisse)

HiSSE <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate.hisse)

Fit 
            lnL             AIC            AICc 
      -1418.002        2856.005        2856.703 
         n.taxa n.hidden.states 
        326.000           2.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A 
1.184479e-01 3.152980e-01 7.502602e-01 7.502602e-01 
       q0A1A        q1A0A        q0A0B        q1A1B 
2.064378e-09 2.061154e-09 3.759845e-02 3.759845e-02 
  turnover0B   turnover1B        eps0B        eps1B 
4.829065e-01 9.047118e-01 7.502602e-01 7.502602e-01 
       q0B1B        q1B0B        q0B0A        q1B1A 
1.769959e-03 2.484168e-03 3.759845e-02 3.759845e-02 



#plotting
hisse_figure <- MarginReconHiSSE(phy=phy, data=trait.dat, f=f,par=HiSSE$solution, hidden.states=1,AIC=HiSSE$AIC)
plot_hisse <- plot.hisse.states(hisse_figure, rate.param = "net.div", show.tip.label = FALSE)


#CID-2 (null HISSE)
turnover <- c(1, 1, 2, 2)
extinction.fraction <- rep(1, 4) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
trans.rates = ParDrop(trans.rates, c(3,5,8,10))

CID2 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)

CID2

Fit 
            lnL             AIC            AICc 
      -1433.807        2879.614        2879.878 
         n.taxa n.hidden.states 
        326.000           2.000 

Model parameters: 

 turnover0A  turnover1A       eps0A       eps1A 
0.340838325 0.340838325 0.657132315 0.657132315 
      q0A1A       q1A0A       q0A0B       q1A1B 
0.002017159 0.001276213 0.001272090 0.001272090 
 turnover0B  turnover1B       eps0B       eps1B 
0.717930550 0.717930550 0.657132315 0.657132315 
      q0B1B       q1B0B       q0B0A       q1B1A 
0.002017159 0.001276213 0.001272090 0.001272090 


#CID-4
turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
extinction.fraction <- rep(1, 8) 
trans.rate <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
CID4 <- hisse(phy=phy, data=trait.dat, f=f, turnover=turnover, 
                     eps=extinction.fraction, hidden.states=TRUE, 
                     trans.rate=trans.rate)

Fit 
            lnL             AIC            AICc 
      -1410.614        2837.228        2837.682 
         n.taxa n.hidden.states 
        326.000           4.000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A 
0.5092584037 0.5092584037 0.4473014127 0.4473014127 
       q0A1A        q1A0A        q0A0B        q0A0C 
0.0019845940 0.0012557483 0.0007747367 0.0007747367 
       q0A0D        q1A1B        q1A1C        q1A1D 
0.0007747367 0.0007747367 0.0007747367 0.0007747367 
  turnover0B   turnover1B        eps0B        eps1B 
0.0820493311 0.0820493311 0.4473014127 0.4473014127 
       q0B1B        q1B0B        q0B0A        q0B0C 
0.0019845940 0.0012557483 0.0007747367 0.0007747367 
       q0B0D        q1B1A        q1B1C        q1B1D 
0.0007747367 0.0007747367 0.0007747367 0.0007747367 
  turnover0C   turnover1C        eps0C        eps1C 
0.2761599701 0.2761599701 0.4473014127 0.4473014127 
       q0C1C        q1C0C        q0C0A        q0C0B 
0.0019845940 0.0012557483 0.0007747367 0.0007747367 
       q0C0D        q1C1A        q1C1B        q1C1D 
0.0007747367 0.0007747367 0.0007747367 0.0007747367 
  turnover0D   turnover1D        eps0D        eps1D 
0.1881450334 0.1881450334 0.4473014127 0.4473014127 
       q0D1D        q1D0D        q0D0A        q0D0B 
0.0019845940 0.0012557483 0.0007747367 0.0007747367 
       q0D0C        q1D1A        q1D1B        q1D1C 
0.0007747367 0.0007747367 0.0007747367 0.0007747367 
