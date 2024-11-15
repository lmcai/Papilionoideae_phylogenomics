# Diversification rate

## BAMM to detect shifts in diversification rates

See R script `bammtools.R` for details

## Trait-dependent diversification

BiSSE and HiSSE models as well as the corresponding trait-independent models (CD-2 and CD-4) can be used to evaluate whether the evolution of certain traits is correlated with shifts in diversification rate.

Tutorials on this topic can be found on the following websites: 

https://speciationextinction.info/articles/hisse-new-vignette.html

https://lukejharmon.github.io/ilhabela/2015/07/05/BiSSE-and-HiSSE/

https://roszenil.github.io/portfolio/SSBworkshop/ For model comparison.

To implement this test, use the script in `HiSSE.R`. The reference phylogeny is provided in `legume_matK.BEAST_treePL.tre` and the trait matrix is `HISSE.csv`. Notably, the following criteria can be used to evaluate hypotheses using AIC values.

* If ΔAIC < 2, then there is substantial support for the model with the lower AIC score.
* If 2 < ΔAIC < 6, then there is less support for the model with the higher AIC score.
* If ΔAIC > 6, then there is essentially no support for the model with the higher AIC score.

