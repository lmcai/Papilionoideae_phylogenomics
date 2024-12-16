# Diversification rate

## BAMM to detect shifts in diversification rates

See R script `bammtools.R` for details. The posterior distribution of rate shift reconstructions is stored in `bamm.event_data.txt`.

To generate Fig. 2, we used `bammtools_plot_script.R` to plot the net diversification rate estimated by BAMM, tip rates estimated by MiSSE, and the keel or non-keel floral trait.

## Trait-dependent diversification

BiSSE and HiSSE models as well as the corresponding trait-independent models (CD-2 and CD-4) can be used to evaluate whether the evolution of certain traits is correlated with shifts in diversification rate. We also implemented the trait-free MiSSE model as well.

Tutorials on this topic can be found on the following websites: 

https://speciationextinction.info/articles/hisse-new-vignette.html

https://lukejharmon.github.io/ilhabela/2015/07/05/BiSSE-and-HiSSE/

https://roszenil.github.io/portfolio/SSBworkshop/ For model comparison.

To implement this test, use the script in `macroevolution_hypo_test.R`. The reference phylogeny is provided in `legume_matK.BEAST_treePL.tre` and the trait matrix is `HISSE.csv`. Notably, the following criteria can be used to evaluate hypotheses using AIC values.

* If ΔAIC < 2, then there is substantial support for the model with the lower AIC score.
* If 2 < ΔAIC < 6, then there is less support for the model with the higher AIC score.
* If ΔAIC > 6, then there is essentially no support for the model with the higher AIC score.

