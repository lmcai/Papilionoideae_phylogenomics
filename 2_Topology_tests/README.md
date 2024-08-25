# Testing alternative placements of clades

## Likelihood-based tests: approximately unbiased (AU), Kishino-Hasegawa (KH), and Shimodaira-Hasegawa (SH) tests

1. Prepare multiple species trees (two competing hypotheses in the example below) by alternating the placements of focal group, and then infer branch length in iqtree, concatenated into one file. Commands:
```
iqtree2 -s sp287.G672.fas -m GTR+F+R4 -g global_best.tre --prefix conTr1_1.DalbergioidGenistoid
iqtree2 -s sp287.G672.fas -m GTR+F+R4 -g conTr1_2.DalbergioidNPAAA.tre --prefix conTr1_2.DalbergioidNPAAA
cat conTr1_1.DalbergioidGenistoid.treefile conTr1_2.DalbergioidNPAAA.treefile >conTr1.trees
```
2. Implement AU, KH, SH tests in IQ-TREE with the following commands:
```
#AU, KH, Sh tests
iqtree2 -s sp287.G672.fas -m GTR+F+R4 -z conTr1.trees -n 0 -zb 1000 -au -zw
```

## Gene and site concordance factor
1. Commands to get gCF and sCF in iqtree:
```
#gene concordance factor
iqtree2 -t iqtree4gCF.tre --gcf G1456.trees --prefix G1456onRAxML.gCF
#site concordance factor
iqtree2 -te iqtree.tre -s sp287.G672.fas --scfl 200 --prefix G672onRaxML.sCF
```

Result explained: The gDF quantifies support for the two nearest-neighbor interchange bipartitions (gDF1 and gDF2) and for all other possible
topologies (gDFP, as these are paraphyletic relative to the species
tree bipartition).
The sDF metrics quantify support among sites for the two possible alternative quartets (sDF1 and sDF2). Low gDF1 and gDF2 values or high gDFP values suggest the gene trees or alignments lack a clear signal, as do sDF values close to 33% (Minh et al., 2020a).

## Gene-wise signals.
Now we want to investigate the cause for topological difference between various phylogenomic analysis. One way is to identify genes contributing most phylogenetic signal towards one tree but not the other as described in Shen et al. (2017, Nat Ecol Evol). How can one do this? We can look at the gene-wise log-likelihood (logL) differences between the two given trees T1 and T2. Those genes having the largest logL(T1)-logL(T2) will be in favor of T1. Whereas genes showing the largest logL(T2)-logL(T1) are favoring T2.

To compute gene-wise log-likelihoods for the two trees, you can use the -wpl option (for writing partition log-likelihoods):
```
iqtree2 -s sp287.G672.fas -p G672.partition -z conTr1.trees -n 0 -wpl --prefix conTr.genewisell.wpl -m GTR+F+R5
```
will write a file `conTr.genewisell.wpl.partlh`, that contains log-likelihoods for all partitions in the original partition file.

Import `conTr.genewisell.wpl.partlh` into MS Excel, R, or any other spreadsheet software to examine the contribution of each locus to the final likelihoods. To visualize the results, use the R script `genewiseLL_plot.R` 

