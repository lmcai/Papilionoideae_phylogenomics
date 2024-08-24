1. Gene and site concordance factor
A. Gene


 iqtree2 -t raxmlng4gCF.tre --gcf G1456.trees --prefix G1456onRAxML.gCF

G1456onRAxML.gCF.cf.tree.nex and  G1456onRAxML.gCF.cf.tree

Result explained: The gDF quantifies support for the two nearest-neighbor interchange bipartitions (gDF1 and gDF2) and for all other possible
topologies (gDFP, as these are paraphyletic relative to the species
tree bipartition).
The sDF metrics quantify support among sites for the two possible alternative quartets (sDF1 and sDF2). Low gDF1 and gDF2 values or high gDFP values suggest the gene trees or alignments lack a clear signal, as do sDF values close to 33% (Minh et al., 2020a).

B. Site
iqtree2 -te raxmlng.tre -s sp287.G672.fas --scfl 200 --prefix G672onRaxML.sCF
Results are in 

3. Gene-wise signals.
Now we want to investigate the cause for such topological difference between trees inferred by single and partition model. One way is to identify genes contributing most phylogenetic signal towards one tree but not the other.
How can one do this? We can look at the gene-wise log-likelihood (logL) differences between the two given trees T1 and T2. Those genes having the largest logL(T1)-logL(T2) will be in favor of T1. Whereas genes showing the largest logL(T2)-logL(T1) are favoring T2.

To compute gene-wise log-likelihoods for the two trees, you can use the -wpl option (for writing partition log-likelihoods):
```
iqtree2 -s sp287.G672.fas -p G672.partition -z conTr1.trees -n 0 -wpl --prefix conTr.genewisell.wpl -m GTR+F+R5
```
will write a file turtle.wpl.partlh, that contains log-likelihoods for all partitions in the original partition file. We use -p turtle.nex.best_scheme.nex here (instead of -p turtle.nex) to avoid doing model selection again.

Import turtle.wpl.partlh into MS Excel, Libre Office Calc, or any other spreadsheet software. You will need to tell the software to treat spaces as delimiters, so that the values are imported into different columns for easy processing (e.g., doing log-likelikehood subtraction as pointed out above).
