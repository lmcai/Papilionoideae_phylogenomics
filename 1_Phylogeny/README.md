# Phylogeny

1. To infer a species tree under the coalescent model, we used the full dataset of 1456 gene in ASTRAL

Astral with ‘informative’ branches (>70 UFBP); low-support branches are collapsed using newick-utils

```
nw_ed  G2299.gene.trees 'i & b<=70' o > G2299.gene.70shalrt.trees
```
Then use ASTRAL-III to infer a species tree under the coalescent model
```
java -jar astral.5.7.8.jar -i G2299.gene.70shalrt.trees -o astral.ufbp70.rooted.tre
```

2. To infer a species tree using the concatenated DNA sequences, we first 
