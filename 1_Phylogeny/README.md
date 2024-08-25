# Phylogeny

## Coalescent analysis

1. To infer a species tree under the coalescent model, we used the full dataset of 1456 gene in ASTRAL

Astral with ‘informative’ branches (>70 UFBP); low-support branches are collapsed using newick-utils

```
nw_ed  G2299.gene.trees 'i & b<=70' o > G2299.gene.70shalrt.trees
```
2. Then use ASTRAL-III to infer a species tree under the coalescent model
```
java -jar astral.5.7.8.jar -i G2299.gene.70shalrt.trees -o astral.ufbp70.rooted.tre
```

## Concatenation analysis

1. To infer a species tree using the concatenated DNA sequences, we first filter the genes based on length and missing data. We only used genes with >227 species (80%) and longer than 900 bp. This result in 672 loci listed in `G672_len900_sp227.list`.

Then we generate a gene concatenation file with `PhyloHerb`
```
python phyloherb.py -m conc -i G672 -o G672.conc -suffix .na.aln.fas
```

2. Run partitionfinder to get the optimal partition scheme with the configureation file `partition_finder.cfg`
```
python PartitionFinder.py -f -p 8 -r --rcluster-percent=10 --rcluster-max=1000 G1553_pf
```
3. With the optimal partition scheme is `G672.partitionfinder.scheme`, run iqtree for species tree inference with 100 bootstrap
```
iqtree2 -s sp287.G672.phy -p G672.partitionfinder.scheme -n AUTO -b 100
```
