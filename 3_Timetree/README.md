# Divergence time estimation

1. Selecte the most clock-like genes using Sortadate

The input include rooted gene trees and one rooted species tree
```
python ~/programs/SortaDate/src/get_var_length.py --flend .rooted.tre --outf root2tip_var --outg FC230_c98,Cercis_canadensis_ge,FC232_c98 geneTr/ --loc ~/programs/phyx/bin/
python ~/programs/SortaDate/src/get_bp_genetrees.py geneTr/ legume.sortadate.sp.tre --flend .rooted.tre --outf bp
python ~/programs/SortaDate/src/combine_results.py ./root2tip_var ./bp --outf comb

python ~/programs/SortaDate/src/get_good_genes.py ./comb.filtered --max 40 --order 3,1,2 --outf clock.gene.list
```

2. Infer a BEAST time tree

This relatively large dataset of 30 genes (large for BEAST to converge) and 12 calibration points **requires** a starting tree, otherwise it struggles to find a proper starting points. To do this, we inferred a treePL tree compliant with all fossil calibration points `round2.spname.rooted.treePL.tre`. Then add this as the starting tree following BEAST [instrcutions](https://www.beast2.org/fix-starting-tree/). BEAST tree search was conducted with fixed topology was fixed, linked tree length among partitions, GTR substitution model, an uncorrelated relaxed molecular clock model, and a birth-death tree prior. The BEAST control file `15k_legume_fixedtopo.xml` contains all alignment, prior, and fossil constraint info.

4. Use the BEAST time tree as secondary calibrations for treePL tree

With a well-curated matK dataset, we infer an ML phylogeny with partially fixed relationships indicated in `round2.topo_constraint.tre`: relationships **between** the major Papilionoid clades were constrained including the ADA, Swartzieae, Cladrastis, Vataireoid, Exostyleae, Dalbergioid s.l., Genistoid s.l., Andira, Baphieae, and NPAAA. Relationships **within** each clade and **between all subfamilies and families** remained **as-is**. 
···
iqtree -s G1088.legumetimetr.fas -st DNA -m GTR -nt AUTO -pre legume2 -te round2.topo_constraint.tre
···

