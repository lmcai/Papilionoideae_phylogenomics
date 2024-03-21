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

3. Use the BEAST time tree as secondary calibrations for treePL tree