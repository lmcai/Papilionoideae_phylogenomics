1. Manually prepare alternative subsampled species trees with alternative topologies

2. Optimize their branch length and infer the likelihood in MPEST

3. Simulate 100 bootstrap sets of 1456 gene trees under the optimum species tree using the `sim.coal.mpest` in the R package `Phybase`

```
library(phybase)
SpTr_text=readLines('G1456.ufbp70.mpestBrLen.tre')
geneTr_sim_text=sim.coal.mpest(SpTr_text,145600)
write(paste(geneTr_sim_text,collapse='\n'),'')
```

4. Calculate the likelihood of the two alternative topologies under the 100 bootstrap gene tree sets. Examine if the difference of the likelihoods are significant.