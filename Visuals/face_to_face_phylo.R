setwd('MyDir/')
library(ape)
#association=read.table('pt_ITS.filtered.tsv')
association=read.csv('astral_raxml.csv',sep=',',header = T)
tree1=read.tree('~/Downloads/raxmlng.rooted.tre')
tree2=read.tree('~/Downloads/astral.rooted.tre')
pdf(file = 'raxml_astral.pdf',width = 150,height = 80)
cophyloplot(tree1, tree2, assoc = association,length.line = 0, space = 120, gap = 0,use.edge.length =F)
dev.off()