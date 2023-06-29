#gene concordance factor
iqtree2 -t raxmlng4gCF.tre --gcf G1456.trees --prefix G1456onRAxML.gCF
#site concordance factor
iqtree2 -te raxmlng.tre -s sp287.G672.fas --scfl 200 --prefix G672onRaxML.sCF


#AU, KH, Sh tests
iqtree -s sp287.G672.fas -m GTR+F+R4 -g raxmlng.tre --prefix conTr1_1.DalbergioidGenistoid -nt 3
iqtree -s sp287.G672.fas -m GTR+F+R4 -g conTr1_2.DalbergioidNPAAA.tre --prefix conTr1_2.DalbergioidNPAAA -nt 3
iqtree -s sp287.G672.fas -m GTR+F+R4 -z conTr1.trees -n 0 -zb 1000 -au -zw

#gene-wise likelihood comparisons fo two topologies.
iqtree2 -s sp287.G672.fas -p G672.partition -z conTr1.trees -n 0 -wpl --prefix conTr1.genewisell.wpl -m GTR+F+R4