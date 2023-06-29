#gene concordance factor
iqtree2 -t raxmlng4gCF.tre --gcf G1456.trees --prefix G1456onRAxML.gCF
#site concordance factor
iqtree2 -te raxmlng.tre -s sp287.G672.fas --scfl 200 --prefix G672onRaxML.sCF


#AU, KH, Sh tests
iqtree -s sp287.G672.fas -m GTR+F+R6 -g raxmlng.tre --prefix conTr1_1.DalbergioidGenistoid -nt 3
iqtree -s sp287.G672.fas -m GTR+F+R6 -g conTr1_2.DalbergioidNPAAA.tre --prefix conTr1_2.DalbergioidNPAAA -nt 3
iqtree -s sp287.G672.fas -m GTR+F+R6 -z conTr1.trees -n 0 -zb 1000 -au -zw

