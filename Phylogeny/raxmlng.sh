/work2/08400/lmcai/stampede2/apps/raxml-ng --msa sp287.G672.phy --tree RAxML.parsimony.tre --threads 32 --model G672.partitionfinder.scheme

##############
#examl
#parse alignments
#/work2/08400/lmcai/stampede2/apps/ExaML-3.0.22/bin/parse-examl -s sp287.G672.phy -n sp287.G672 -m DNA -q G672.partitionfinder.scheme
#generate starting tree
#raxmlHPC-SSE3 -y -m GTRGAMMA -s sp287.G672.phy -n BP -p 3256179
#run examl
#mpirun -np 32 /work2/08400/lmcai/stampede2/apps/ExaML-3.0.22/bin/examl-OMP-AVX -S -s sp287.G672.binary -m GAMMA -n sp287.G672 -t RAxML.parsimony.tre

#generate BP alignments
#raxmlHPC-SSE3 -# 100 -b 12345 -f j -m GTRGAMMA -s sp287.G672.phy -n BP
#run BP
#/work2/08400/lmcai/stampede2/apps/ExaML-3.0.22/bin/parse-examl -s sp287.G672.phy.BS0 -n BP -m DNA
#rm RAxML_info.BP
#mpirun -np 32 examl-OMP-AVX -S -s BP.binary -m GAMMA -n BP -t RAxML.parsimony.tre

