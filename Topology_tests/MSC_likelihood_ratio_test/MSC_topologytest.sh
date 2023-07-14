#based on the optimum astral tree, create alternative topologies by alternating the position of focal clade
#then optimize branch length in astral tree
java -jar astral.5.7.8.jar -q concTr1_1.GenistoidsDalbergioidsmono.tre -i G1456.trees -o concTr1_1.astral.tre


#infer likelihood in MPE-EST. This is because ASTRAL does not output likelihood scores
