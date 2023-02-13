#!/bin/bash
#SBATCH -J bowtie           # Job name
#SBATCH -o bowtie.o%j       # Name of stdout output file
#SBATCH -e bowtie.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 24:00:00        # Run time (hh:mm:ss)

#module load intel/17.0.4
#module load bowtie/2.3.2

#module load samtools/1.10
#source activate getorg

bowtie2 --very-sensitive-local -p 32 -x G1559 -1 ../jansen_lab_genome_skimming/$1/$2 -2 ../jansen_lab_genome_skimming/$1/$3 -S $1\_G1559.sam
#get mapped reads
samtools view -b -F 4 $1\_G1559.sam | samtools bam2fq >$1\.G1559_mapped.fq

spades.py -t 32 --phred-offset 33 --s1 $1\.G1559_mapped.fq -k 21,55,85,115 -o $1\_spades
