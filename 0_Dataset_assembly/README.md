# Dataset assembly

1. Transcriptome assembly

RNA extraction and transcriptome assembly pipeline followed Zhang et al. (2013) "Comparative analyses of two Geraniaceae transcriptomes using next-generation sequencing".


2. Extracting nuclear loci from genome skimming data using PhyloHerb

a. A reference fasta file `ref.fas` was prepared to include one reference sequence per locus from Zhao et al. (2019).

b. Execute the following PhyloHerb command to extract these loci and remove species with larger than 70% missings data.

```
python phyloherb.py -m assemb -r1 R1.fq -r2 R2.fq -ref ref.fas -prefix sp_A -n 8
python phyloherb.py -m ortho -i assemblies/ -o seq_extract_output -ref ref.fas -nuc -missing 0.7

```

3. Build a quick and dirty ML tree to identify paralogs with MAFFT and IQTREE

```
mafft input.fas >output.aln.fas
iqtree2 -s output.aln.fas
```

4. Prune exceptionally long branches using `prune_long_branches.py` and applied the Yang and Smith (2016) pipeline to remove paralogs or loci with less than 199 (70%) species.

```
python prune_paralogs_RT.py input_gene_trees/ *.treefile output_dir 199 ingroup_outgroup_def.txt
```

5. Realign with MAFFT-ensi, trim with trimal, then remove gaps unique to Zhao et al's data, finally infer a gene tree.

```
mafft --genafpair --maxiterate 1000 input.fas >output.aln.fas
trimal -gt 0.7 output.aln.fas >output.aln.fas
```

Remove gappy sites unique to Zhao's data using `SeqTrimmingByConsistency.R`

Then infer gene tree with IQ-TREE
```
iqtree2 -s output.aln.fas -B 1000 -bnni
```
