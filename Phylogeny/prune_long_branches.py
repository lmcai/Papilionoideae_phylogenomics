from ete3 import Tree
import os
from scipy import stats

treefiles=os.listdir('.')
treefiles=[i for i in treefiles if i.endswith('.treefile')]

jdata=open('jansen_lab_data.txt').readlines()
jdata=[i.strip() for i in jdata]

duptaxa=open('dupsp2remove.txt').readlines()
duptaxa=[i.strip() for i in duptaxa]

def prune_long_br(tr,length_limit):
	#calculate trimmed mean branch length
	br_len=[node.dist for node in tr.traverse("postorder")]
	mean_br = stats.trim_mean(br_len, 0.05)
	#first round set to >15 times, second time set to >10 times
	tips2remove=[]
	for node in tr.traverse("postorder"):
		if node.dist>length_limit*mean_br:tips2remove=tips2remove+[j.name for j in node]
	tips2remove=list(set(tips2remove))
	#do not remove anything in our data yet
	#tips2remove=[j for j in tips2remove if not j in jdata]
	#remove duplicate taxa
	tips2remove=tips2remove+[node.name for node in tr if node.name in duptaxa]
	if len(tips2remove)>0:
		taxa=[node.name for node in tr if not node.name in tips2remove]
		tr.prune(taxa,preserve_branch_length =True)
	return(tr)

paralog_genes={}
for i in treefiles:
	tree=Tree(i)
	num_taxa_before=len(tree)
	#two rounds of long branch pruning; first remove branches longer than 15 times then average; second round 10 times
	tree = prune_long_br(tree,15)
	#tree = prune_long_br(tree,10)
	num_taxa_after=len(tree)
	paralog_genes[i]=num_taxa_before-num_taxa_after
	#tree.write(outfile=i.split('.')[0]+'.nolongbr.tre')
	tree.write(outfile=i.split('.')[0]+'.nolongbrnodup.tre')


###################
#get summary statistics
output=open('num_paralog_per_gene.tsv','a')
for key in paralog_genes.keys():
	output.write(key+'\t'+str(paralog_genes[key])+'\n')

output.close()

treefiles=os.listdir('.')
treefiles=[i for i in treefiles if i.endswith('.nolongbrnodup.tre')]
a={}
b={}
for i in treefiles:
    t=Tree(i)
    a[i]=len(t)
    for node in t:
            try:
                    b[node.name]=b[node.name]+1
            except KeyError:
                    b[node.name]=1

out1=open('num_sp_per_gene_after_pruning.tsv','a')
out2=open('number_gene_per_sp.tsv','a')
for key in a.keys():
    d=out1.write(key+'\t'+str(a[key])+'\n')

out1.close()

for key in b.keys():
    d=out2.write(key+'\t'+str(b[key])+'\n')

out2.close()

#################
#extract seq according to 
from Bio import SeqIO

treefiles=os.listdir('.')
#treefiles=[i for i in treefiles if i.endswith('.nolongbr.tre')]
treefiles=[i for i in treefiles if i.endswith('.nolongbrnodup.tre')]

for i in treefiles:
	tr=Tree(i)
	taxa=[j.name for j in tr]
	#recs=SeqIO.parse(i.split('.')[0]+'.aln.fas','fasta')
	recs=SeqIO.parse(i.split('.')[0]+'.nolongbr.Rtrim.aln.fas','fasta')
	#out=open(i.split('.')[0]+'.nolongbr.aln.fas','a')
	out=open(i.split('.')[0]+'.nolongbr.Rtrim.nodup.fas','a')
	for rec in recs:
		if rec.id in taxa:
			d=SeqIO.write(rec,out,'fasta')
	out.close()

