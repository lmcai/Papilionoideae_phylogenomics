from Bio import SeqIO

x=open('G672_len900_sp227.list').readlines()
a=[]
for i in x:
	y=open('../../4_cleaned_data_for_phylo/aln/'+i.strip()).readlines()
	y=[j[1:].strip() for j in y if j.startswith('>')]
	a=a+y


a=list(set(a))
seq={}
for i in a:
	seq[i]=''

b=1
out1=open('sp287.G672.fas','a')
out2=open('G672.partition','a')

for i in x:
	recs=SeqIO.index('../../4_cleaned_data_for_phylo/aln/'+i.strip(),'fasta')
	for rec in recs:
		seq_len=len(recs[rec].seq)
		break
	for j in a:
		try:seq[j]=seq[j]+str(recs[j].seq)
		except KeyError:seq[j]=seq[j]+'-'*seq_len
	out2.write(i.strip() + '='+str(b)+'-'+str(b+seq_len-1)+';\n')
	b=b+seq_len

for i in a:
	out1.write('>'+i+'\n'+seq[i]+'\n')