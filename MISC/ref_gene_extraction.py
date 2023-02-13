from Bio import SeqIO
recs=SeqIO.index('Lj2.5_cds.fas','fasta')
x=open('Gmax_Lotus_floral_genes.txt').readlines()
out=open('legume_floral_dev.ref.fas','a')
for l in x[:16]:
	y=l.split()
	for i in y[1:]:
			try:
				d=out.write('>'+y[0]+'_'+i+'\n')
				d=out.write(str(recs[i].seq)+'\n')
			except KeyError:
				print(i)

recs=SeqIO.index('Gmax_189_cds_primaryTranscriptOnly.fa','fasta')
for l in x[16:35]:
	y=l.split()
	for i in y[1:]:
			try:
				d=out.write('>'+y[0]+'_'+i+'\n')
				d=out.write(str(recs[i].seq)+'\n')
			except KeyError:
				print(i)

from Bio import SeqIO
recs=SeqIO.index('Mtruncatula_285_Mt4.0v1.transcript_primaryTranscriptOnly.fa','fasta')
x=open('test.blast').readlines()
a={}
for l in x:
	gene=l.split()[0]
	gene=gene.split('_')[0]
	try:a[gene].append(l.split()[1])
	except KeyError:a[gene]=[l.split()[1]]
	
for k in a.keys():
	a[k]=list(set(a[k]))
	out=open(k+'.fas','a')
	for i in a[k]:
		d=SeqIO.write(recs[i],out,'fasta')

