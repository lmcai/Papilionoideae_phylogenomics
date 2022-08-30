from Bio import SeqIO
recs=SeqIO.index('Lj2.5_cds.ffn','fasta')
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

