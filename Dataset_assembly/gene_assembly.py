from Bio import SeqIO
import pandas as pd
import sys

def filter_rec_per_query(blast_results):
	pergene_rec={}
	filtered_hits=[]
	#sort blast hits by gene
	for l in blast_results:
		try:pergene_rec[l.split()[0]].append(l)
		except KeyError:pergene_rec[l.split()[0]]=[l]
	#get one best rec per gene
	for g in pergene_rec.keys():
		blast_lst=[l.strip().split() for l in pergene_rec[g]]
		blast_table= pd.DataFrame(blast_lst, columns =['V'+str(i) for i in list(range(12))]) 
		#sort dataframe
		#blast_table["V2"] = pd.to_numeric(blast_table["V2"])
		#blast_table=blast_table.sort_values(by='V2', ascending=False)
		blast_table['ref_start'] = [0] *len(blast_table)
		blast_table['ref_end'] = [0] *len(blast_table)
		for i in range(0,len(blast_table)):
			blast_table.iloc[i,12]=min(int(blast_table.iloc[i,8]),int(blast_table.iloc[i,9]))
			blast_table.iloc[i,13]=max(int(blast_table.iloc[i,8]),int(blast_table.iloc[i,9]))
		#blast_table=blast_table.sort_values(by='ref_start')
		blast_table["V6"] = pd.to_numeric(blast_table["V6"])
		blast_table["V7"] = pd.to_numeric(blast_table["V7"])
		blast_table=blast_table.sort_values(by='V6')
		cur_rec=blast_table.iloc[0]
		#consolidate ranges and preserve multiple hits on one sequences
		for i in range(0,len(blast_table)):
			if blast_table.iloc[i,6]<=cur_rec.V7:
				#extend range
				if blast_table.iloc[i,7]>cur_rec.V7:
					cur_rec.V7=blast_table.iloc[i,7]
					if blast_table.iloc[i,13]>cur_rec.ref_end:
						cur_rec.ref_end=blast_table.iloc[i,13]
			else:
				#different region
				cur_rec.V6=str(cur_rec.V6); cur_rec.V7=str(cur_rec.V7); cur_rec.ref_start=str(cur_rec.ref_start); cur_rec.ref_end=str(cur_rec.ref_end)
				filtered_hits.append('\t'.join(cur_rec))
				cur_rec=blast_table.iloc[i]
		#append the last record
		cur_rec.V6=str(cur_rec.V6); cur_rec.V7=str(cur_rec.V7); cur_rec.ref_start=str(cur_rec.ref_start); cur_rec.ref_end=str(cur_rec.ref_end)
		filtered_hits.append('\t'.join(cur_rec))
	#return list of strings as best hits
	return(filtered_hits)


def gene_assem(blast_results,assemb,species_name):
	#initiate values
	cur_gene=blast_results[0].split('\t')[1].split('_')[0]
	if '_old' in blast_results[0]:cur_gene=cur_gene+'_old'
	blast_syntax={}
	blast_syntax[cur_gene]=[blast_results[0]]
	for l in blast_results:
		gene=l.split('\t')[1].split('_')[0]
		if '_old' in l:
			gene=gene+'_old'
		if gene==cur_gene:
			blast_syntax[gene].append(l)
		else:
			#the blast syntax for last gene is complete, process them to extract sequences
			hits=filter_rec_per_query(blast_syntax[cur_gene])
			hits=[l.split() for l in hits]
			hit_table=pd.DataFrame(hits, columns =['V'+str(i) for i in list(range(14))]) 
			#hit_table["V2"] = pd.to_numeric(hit_table["V2"])
			hit_table["V10"] = pd.to_numeric(hit_table["V10"])
			hit_table['V12'] = pd.to_numeric(hit_table["V12"])
			hit_table['V13'] = pd.to_numeric(hit_table["V13"])
			hit_table=hit_table.sort_values(by='V12')
			if len(hit_table)==1:
			#if only one rec, output it
				if int(hit_table.iloc[0,8])>int(hit_table.iloc[0,9]):
					seq_str = str(assemb[hit_table.iloc[0,0]].seq[(int(hit_table.iloc[0,6])-1):int(hit_table.iloc[0,7])].reverse_complement())
				else:
					seq_str = str(assemb[hit_table.iloc[0,0]].seq[(int(hit_table.iloc[0,6])-1):int(hit_table.iloc[0,7])])
			#print(hit_table)
			#examine if the boundaries of exon hits are overlapping
			else:
			#multiple regions
				cur_region = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
				cur_region = cur_region.append(hit_table.iloc[0])
				right_end = cur_region.iloc[0,13]
				output_recs = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
				for i in range(0,len(hit_table)):
					if hit_table.iloc[i,12]-right_end>-10:
						#move to the next block, get optimum hit for the current region
						cur_region = cur_region.sort_values(by='V10')
						#print(cur_region)
						output_recs = output_recs.append(cur_region.iloc[0])
						cur_region = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
						cur_region = cur_region.append(hit_table.iloc[i])
						right_end = cur_region.iloc[0,13]
					else:
						#left end overlap
						cur_region = cur_region.append(hit_table.iloc[i])
						if hit_table.iloc[i,13]>right_end:right_end=hit_table.iloc[i,13]
				#add the last rec
				cur_region = cur_region.sort_values(by='V10')
				output_recs = output_recs.append(cur_region.iloc[0])
				print(output_recs)
				seq_str=''
				for i in range(0,len(output_recs)):
					#reverse compliment the sequence if needed
					if int(output_recs.iloc[i,8])>int(output_recs.iloc[i,9]):
						seq_str=seq_str+'NNNNN'+str(assemb[output_recs.iloc[i,0]].seq[(int(output_recs.iloc[i,6])-1):int(output_recs.iloc[i,7])].reverse_complement())
					else:
						seq_str=seq_str+'NNNNN'+str(assemb[output_recs.iloc[i,0]].seq[(int(output_recs.iloc[i,6])-1):int(output_recs.iloc[i,7])])
			#write sequence to output file
			out=open(cur_gene+'.fas','a')
			out.write('>'+species_name+'\n'+seq_str+'\n')
			out.close()
			#update the name of the working gene
			cur_gene=gene
			blast_syntax[gene]=[l]


x=open(sys.argv[1]).readlines()
x=x+['end\tend\t86.111\t324\t37\t2\t60\t383\t1315\t1000\t4.70e-106\t379\n']
assemb=SeqIO.index('scaffolds.fasta','fasta')
gene_assem(x,assemb,sys.argv[2])

##########################
#			ALL_exons=1
#			for i in range(1,len(hit_table)):
#				if int(hit_table.iloc[i,12])-int(hit_table.iloc[i-1,13])<-10:ALL_exons=0
#			#output to sequences
#			if ALL_exons:
#				#write multiple exons in one sequence
#				seq_str=''
#				for i in range(0,len(hit_table)):
#					#reverse compliment the sequence if needed
#					if int(hit_table.iloc[i,8])>int(hit_table.iloc[i,9]):
#						seq_str=seq_str+'NNNNN'+str(assemb[hit_table.iloc[i,0]].seq[(int(hit_table.iloc[i,6])-1):int(hit_table.iloc[i,7])].reverse_complement())
#					else:
#						seq_str=seq_str+'NNNNN'+str(assemb[hit_table.iloc[i,0]].seq[(int(hit_table.iloc[i,6])-1):int(hit_table.iloc[i,7])])
#			else:
#				#overlapping blast result, difficult to determine synteny, only output best hit
#				print('synteny error: '+cur_gene)
#				hit_table=hit_table.sort_values(by='V10')
#				#print(hit_table)
#				seq_str = str(assemb[hit_table.iloc[0,0]].seq)
