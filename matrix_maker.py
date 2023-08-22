from Bio.Seq import Seq
from Bio import SeqIO
import time

# we know which number corresponds to which loci and genome based on the annotations using blast against GenBank data
mtDNA=[0,1,2,3,4,7,10,11,14,17,18,24,48,49,62,71,92,174]
nDNA=[19,21,28,39,42,46,50,57,59,96,97,98,134,138,139,140,159,160,170,171,176,177,178,179,190,202,203,216,217,218,219]

sample_list=[]

# iterates through all fastas and stores all unique sample IDs
for seq in mtDNA:
	for seq_record in SeqIO.parse("./mito/sequences_%s.fasta"%(seq), "fasta"):
		if seq_record.id not in sample_list:
			sample_list.append(seq_record.id)
		else:
			pass

for seq in nDNA:
	for seq_record in SeqIO.parse("./nuc/sequences_%s.fasta"%(seq), "fasta"):
		if seq_record.id not in sample_list:
			sample_list.append(seq_record.id)
		else:
			pass

# iterates through all loci files and creates a dictionary of sample=list of sequences throughout all loci
# the list of sequences are ordered equal
spp_seq_mit={}
for seq in mtDNA:
	seq_list=[]
	for spp in sample_list:
		for seq_record in SeqIO.parse("./mito/sequences_%s.fasta"%(seq), "fasta"):
			if spp == seq_record.id:
				seq_list.append(seq_record.seq)
		spp_seq_mit[spp]=seq_list

spp_seq_nuc={}
for seq in nDNA:
	seq_list=[]
	for spp in sample_list:
		for seq_record in SeqIO.parse("./nuc/sequences_%s.fasta"%(seq), "fasta"):
			if spp == seq_record.id:
				seq_list.append(seq_record.seq)
		spp_seq_nuc[spp]=seq_list

# testing print
print(spp_seq_mit['Ctenomys_goodfellowi'])
print(spp_seq_nuc['Ctenomys_goodfellowi'])

# opens mito matrix file
mtDNA_f = open('./supermatrix_mt.fasta','a+')
nDNA_f = open('./supermatrix_nuc.fasta','a+')

#spp_mit_concat={}
for key,value in spp_seq_mit.items():
	#concat=Seq("")
	concat=''
	for seq in value:
		concat += seq
	#spp_mit_concat[key]=concat
	mtDNA_f.write('>%s\n'%(key))
	mtDNA_f.write('%s\n'%(concat))

for key,value in spp_seq_nuc.items():
	concat=''
	for seq in value:
		concat += seq
	nDNA_f.write('>%s\n'%(key))
	nDNA_f.write('%s\n'%(concat))

mtDNA_f.close()
nDNA_f.close()

#print(spp_mit_concat['Ctenomys_goodfellowi'])
#print(len(spp_mit_concat['Ctenomys_goodfellowi']))
