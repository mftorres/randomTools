import re
import os
import glob
import sys

aln_paths='./alns/'
fileSuf='*.fas'

sample_names=[]
genes=[]
gene_seqfiles_dict={}
notmissing_fromfile={}
geneseq_length={}

for filename in glob.glob(os.path.join(aln_paths,'%s'%(fileSuf))): # files have alignments
    gene=filename.split('\\')[-1]
    genes.append(gene.strip('.fas'))
    # binds sequences and samples
    sampleseq_infile={}
    # tracks missing samples in the file
    sampleinfile_list=[]
    for line in open(filename,'r'):
        matchsample=re.search('(>[A-Za-z0-9]+)',line)
        matchseq=re.search('^([AGTCNagtcn]+)',line)
        if matchsample:
            sample=matchsample.group(1)
            if sample not in sample_names:
                sample_names.append(sample)
            if sample not in sampleinfile_list:
                sampleinfile_list.append(sample)
        if matchseq:
            seq=matchseq.group(1)
            seq_lenght=len(seq)
        sampleseq_infile[sample]=seq
    gene_seqfiles_dict[gene.strip('.fas')]=sampleseq_infile
    notmissing_fromfile[gene.strip('.fas')]=sampleinfile_list
    geneseq_length[gene.strip('.fas')]=seq_lenght

wdummies_path='./wdummies/'
for gene in genes:
    new_fasta=open(os.path.join(wdummies_path,'%s_wd.fasta'%(gene)),'a')
    for sample in sample_names:
        if sample in gene_seqfiles_dict[gene].keys():
            new_fasta.write('%s\n%s\n'%(sample,gene_seqfiles_dict[gene][sample]))
        elif sample not in gene_seqfiles_dict[gene].keys():
            new_fasta.write('%s\n%s\n'%(sample,'N'*geneseq_length[gene]))
    new_fasta.close()
    print('Done')
