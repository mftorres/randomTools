%%bash
perl -p -i -e 's/\n/@@/g' extracted_target_contigs_all_samples.fasta;
perl -p -i -e 's/@@>/\n>/g' extracted_target_contigs_all_samples.fasta;
perl -p -i -e 's/([0-9]+)@@/\1\n/g' extracted_target_contigs_all_samples.fasta;
perl -p -i -e 's/@@//g' extracted_target_contigs_all_samples.fasta;

%%python
import re
import os
import glob
import sys
import time

print(os.getcwd())
path=os.getcwd()
file='extracted_target_contigs_all_samples.fasta'

Loci=[]
Sample=[]
Sequence=[]
for line in open(os.path.join(path,'%s'%(file)),'r'):
    # header
    if '>' in line:
        tempsearch=re.search('>([0-9]+)_(.+)\s.+',line)
        Loci.append(tempsearch.group(1) if tempsearch else None)
        Sample.append(tempsearch.group(2) if tempsearch else None)
    else:
        Sequence.append(line)
for i in range(0,len(Loci),1):
    file_fasta=open(os.path.join(path,'Loci_%s_extracted.fasta'%(Loci[i])),'a')
    file_fasta.write('>%s\n%s'%(Sample[i],Sequence[i]))
    file_fasta.close()
print('Done')
