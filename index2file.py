#######################################################################################################################################
# This script recovers the index information from the fastq headers and prints a configuration file compatible with SECAPR
# See SECAPR here: https://github.com/AntonelliLab/seqcap_processor
# This script is written by Maria Fernanda Torres
# Copyright 2020 Maria Ferndanda Torres
# Licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
#######################################################################################################################################

import re
import os
import glob
import sys
import time
from itertools import islice # read range of lines after a coditional is True

# current working directory
print('Working directory: ',os.getcwd())

input_path=os.getcwd()
output_path=os.getcwd()
fileSuf='*fastq'
outputfile='%s/config_file.txt'%(os.getcwd())
print('Path to output: ',outputfile)


adapter7i='i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG'
adapter5isingle='i5:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
adapter5idouble='i5:AATGATACGGCGACCACCGAGATCTACAC*ACACTCTTTCCCTACACGACGCTCTTCCGATCT'

# asigns all files in input_path with fileSuf suffix:
file_configfile=open(outputfile,'w')
file_configfile.write('[adapters]\n')
file_configfile.write('%s\n'%(adapter7i))

adapterswitch='single'
dict_sample_ID={}
for filename in glob.glob(os.path.join(input_path,'%s'%(fileSuf))):
    samplefile=filename.split('/')[-1]
    tempsample=re.search('(.+).(R..fastq)',samplefile)
    sampleID=tempsample.group(1) if tempsample and tempsample not in dict_sample_ID.keys() else None
    indices=[]
    with open(filename) as f:
        for line in islice(f,0,140,4):
            index7i=''
            index5i=''
            temp=re.search('\s\d:.:\d+:(\w+.\w+)',line)
            if '+' in temp.group(1):
                adapterswitch='double'
                temp1=re.search('\s\d:.:\d+:(\w+).(\w+)',line)
                if 'N' not in temp1.group(1):
                # put them on dict with file name, print at the end.
                    index7i=temp1.group(1)
                if 'N' not in temp1.group(2):
                    index5i=temp1.group(2)
            else:
                adapterswitch='single'
                temp1=re.search('\s\d:.:\d+:(\w+)',line)
                if 'N' not in temp1.group(1):
                # put them on dict with file name, print at the end.
                    index7i=temp1.group(1)
                    index5i='nan'
    indices.append([index7i,index5i])
    dict_sample_ID[sampleID]=indices
file_configfile.write('%s\n'%(adapter5isingle)) if adapterswitch == 'single' else file_configfile.write('%s\n'%(adapter5idouble))

file_configfile.write('\n[names]\n')
for key in dict_sample_ID.keys():
    file_configfile.write('%s:_\n'%(key))

file_configfile.write('\n[barcodes]\n')
if adapterswitch == 'double':
    for key,value in dict_sample_ID.items():
        for i in value:
            file_configfile.write('i7-%s:%s\n'%(key,i[0]))
    for key,value in dict_sample_ID.items():
        for i in value:
            file_configfile.write('i5-%s:%s\n'%(key,i[1]))
else:
    for key,value in dict_sample_ID.items():
        for i in value:
            file_configfile.write('i7-%s:%s\n'%(key,i[0]))

file_configfile.close()

print('Done')
