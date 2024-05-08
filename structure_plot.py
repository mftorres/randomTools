# nano structure_plot.py
# run in angsd environment
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import time
import glob
import os
from os import path
from matplotlib import gridspec
import argparse


# you have this file, from previous lectures
samples = {}
with open('./Heliconius_localities.csv') as file:
    next(file)
    for line in file:
        # samples[line.split('\t')[0]] = line.split('\t')[28]
        samples[line.split(',')[0]] = line.split(',')[9]

# print('samples dict')
# print(samples)

sample_VCForder = {}
# temp.ind is a list of samples, split by a dot, no spaces
# in the order they would be plotted, left to right

with open('./temp.ind') as file:
    for line in file:
        temp = line.strip('\n').split('.')
#         print(temp)
        sample_VCForder = dict(zip(range(0,len(temp)), temp))

def plot_structure_results():
    cwd = os.getcwd()
    ext ='%s*meanQ'%(args.file_prefix)
    max_k = args.max_k
    print('Max k: ',max_k)
    fig, ax = plt.subplots(figsize=(10, 35), facecolor='w')
    G = gridspec.GridSpec(max_k, 1, hspace=0.2, wspace=1)

    cmap = mpl.cm.tab20c
    pop_list = ['k_%s'%(x+1) for x in range(0,max_k,1)]
    #col_list = [cmap(float(pop)/(max_k-1)) for pop in range(0,max_k,1)]
    col_list = [cmap(float(pop)/(max_k)) for pop in range(0,max_k,1)]
    pop_col_dict = dict(zip(pop_list,col_list))
    print(pop_col_dict)

    for filename in glob.glob(os.path.join(cwd, ext)):
#        print('---> Filename: ',filename)
        knum = int(filename.split('.')[-2])
#        print('---> knum: ',knum)
        names = ['k_%s'%(x+1) for x in range(0,knum,1)]
#         print('names poplist',names, pop_list)
        temp = pd.read_csv(filename, sep = '\s+', header = None, index_col = False,
                              names = names)
        results = temp.sort_values(by = list(temp.columns), ascending = False)
        results.reset_index(inplace = True)

        print(results.head(2))

        ax = plt.subplot(G[knum-1, 0], facecolor='w')
        x = range(0,len(results),1)
        bottom = 0
        for i in range(0,knum,1):
            width = 1
            height = results['k_%s'%(i+1)]
            ax.bar(x, height, width, bottom, edgecolor = 'w', linewidth = 1, color = pop_col_dict['k_%s'%(i+1)])
            bottom = bottom + results['k_%s'%(i+1)]
            plt.ylabel('k_%s'%(i+1))

        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_xlim(-1, len(results)+1)
        ax.set_ylim(0, 1)
        [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'bottom']]

    for key in sample_VCForder.keys():
        ax.text(x = key, y = 0, s = samples[sample_VCForder[key]],
        rotation = 90, va = 'top', fontsize=15)


    plt.savefig("{}_plot_test.pdf".format(args.file_prefix))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Python3 code for producing population structure plots")
    parser.add_argument("file_prefix", type = str, help="File prefix for all K.meanQ files with path, e.g. path/msp_0.31_fs.out, string")
    parser.add_argument("max_k", type = int, help="Maximum number of K used for structure analyses, e.g. 4, integer")

    args = parser.parse_args()

    plot_structure_results()

# run as python3 structure_plot.py file_prefix max_k
