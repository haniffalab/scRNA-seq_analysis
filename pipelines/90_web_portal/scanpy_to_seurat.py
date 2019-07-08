#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 12:00:47 2018

@author: doru
"""

# gene expression should be at data.X
# gene symbols should be at data.var.GeneName
# DRs should be at data.obsm
# categoeis should be columns at data.obs

import sys
option_files = sys.argv[1]

#option_files = "thymus_web_portal_options.txt"

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc
from os.path import join
from scipy.io import mmwrite
from scipy.sparse import hstack
import numpy as np
import random
import pandas as pd

# extract options
options_fobj = open(option_files, 'r')
options      = options_fobj.readlines()
options_fobj.close()
options      = [op.strip() for op in options]
# make file_name variable
file_name    = options[0].split(':')[1].strip()
# make output_folder variable
output_folder = options[1].split(':')[1].strip()
# make dr_coordinates
dr_coordinates = options[2].split(':')[1].strip()
dr_coordinates = dr_coordinates.split(";")
dr_coordinates = [dr.split('->')[1] for dr in dr_coordinates if dr != '']
# extract categories
categories = options[3].split(';')
categories = [cat.split('->') for cat in categories if cat != '']

data = sc.read(file_name)

# save expression data
expression_data_filename = join(output_folder, 'expression')
nUMI = np.array(data.obs.n_UMIs)
nUMI = nUMI.reshape((nUMI.shape[0], 1))
nGene = np.array(data.obs.n_genes)
nGene = nUMI.reshape((nGene.shape[0], 1))
gene_expression = hstack([data.X, nUMI, nGene])
mmwrite(expression_data_filename, gene_expression)

# save gene names
gene_names = [g for g in data.var.GeneName]
gene_names.append('nUMI')
gene_names.append('nGene')
gene_names = '\n'.join(gene_names)
gene_names_fobj = open(join(output_folder, 'gene_names.txt'), 'w')
gene_names_fobj.write(gene_names)
gene_names_fobj.close()

# get DR coordinates
for index, dr_coor in enumerate(dr_coordinates):
    if index == 0:
        dr_data = data.obsm[dr_coor]
    else:
        dr_data = np.concatenate([dr_data, data.obsm[dr_coor]], axis = 1)
np.savetxt(join(output_folder, "dr.csv"), dr_data, delimiter=",")
    
# get categories
for index, category in enumerate(categories):
    slot          = category[1]
    category_col  = category[2] 
    category_data = data.obs[slot]
    if category_col == 'null':
        category_unique = category_data.values.categories
        r = lambda: random.randint(0,255)
        category_col = {}
        for unique in category_unique:
            category_col[unique] = '#%02X%02X%02X' % (r(),r(),r())
        category_col = category_data.map(category_col)
    else:
        category_col = pd.read.csv('when colour keys are ready continue script here')
    categories_data = pd.concat([category_data, category_col], axis = 1)
    categories_data.columns = [category[0], '{c}_Colour'.format(c = category[0])]
    if index == 0:
        csv_data = categories_data
    else:
        csv_data = pd.concat([csv_data, categories_data], axis = 1)
categories_fname = join(output_folder, 'categories.csv')
csv_data.to_csv(categories_fname)
