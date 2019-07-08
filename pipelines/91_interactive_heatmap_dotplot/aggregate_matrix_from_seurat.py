#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 12:00:47 2018

@author: doru
"""

# gene expression should be at data.X
# DRs should be at data.obsm
# categoeis should be columns at data.obs

import sys
file_name     = sys.argv[1]
category      = sys.argv[2]
output_folder = sys.argv[3]

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc
from os.path import join
from scipy.io import mmwrite

data = sc.read(file_name)

# save expression data
expression_data_filename = join(output_folder, 'expression')
gene_expression = data.X
mmwrite(expression_data_filename, gene_expression)

# save gene names
gene_names = "\n".join(data.var_names.tolist())
with open(join(output_folder, 'gene_names.txt'), "w") as gene_names_fobj:
    gene_names_fobj.write(gene_names)
    
# get categories
categories = "\n".join(data.obs[category].tolist())
with open(join(output_folder, "cell_types.txt"), "w") as categories_file:
    categories_file.write(categories)
