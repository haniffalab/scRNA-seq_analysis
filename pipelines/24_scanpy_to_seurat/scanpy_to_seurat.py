#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 12:00:47 2018

@author: doru
"""

# gene expression should be at data.X
# gene names should be at data.var_names
# cell names should be in data.obs_names and must be unique
# categories (the meta data) should be columns at data.obs
# any exception from the above assumptions will result in errors

import sys
file_name     = sys.argv[1]
output_folder = sys.argv[2]
print("printing file_name")
print(file_name)
print("printing output folder")
print(output_folder)

import matplotlib; matplotlib.use('Agg');
import scanpy.api as sc
from os.path import join
from scipy.io import mmwrite

data = sc.read(file_name)
print("printing data")
print(data)

# save expression data
expression_data_filename = join(output_folder, 'expression')
gene_expression = data.X
mmwrite(expression_data_filename, gene_expression)

# save gene names
gene_names_file = join(output_folder, "gene_names.csv")
data.var_names.to_series().to_csv(gene_names_file)

# save cell names
cell_names_file = join(output_folder, "cell_names.csv")
data.obs_names.to_series().to_csv(cell_names_file)
 
# get categories
meta_data = join(output_folder, "meta_data.csv")
data.obs.to_csv(meta_data)

