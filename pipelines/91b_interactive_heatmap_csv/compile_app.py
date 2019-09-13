#!/usr/bin/env python3
import sys
args = sys.argv
options_file = args[1]

with open(options_file, "r") as options:
    options_fields = options.readlines()
    csv_filename   = options_fields[0].strip()
    output_dir     = options_fields[1].strip()
    save_to        = options_fields[2].strip()
    data_name      = options_fields[3].strip()

# function to truncated floats to 3 digits - for memory efficiency
def truncateFloat(val):
    return round(val, 2)

# join save_to to output_dir and create file name for expression data csv file
from os.path import join
save_to = join(output_dir, save_to)

# open the required csv file
import pandas as pd
expression_data = pd.read_csv(csv_filename, index_col = 0, header = 0)
expression_data = expression_data.apply(truncateFloat)

# populate data
gene_names = expression_data.columns.values
cell_names = expression_data.index.values
data = []
for gene_name in gene_names:
    indata = expression_data[gene_name][cell_names].values
    indata = [str(d) for d in indata]
    indata = ','.join(indata)
    indata = 'expression_data["{gene_name}"] = [{values}]'.format(gene_name = gene_name, values = indata)
    data.append(indata)
data = '\n'.join(data)
cell_names = ['"{cell_name}"'.format(cell_name = cell_name) for cell_name in cell_names]
cell_names = ','.join(cell_names)
cell_names = 'cell_names = [{values}]'.format(values = cell_names)

data_name_var = "dataset_name = '{dataname}'".format(dataname = data_name)
data = '\n'.join([data, cell_names, data_name_var])

# open template and insert data
template_addr = 'template.html'
template_fobj = open(template_addr, 'r')
template = template_fobj.read()
template_fobj.close()

with open(save_to, 'w') as save_to_fobj:
    template = template.replace('// insert data here', data)
    save_to_fobj.write(template)
