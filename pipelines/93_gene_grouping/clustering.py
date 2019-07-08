#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:52:21 2019

@author: doru
"""

import sys
from os.path import join
import pandas as pd
import numpy as np

args = sys.argv
output_folder = args[1]
no_clusters = int(args[2])

expression_file = join(output_folder, "expression.csv")
expression_df = pd.read_csv(expression_file, index_col = 0)
expression      = np.transpose(expression_df.values)

from sklearn.mixture import GaussianMixture
clustering = GaussianMixture(n_components = no_clusters, random_state = 19).fit(expression)
clustering = clustering.predict(expression)

# save the output
gene_names = list(expression_df.head(0))
df = {"GeneNames": gene_names, "Cluster": clustering}
df = pd.DataFrame.from_dict(df)

df.to_csv(join(output_folder, "clustering.csv"))
