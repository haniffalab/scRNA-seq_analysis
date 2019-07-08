#!/bin/bash

#$ -cwd
#$ -N seurat_to_interactive_gene_expression
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript seurat_to_interactive_gene_expression.R $1

echo "End on `date`"
