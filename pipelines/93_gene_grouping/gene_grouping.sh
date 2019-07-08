#!/bin/bash

#$ -cwd
#$ -N gene_grouping
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=100G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript gene_grouping.R $1

echo "End on `date`"