#!/bin/bash

#$ -cwd
#$ -N violin_plots
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=200G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript violin_plots.R $1

echo "End on `date`"
