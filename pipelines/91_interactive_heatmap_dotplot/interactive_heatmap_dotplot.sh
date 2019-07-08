#!/bin/bash

#$ -cwd
#$ -N interactive_heatmap_dotplot
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=100G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript interactive_heatmap_dotplot.R $1

echo "End on `date`"