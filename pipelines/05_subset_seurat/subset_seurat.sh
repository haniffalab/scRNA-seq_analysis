#!/bin/bash

#$ -cwd
#$ -N subset_seurat
#$ -V
#$ -l h_rt=47:59:59
#$ -l h_vmem=300G

Rscript subset_seurat.R

echo "End on `date`"
