#!/bin/bash

#$ -cwd
#$ -N scanpy_to_seurat
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=100G

Rscript scanpy_to_seurat.R $1

echo "End on `date`"
