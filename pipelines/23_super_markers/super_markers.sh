#!/bin/bash

#$ -cwd
#$ -N super_markers
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=400G

Rscript super_markers.R $1

echo "End on `date`"
