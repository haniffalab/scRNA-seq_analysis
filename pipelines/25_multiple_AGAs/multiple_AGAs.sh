#!/bin/bash

#$ -cwd
#$ -N multiple_AGAs
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=100G

Rscript multiple_AGAs.R $1

echo "End on `date`"
