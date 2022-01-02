#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --mem=5GB
#SBATCH --partition=shared





LDREF_dir="$1"
one_KG_snps_dir="$2"


python extract_1kg_snps.py $LDREF_dir $one_KG_snps_dir