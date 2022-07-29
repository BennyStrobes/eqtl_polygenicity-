#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-60:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



data_set_name="$1"
gene_info_file="$2"
output_root="$3"
gene_heritability_normalization="$4"


python3 run_imputed_polygenecity_analysis.py $data_set_name $gene_info_file $output_root $gene_heritability_normalization