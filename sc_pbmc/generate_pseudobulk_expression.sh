#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




processed_sc_expression_dir="$1"
input_h5py_file="$2"
processed_genotype_dir="$3"
pseudobulk_expression_dir="$4"
gene_annotation_file="$5"

if false; then
python3 generate_pseudobulk_expression.py $processed_sc_expression_dir $input_h5py_file $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file
fi



if false; then
python3 generate_cell_type_matched_pseudobulk_expression.py $processed_sc_expression_dir $input_h5py_file $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file
fi





python3 generate_gene_info_file_and_ordered_individual_file.py $pseudobulk_expression_dir $gene_annotation_file























#############
# OLD
################




if false; then
python get_pseudobulk_data_in_plink_format.py $pseudobulk_expression_dir $gene_annotation_file $plink_input_dir
fi