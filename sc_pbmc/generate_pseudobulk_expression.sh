#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --nodes=1




processed_sc_expression_dir="$1"
processed_genotype_dir="$2"
pseudobulk_expression_dir="$3"
gene_annotation_file="$4"
plink_input_dir="$5"



if false; then
python generate_pseudobulk_expression.py $processed_sc_expression_dir $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file
fi

python get_pseudobulk_data_in_plink_format.py $pseudobulk_expression_dir $gene_annotation_file $plink_input_dir
