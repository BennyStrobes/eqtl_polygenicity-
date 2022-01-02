#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --mem=200GB
#SBATCH --partition=shared



input_h5py_file="$1"
input_h5py_scran_file="$2"
processed_sc_expression_dir="$3"
visualize_processed_sc_expression_dir="$4"
processed_genotype_dir="$5"


module load R/3.6.1
module load python/3.7.4-anaconda
if false; then
python process_sc_expression.py $input_h5py_file $input_h5py_scran_file $processed_sc_expression_dir
fi

if false; then
python generate_single_cell_clusters_at_various_resolutions.py $processed_sc_expression_dir
fi

if false; then
python make_list_of_individuals_with_rna_seq_and_of_european_ancestry.py $processed_sc_expression_dir $processed_genotype_dir
fi


module load R/3.5.1
Rscript visualize_processed_sc_expression.R $processed_sc_expression_dir $visualize_processed_sc_expression_dir

# Make sure same cells/same order in _scran_normalized_all_genes.h5ad vs _2.h5ad
# Come up with some basic visualizations of SC expression (ie umap expression colored by cell type, etc)
# Create matrix of dimension cellsXclutering methods where element represents cluster assignment of cell according to that method
# Visualize 