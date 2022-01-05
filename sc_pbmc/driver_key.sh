




######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/expression/Lupus_study_adjusted.h5ad"
# Scran normalized h5ad file (generated here: https://github.com/BennyStrobes/qtl_factorization/blob/main/ye_lab_single_cell/preprocess_scran_single_cell_expression.py)
input_h5py_scran_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/processed_expression/scran_normalization.h5ad"

# Gene annotation file (hg19)
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Directory containing genotype data
genotype_data_dir="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/"

# VCF file of genotype
# Previously generated here: https://github.com/BennyStrobes/qtl_factorization/blob/main/ye_lab_single_cell/preprocess_genotype.sh
genotype_input_vcf_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/processed_genotype/clues_immvar_merged.vcf.gz"

# File containing mapping from old genotype ids to new genotype ids
# Basically rna-seq had old genotype ids and genotype data had new genotype ids
# We will convert rna-seq to new genotype ids
genotype_id_mapping_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/sample_mappings_ye_SLE.txt"

# File containing genotyped individuals
genotyped_individuals_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/genotyped_individuals.txt"

# Gene annotation file (hg19)
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Directory containing fusion code
fusion_code_dir="/work-zfs/abattle4/bstrober/tools/fusion_twas-master/"

# LDREF directory (downloaded from https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2 on Dec. 13 2021)
LDREF_dir="/work-zfs/abattle4/bstrober/eqtl_polygenicity/input_data/LDREF/"

######################
# Output data
######################
# Root output directory
output_root="/work-zfs/abattle4/bstrober/eqtl_polygenicity/sc_pbmc/"
# Processed single cell expression dir
processed_sc_expression_dir=$output_root"processed_sc_expression/"
# Visualize processed single cell expression dir
visualize_processed_sc_expression_dir=$output_root"visualize_processed_sc_expression/"
# Processed single cell expression dir
processed_genotype_dir=$output_root"processed_genotype/"
# Processed single cell pseudobulk expression
pseudobulk_expression_dir=$output_root"pseudobulk_expression/"
# PLINK input dir
one_KG_snps_dir=$output_root"one_KG_snps/"
# PLINK input dir
plink_input_dir=$output_root"plink_input/"
# FUSION INPut data dir
fusion_input_dir=$output_root"fusion_input_1kg/"
# Fusion results dir
fusion_output_dir=$output_root"fusion_output_1kg/"
# Fusion results dir
fusion_processed_output_dir=$output_root"fusion_processed_ouput/"
# Visualize fusion results
fusion_visualization_dir=$output_root"fusion_visualization/"


# Process SC expression
if false; then
sh process_sc_expression.sh $input_h5py_file $input_h5py_scran_file $processed_sc_expression_dir $visualize_processed_sc_expression_dir $processed_genotype_dir
fi

# Process genotype
if false; then
sh process_genotype.sh $processed_genotype_dir $genotype_input_vcf_file
fi

# Generate pseudobulk expression
if false; then
sh generate_pseudobulk_expression.sh $processed_sc_expression_dir $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file $plink_input_dir
fi

# Extract 1KG snps
if false; then
sh extract_1kg_snps.sh $LDREF_dir $one_KG_snps_dir
fi

# Each line in this file is a data set. loop through it
data_set_summary_file=$plink_input_dir"pseudobulk_data_set_summary_plink_ready.txt"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file plink_gene_summary_file; do
	sbatch run_fusion_on_pseudobulk_data_set.sh $data_set_name $plink_gene_summary_file $fusion_input_dir $fusion_output_dir $processed_genotype_dir $fusion_code_dir $one_KG_snps_dir
done
fi





# Organize fusion output
if false; then
sh organize_fusion_results.sh $data_set_summary_file $fusion_output_dir $fusion_processed_output_dir
fi

module load R/3.5.1
Rscript visualize_fusion_results.R $data_set_summary_file $fusion_processed_output_dir $pseudobulk_expression_dir $fusion_visualization_dir



