




######################
# Input data
######################
# Input single cell expression data 
input_h5py_file="/n/groups/price/scdata/Perez_Science_2021/lupus.h5ad"
# pseudobulk h5ad file
input_h5py_pseudobulk_file="/n/groups/price/scdata/Perez_Science_2021/lupus-pseudobulk.h5ad"

# Gene annotation file (hg19)
gene_annotation_file="/n/groups/price/ben/reference_data/gene_annotation_files/gencode.v19.annotation.gff3"

# Directory containing genotype data
genotype_data_dir="/n/groups/price/scdata/Perez_Science_2021/genotype/plink_maf10/"

# individual file
individual_info_file="/n/groups/price/scdata/Perez_Science_2021/genotype/indv_info_from_scdata.tsv"



######################
# Output data
######################
# Scratch output root
scratch_output_root="/n/scratch3/users/b/bes710/eqtl_polygenecity/"
# regular output root
output_root="/n/groups/price/ben/eqtl_polygenecity/"
# Processed single cell expression dir
processed_sc_expression_dir=$scratch_output_root"processed_sc_expression/"
# Visualize processed single cell expression dir
visualize_processed_sc_expression_dir=$output_root"visualize_processed_sc_expression/"
# Processed single cell expression dir
processed_genotype_dir=$scratch_output_root"processed_genotype/"
# Processed single cell pseudobulk expression
pseudobulk_expression_dir=$output_root"pseudobulk_expression/"
# Processed single cell pseudobulk expression
polygenecity_analysis_input_dir=$scratch_output_root"polygenecity_analysis_input/"
# Processed single cell pseudobulk expression
polygenecity_analysis_results_dir=$output_root"polygenecity_analysis_results/"




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
sbatch process_sc_expression.sh $input_h5py_file $input_h5py_pseudobulk_file $processed_sc_expression_dir $individual_info_file $genotype_data_dir $visualize_processed_sc_expression_dir $processed_genotype_dir 
fi

# Process genotype
filtered_sample_info_file=$processed_sc_expression_dir"individual_info_european_rna_and_dna.txt"
if false; then
sh process_genotype.sh $processed_genotype_dir $genotype_data_dir $filtered_sample_info_file
fi

# Generate pseudobulk expression
if false; then
sbatch generate_pseudobulk_expression.sh $processed_sc_expression_dir $input_h5py_file $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file
fi




# Each line in this file is a data set. loop through it
data_set_summary_file=$pseudobulk_expression_dir"pseudobulk_data_set_summary_final.txt"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	echo $data_set_name
	sh prepare_pseudobulk_data_set_for_polygenecity_analysis.sh $data_set_name $pseudobulk_expression_file $plink_covariate_file $gene_info_file $plink_individual_file $polygenecity_analysis_input_dir $processed_genotype_dir
done
fi

# Each line in this file is a data set. loop through it
data_set_summary_file=$pseudobulk_expression_dir"cell_type_matched_pseudobulk_data_set_summary_final.txt"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	echo $data_set_name
	sbatch prepare_pseudobulk_data_set_for_polygenecity_analysis.sh $data_set_name $pseudobulk_expression_file $plink_covariate_file $gene_info_file $plink_individual_file $polygenecity_analysis_input_dir $processed_genotype_dir
done
fi





#####################
# eqtl polygenecity
#####################
data_set_summary_file=$pseudobulk_expression_dir"pseudobulk_data_set_summary_final.txt"
gene_heritability_normalization="True"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	data_set_gene_info_file=$polygenecity_analysis_input_dir$data_set_name"_organized_gene_info_non_negative_h2.txt"
	output_root=$polygenecity_analysis_results_dir$data_set_name"_polygenecity_results_non_negative_h2_gene_h2_norm_"$gene_heritability_normalization"_"
	sbatch run_polygenecity_analysis.sh $data_set_name $data_set_gene_info_file $output_root $gene_heritability_normalization
done
fi

data_set_summary_file=$pseudobulk_expression_dir"pseudobulk_data_set_summary_final.txt"
gene_heritability_normalization="True_estimate"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	data_set_gene_info_file=$polygenecity_analysis_input_dir$data_set_name"_organized_gene_info_non_negative_h2.txt"
	output_root=$polygenecity_analysis_results_dir$data_set_name"_polygenecity_results_non_negative_h2_gene_h2_norm_"$gene_heritability_normalization"_"
	sbatch run_polygenecity_analysis.sh $data_set_name $data_set_gene_info_file $output_root $gene_heritability_normalization
done
fi



data_set_summary_file=$pseudobulk_expression_dir"pseudobulk_data_set_summary_final.txt"
gene_heritability_normalization="False"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	data_set_gene_info_file=$polygenecity_analysis_input_dir$data_set_name"_organized_gene_info_non_negative_h2.txt"
	output_root=$polygenecity_analysis_results_dir$data_set_name"_polygenecity_results_non_negative_h2_gene_h2_norm_"$gene_heritability_normalization"_"
	sbatch run_polygenecity_analysis.sh $data_set_name $data_set_gene_info_file $output_root $gene_heritability_normalization
done
fi

#####################
# eqtl polygenecity (imputed version)
#####################
gene_heritability_normalization="False"
if false; then
sed 1d $data_set_summary_file | while read data_set_name cluster_method_name cluster_name pseudobulk_expression_file covariate_file num_donors num_genes num_cells_per_individual_file gene_info_file plink_individual_file plink_covariate_file; do
	data_set_gene_info_file=$polygenecity_analysis_input_dir$data_set_name"_organized_gene_info_non_negative_h2.txt"
	output_root=$polygenecity_analysis_results_dir$data_set_name"_imputed_polygenecity_results_non_negative_h2_gene_h2_norm_"$gene_heritability_normalization"_"
	sbatch run_imputed_polygenecity_analysis.sh $data_set_name $data_set_gene_info_file $output_root $gene_heritability_normalization
done
fi






























#####################
# OLD
#####################



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

if false; then
module load R/3.5.1
Rscript visualize_fusion_results.R $data_set_summary_file $fusion_processed_output_dir $pseudobulk_expression_dir $fusion_visualization_dir
fi


