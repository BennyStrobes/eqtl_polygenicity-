#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --mem=5GB
#SBATCH --partition=shared



data_set_name="$1"
data_set_gene_summary_file="$2"
fusion_input_dir="$3"
fusion_output_dir="$4"
processed_genotype_dir="$5"
fusion_code_dir="$6"
LDREF_dir="$7"

PRE_GENO=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_plink_chr_"

module load R/3.5.1

cis_window_size="500000"
# Loop through lines (genes) of gene summary file while skipping header
sed 1d $data_set_gene_summary_file | while read gene_id chrom_num tss gene_pheno_file covariate_file; do

	echo $gene_id
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	# Set OUT directory for this gene
	OUT=$fusion_input_dir$data_set_name"_"$gene_id"_1KG_only_fusion_input"

	# Run PLINK to set up gene for fusion
	plink --bfile $PRE_GENO$chrom_num --pheno $gene_pheno_file --make-bed --out $OUT --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1 --extract ${LDREF_dir}1000G.EUR.$chrom_num.snp_ids
	# Set FINAL_OUT directory for this gene
	FINAL_OUT=$fusion_output_dir$data_set_name"_"$gene_id"_1KG_only_fusion_output"
	TMP=$fusion_output_dir$data_set_name"_"$gene_id"_1KG_only_fusion_temp_output"
	Rscript "FUSION_EDITED.compute_weight_sparsity.R" --bfile $OUT --hsq_p 0.01 --tmp $TMP --covar $covariate_file --out $FINAL_OUT --verbose 1 --crossval 0 --save_hsq --PATH_gcta ${fusion_code_dir}"gcta_nr_robust" --models lasso_fixed_lambda
done









