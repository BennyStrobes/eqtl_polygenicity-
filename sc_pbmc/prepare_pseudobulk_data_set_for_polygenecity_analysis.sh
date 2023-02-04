#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-60:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=25GB                         # Memory total in MiB (for all cores)



data_set_name="$1"
pseudobulk_expression_file="$2"
covariate_file="$3"
gene_info_file="$4"
plink_individual_file="$5"
polygenecity_analysis_input_dir_root="$6"
processed_genotype_dir="$7"


cis_window_size="100000"
gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"

module load R/3.5.1


polygenecity_analysis_input_dir=$polygenecity_analysis_input_dir_root$data_set_name"/"

mkdir $polygenecity_analysis_input_dir



sed 1d $gene_info_file | while read gene_id chrom_num tss strand gene_pheno_file; do

	echo $gene_id
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	# Set OUT directory for this gene
	OUT=$polygenecity_analysis_input_dir$data_set_name"_"$gene_id"_plink_input_tmp"
	# Run PLINK to set up gene for fusion
	plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --threads 1 --pheno $gene_pheno_file --make-bed --out $OUT --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1 

	FINAL_OUT=$polygenecity_analysis_input_dir$data_set_name"_"$gene_id"_gene_heritability"
	TMP=$polygenecity_analysis_input_dir$data_set_name"_"$gene_id"_tmp_output"
	Rscript "FUSION_EDITED.expr_heritability.R" --bfile $OUT --tmp $TMP --covar $covariate_file --out $FINAL_OUT --verbose 1 --save_hsq --PATH_gcta ${gcta_path} 

	gene_output_root=$polygenecity_analysis_input_dir$data_set_name"_"$gene_id"_"
	python3 prepare_pseudobulk_data_set_gene_for_polygenecity_analysis.py $TMP $gene_pheno_file $covariate_file $gene_output_root

	# Clean up
	rm $OUT*
	rm $TMP*

done



python3 organize_pseudobulk_data_set_for_polygenecity_analysis.py $polygenecity_analysis_input_dir $data_set_name




