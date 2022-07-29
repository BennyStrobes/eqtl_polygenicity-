#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --nodes=1



processed_genotype_dir="$1"
genotype_data_dir="$2"
filtered_sample_info_file="$3"


plink_format_filter_sample_info_file=$processed_genotype_dir"plink_formatted_filtered_sample_info.txt"
python3 extract_ordered_individuals_in_plink_format.py $filtered_sample_info_file $plink_format_filter_sample_info_file


plink2 --pfile $genotype_data_dir"inds_v2_header.maf10" --threads 1 --keep $plink_format_filter_sample_info_file --indiv-sort f $plink_format_filter_sample_info_file --maf .1 --make-pgen --out $processed_genotype_dir"inds_v2_sample_filter"


plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --threads 1 --freq --out $processed_genotype_dir"inds_v2_sample_filter_frq"


plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --indep-pairwise 200 50 0.25 --out $processed_genotype_dir"inds_v2_sample_filter_independent_snps"


plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --extract $processed_genotype_dir"inds_v2_sample_filter_independent_snps.prune.in" --pca 3 --out $processed_genotype_dir"inds_v2_sample_filter_ind_snp_pcs"



















################
# OLD
################



#########################
# Filter vcf file to individuals we have sc-rna-seq for
########################
donor_filtered_merged_vcf=$processed_genotype_dir"clues_immvar_donor_filter_merged.vcf.gz"
if false; then
vcftools --gzvcf $genotype_input_vcf_file --keep $sc_rna_seq_individual_file --recode --stdout | gzip -c > $donor_filtered_merged_vcf
fi

#########################
# Filter vcf file to sites with maf .1 and no missing
########################
donor_site_filtered_merged_vcf=$processed_genotype_dir"clues_immvar_donor_site_filter_merged.vcf.gz"
if false; then
vcftools --gzvcf $donor_filtered_merged_vcf --remove-filtered-all --max-missing 1 --maf .1 --recode --stdout | gzip -c > $donor_site_filtered_merged_vcf
fi


#########################
# Extract list of sample names that pass genotype filtering
########################
sample_names=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_samples_that_pass_filters.txt"
if false; then
bcftools query -l $donor_site_filtered_merged_vcf > $sample_names
fi

#########################
# Convert to Plink
#########################
donor_site_filtered_merged_plink_stem=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_plink"
# Convert from to plink format
if false; then
plink --vcf $donor_site_filtered_merged_vcf --out $donor_site_filtered_merged_plink_stem
fi

# Make seperate plink file for each chromosome
if false; then
for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	donor_site_filtered_merged_chrom_specific_plink_stem=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_plink_chr_"$chrom_num
	plink --vcf $donor_site_filtered_merged_vcf --chr $chrom_num  --out $donor_site_filtered_merged_chrom_specific_plink_stem
done
fi




#########################
# Genotype PCs
#########################


# Filter to independent sites
donor_site_filtered_merged_independent_plink_stem=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_independent_plink"
if false; then
plink --bfile $donor_site_filtered_merged_plink_stem --indep-pairwise 200 50 0.25 --out $donor_site_filtered_merged_independent_plink_stem
fi

# Compute Genotype PCs
if false; then
plink --bfile $donor_site_filtered_merged_plink_stem --extract $donor_site_filtered_merged_independent_plink_stem".prune.in" --pca 7 --out $donor_site_filtered_merged_independent_plink_stem
fi




