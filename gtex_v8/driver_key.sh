#############################################################
# Exploratory analysis of eQTL polygenecity in GTEx v8 data
#############################################################




#######################
# Input data
#######################

# Directory containing all eQTL associations (one file for each tissue)
gtex_all_associations_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/"

# File containing gtex tissues and their sample sizes
gtex_tissue_sample_size_file="/work-zfs/abattle4/bstrober/eqtl_polygenicity/input_data/gtex_tissue_sample_size.txt"


# Directory containing ldsc source code
ldsc_source_code_dir="/work-zfs/abattle4/bstrober/tools/ldsc/"

# Directory containing sldsc input data
sldsc_input_data_dir="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/input_data/sldsc_input_data/"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/work-zfs/abattle4/bstrober/tools/liftOver_x86/"


#######################
# Output data
#######################
# Root directory for all data generated in this analysis
output_root="/work-zfs/abattle4/bstrober/eqtl_polygenicity/gtex_v8/"

# Directory containing all gtex v8 associations lifted over to hg18
gtex_hg19_all_associations_dir=$output_root"gtex_v8_all_associations_hg19/"

# Directory containing processed sldsc input data
gtex_tissue_sldsc_processed_input_dir=$output_root"gtex_tissue_sldsc_processed_input/"

# Directory containing sldsc results when applied to each tissues eqtls
gtex_tissue_sldsc_results_dir=$output_root"gtex_tissue_sldsc_results/"





########################################
# Liftover gtex associations from hg38 to hg19
########################################
if false; then
sbatch liftover_gtex_associations.sh $gtex_tissue_sample_size_file $gtex_all_associations_dir $liftover_directory $gtex_hg19_all_associations_dir
fi


########################################
# Run LDSC for each tissue's eqtls
########################################
tissue_name="Adipose_Subcutaneous"
tissue_sample_size="581"
tissue_sumstats_file=$gtex_hg19_all_associations_dir$tissue_name".hg19.allpairs.txt"
sh run_sldsc_in_a_tissue.sh $tissue_name $tissue_sumstats_file $tissue_sample_size $ldsc_source_code_dir $sldsc_input_data_dir $gtex_tissue_sldsc_processed_input_dir $gtex_tissue_sldsc_results_dir


