#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=5GB




tissue_name="$1"
tissue_sumstats_file="$2"
tissue_sample_size="$3"
ldsc_source_code_dir="$4"
sldsc_input_data_dir="$5"
gtex_tissue_sldsc_processed_input_dir="$6"
gtex_tissue_sldsc_results_dir="$7"


if false; then
python prepare_sldsc_data.py $tissue_name $tissue_sumstats_file $tissue_sample_size $sldsc_input_data_dir $gtex_tissue_sldsc_processed_input_dir
fi


module load python/2.7-anaconda
study_file=$gtex_tissue_sldsc_processed_input_dir"Adipose_Subcutaneous_ENSG00000000457.13_sumstats"
sldsc_processed_data_dir="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/sldsc_processed_data/"
python ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${gtex_tissue_sldsc_results_dir}"testing."

if false; then
python ldsc.py --h2 ${study_file} --ref-ld ${sldsc_processed_data_dir}"baselineLD.1" --w-ld ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.1" --overlap-annot --print-coefficients --frqfile ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC.1" --out ${gtex_tissue_sldsc_results_dir}"testing.2"
fi