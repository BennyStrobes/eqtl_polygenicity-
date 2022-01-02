#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=5GB



gtex_tissue_sample_size_file="$1"
gtex_all_associations_dir="$2"
liftover_directory="$3"
gtex_hg19_all_associations_dir="$4"



tissue_name="Adipose_Subcutaneous"

if false; then
python liftover_gtex_association_in_specific_tissue.py $tissue_name $gtex_all_associations_dir $liftover_directory $gtex_hg19_all_associations_dir
fi



while read tissue_name sample_size; do
	if [ "$tissue_name" != "Adipose_Subcutaneous" ]; then
	if [ "$tissue_name" != "tissue_name" ]; then
		python liftover_gtex_association_in_specific_tissue.py $tissue_name $gtex_all_associations_dir $liftover_directory $gtex_hg19_all_associations_dir
	fi
	fi
done <$gtex_tissue_sample_size_file