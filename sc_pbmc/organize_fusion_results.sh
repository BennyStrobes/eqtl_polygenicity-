#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared



data_set_summary_file="$1"
fusion_output_dir="$2"
fusion_processed_output_dir="$3"



python organize_fusion_results.py $data_set_summary_file $fusion_output_dir $fusion_processed_output_dir