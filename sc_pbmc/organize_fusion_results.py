import numpy as np 
import os
import sys
import pdb


def extract_gene_names_from_gene_summary_file(gene_summary_file):
	f = open(gene_summary_file)
	head_count = 0
	genes = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes.append(data[0])
	f.close()
	return np.asarray(genes)

def extract_fields_from_hsq_file(hsq_file):
	data = np.loadtxt(hsq_file,dtype=str)
	return float(data[1]), float(data[2]), float(data[3])

def convert_gene_fit_file_from_rdata_to_csv(gene_fit_file, new_gene_fit_file):
	os.system('Rscript convert_fusion_weights_from_rdata_to_tsv.R ' + gene_fit_file + ' ' + new_gene_fit_file)

def extract_fraction_non_zero_from_gene_fit(gene_fit_file, row_name):
	f = open(gene_fit_file)
	head_count = 0
	fraction_non_zero = -10
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# error checking
		if len(data) != 2:
			print('assumption eroror')
			pdb.set_trace()
		if data[0] == row_name:
			fraction_non_zero = float(data[1])
	f.close()
	if fraction_non_zero == -10:
		print('assumptino eroror')
		pdb.set_trace()
	return fraction_non_zero


def organize_fusion_results_in_single_data_set(data_set_name, gene_summary_file, fusion_output_dir, fusion_processed_output_dir):
	gene_names = extract_gene_names_from_gene_summary_file(gene_summary_file)

	# Keep track of which genes fusion has been run on
	fusion_run_on_this_gene_boolean = []

	# Make fusion data set summary file
	fusion_data_set_summary_file = fusion_processed_output_dir + data_set_name + '_1KG_only_fusion_results_summary.txt'
	t = open(fusion_data_set_summary_file,'w')
	# Last column is NA if not gene model was fit (h2_p > .01)
	t.write('gene_name\th2\th2_err\th2_p\tfraction_non_zero_weights\n')

	sig_h2 = []
	sig_fraction_non_zero = []

	# Loop through genes
	for gene_index, gene_name in enumerate(gene_names):
		# Stem of fusion results corresponding to this gene
		fusion_stem = fusion_output_dir + data_set_name + '_' + gene_name + '_1KG_only_fusion_output'
		# File containing hsq
		hsq_file = fusion_stem + '.hsq'
		# File containing gene fit model (only run if hsq is sig)
		gene_fit_file = fusion_stem + '_fraction_nonzero_weights.txt'

		# First check if fusion was even run on this gene
		if os.path.exists(hsq_file) == False:
			fusion_run_on_this_gene_boolean.append(False)
			continue
		fusion_run_on_this_gene_boolean.append(True)

		# Now get hsq
		h2, h2_err, h2_p = extract_fields_from_hsq_file(hsq_file)

		# Gene fit model was run
		if os.path.exists(gene_fit_file):
			'''
			# old
			# First convert RData weights file to tsv
			new_gene_fit_file = fusion_stem + 'wgt.tsv'
			convert_gene_fit_file_from_rdata_to_csv(gene_fit_file, new_gene_fit_file)
			# Load in gene weights
			gene_weights = (np.loadtxt(new_gene_fit_file,dtype=str)[:,1]).astype(float)
			# Fraction non-zero weights
			fraction_non_zero = sum(gene_weights != 0.0)/len(gene_weights)
			'''
			fraction_non_zero = extract_fraction_non_zero_from_gene_fit(gene_fit_file, '0.1')
			# Write to ouptut file
			t.write(gene_name + '\t' + str(h2) + '\t' + str(h2_err) + '\t' + str(h2_p) + '\t' + str(fraction_non_zero) + '\n')
			# Append to array
			sig_h2.append(h2)
			sig_fraction_non_zero.append(fraction_non_zero)
		else:
			t.write(gene_name + '\t' + str(h2) + '\t' + str(h2_err) + '\t' + str(h2_p) + '\tNA\n')

	t.close()
	fusion_run_on_this_gene_boolean = np.asarray(fusion_run_on_this_gene_boolean)
	sig_h2 = np.asarray(sig_h2)
	sig_fraction_non_zero = np.asarray(sig_fraction_non_zero)
	print(data_set_name)
	print(np.mean(sig_h2))
	print(np.mean(sig_fraction_non_zero))
	print('\n')









data_set_summary_file = sys.argv[1]
fusion_output_dir = sys.argv[2]
fusion_processed_output_dir = sys.argv[3]


# Loop through data sets
head_count = 0
f = open(data_set_summary_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	data_set_name = data[0]
	num_donors = int(data[5])
	num_genes = int(data[6])
	gene_summary_file = data[8]
	organize_fusion_results_in_single_data_set(data_set_name, gene_summary_file, fusion_output_dir, fusion_processed_output_dir)

f.close()