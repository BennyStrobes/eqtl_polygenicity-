import numpy as np 
import os
import sys
import pdb



# Get mapping from genes to (chrom_num, position)
def get_mapping_from_gene_to_chromosome_position(gene_annotation_file):
	# Convert gene array to dictionary
	gene_mapping = {}
	#for gene in genes:
	#	gene_mapping[gene] = [0,0]
	# Fill in gene positions
	f = open(gene_annotation_file)
	used_genes = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Simple error check
		if len(data) != 8 and len(data) != 9:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		gene_id = data[0].split('.')[0]
		# Some more error checks
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		if data[6] != 'protein_coding' or data[7] != 'KNOWN':
			continue
		if data[1] == 'chrX' or data[1] == 'chrY' or data[1] == 'chrM':
			continue
		if gene_id not in gene_mapping:
			gene_mapping[gene_id] = [0,0]
		# Extract relevent info on gene: Chrom num and TSS
		chrom_num = int(data[1].split('hr')[1])
		strand = data[4]
		#gene_id = data[5]
		if strand == '+':
			tss = start
		elif strand == '-':
			tss = end
		else:
			print('assumption error while processing gene annotation file')
		# Add info to gene dictionary
		if gene_mapping[gene_id][0] != 0:
			if gene_mapping[gene_id][0] != chrom_num and gene_mapping[gene_id][1] != tss:
				print('odd event')
				pdb.set_trace()
				gene_mapping.pop(gene_id)
		else:
			gene_mapping[gene_id] = (chrom_num, tss)
	f.close()
	return gene_mapping	

def convert_from_sample_name_to_family_and_within_family_ids(sample_names):
	family_ids = []
	within_family_ids = []
	for sample_name in sample_names:
		temp = sample_name.split('_')[0]
		family_ids.append(temp)
		within_family_ids.append(temp)
	return np.asarray(family_ids), np.asarray(within_family_ids)

def make_gene_pheno_file(gene_pheno_file_name, family_ids, within_family_ids, gene_expression_vector):
	# Open file handle for gene pheno file
	t3 = open(gene_pheno_file_name,'w')

	# QUick erorr checking
	if len(family_ids) != len(gene_expression_vector):
		print('assumption eroror')
		pdb.set_trace()

	# Get number of smaples we have measurement for this gene
	num_samples = len(gene_expression_vector)

	# Loop through samples and print sample measurement to pheno file
	for sample_num in range(num_samples):
		t3.write(family_ids[sample_num] + '\t' + within_family_ids[sample_num] + '\t' + gene_expression_vector[sample_num] + '\n')
	# CLose filehandle
	t3.close()

def make_covariate_file(study_covariate_file, family_ids, within_family_ids, pb_covariate_file):
	# Open file handle for gene pheno file
	t3 = open(study_covariate_file,'w')

	# Load in covariates
	cov_data = np.loadtxt(pb_covariate_file, dtype=str, delimiter='\t')
	cov_data = cov_data[1:, 1:]

	num_samples = cov_data.shape[0]

	if num_samples != len(family_ids):
		print('assumption errror')
		pdb.set_trace()

	for sample_num in range(num_samples):
		t3.write(family_ids[sample_num] + '\t' + within_family_ids[sample_num] + '\t' + '\t'.join(cov_data[sample_num,:]) + '\n')

	# Close filehandle
	t3.close()

# Prepare gene expression for plink in single data set
def prepare_gene_expression_for_plink_in_single_data_set(data_set_name, pb_expression_file, pb_covariate_file, gene_to_tss, plink_input_dir, plink_gene_summary_file):
	# Open file handle for plink summary file
	t2 = open(plink_gene_summary_file,'w')
	# Header
	t2.write('Gene_id\tchrom_num\tTSS\tgene_pheno_file\tcovariate_file\n')

	# Load in gene expression data
	gene_expression_data_raw = np.loadtxt(pb_expression_file, dtype=str, delimiter='\t')
	sample_names = gene_expression_data_raw[0,1:]
	gene_names = gene_expression_data_raw[1:,0]
	gene_expression = gene_expression_data_raw[1:,1:]

	# Convert from sample names to family_id and within_family_id
	family_ids, within_family_ids = convert_from_sample_name_to_family_and_within_family_ids(sample_names)

	# Create Pheno file for this gene
	study_covariate_file = plink_input_dir + data_set_name + '_cov'
	make_covariate_file(study_covariate_file, family_ids, within_family_ids, pb_covariate_file)

	used_genes = {}
	# Loop through genes
	for gene_index, gene_name in enumerate(gene_names):
		# Quick error checking
		if gene_name in used_genes:
			print('gene repeat assumption erorr')
			pdb.set_trace()
		used_genes[gene_name] = 1
		
		# Chrom num and TSS of this gene
		chrom_num, gene_tss = gene_to_tss[gene_name]

		# Create Pheno file for this gene
		gene_pheno_file_name = plink_input_dir + data_set_name + '_' + gene_name + '_pheno'
		make_gene_pheno_file(gene_pheno_file_name, family_ids, within_family_ids, gene_expression[gene_index,:])

		# Update plink gene summary file with new gene
		t2.write(gene_name + '\t' + str(chrom_num) + '\t' + str(gene_tss) + '\t' + gene_pheno_file_name + '\t' + study_covariate_file + '\n')
	t2.close()


pseudobulk_expression_dir = sys.argv[1]
gene_annotation_file = sys.argv[2]
plink_input_dir = sys.argv[3]

# Threshold for deciding which studies pass
min_samples = 120

# Create mapping from gene to TSS
gene_to_tss = get_mapping_from_gene_to_chromosome_position(gene_annotation_file)

# File summarizing data sets
pseudobulk_data_set_summary_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary.txt'

# New file summarizing data sets
new_pseudobulk_data_set_summary_file = plink_input_dir + 'pseudobulk_data_set_summary_plink_ready.txt'

# Open input and output file handles
f = open(pseudobulk_data_set_summary_file)
t = open(new_pseudobulk_data_set_summary_file,'w')

# Loop through pseudobulk data sets
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Header
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\tplink_gene_summary_file\n')
		continue
	# Now we are at a line for the given data set
	# Extract relevent fields
	data_set_name = data[0]
	clustering_method = data[1]
	cluster_assignment_name = data[2]
	pb_expression_file = data[3]
	pb_covariate_file = data[4]
	num_samples = int(data[5])
	num_genes = int(data[6])

	# Filter out data sets with too few samples
	if num_samples < min_samples:
		continue

	# File_name of plink_gene_summary_file for this data set
	plink_gene_summary_file = plink_input_dir + data_set_name + '_gene_summary_file.txt'

	# Write to data set summary file
	t.write(line + '\t' + plink_gene_summary_file + '\n')

	# Prepare gene expression for plink in single data set
	prepare_gene_expression_for_plink_in_single_data_set(data_set_name, pb_expression_file, pb_covariate_file, gene_to_tss, plink_input_dir, plink_gene_summary_file)

# Close input and output file handles
f.close()
t.close()




