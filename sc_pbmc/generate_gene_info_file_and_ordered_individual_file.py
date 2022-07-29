import numpy as np 
import os
import sys
import pdb



def get_gene_names_from_gene_expression_file(gene_expression_data_file):
	raw = np.loadtxt(gene_expression_data_file, dtype=str, delimiter='\t')
	gene_names = raw[1:,0]
	# Quick error checking
	if len(gene_names) != len(np.unique(gene_names)):
		print('assumption error')
		pdb.set_trace()
	gene_name_to_expression = {}
	for ii, gene_name in enumerate(gene_names):
		gene_name_to_expression[gene_name] = raw[(ii+1),1:]
	individual_names = raw[0,1:]
	return gene_names, gene_name_to_expression, individual_names

def get_gene_name_to_to_chrom_num_tss_and_strand(gene_names, gene_annotation_file):
	dicti = {}
	for gene_name in gene_names:
		dicti[gene_name] = ('chrom', 'tss', 'strand')

	f = open(gene_annotation_file)
	head_count = 0
	used = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if line.startswith('#'):
			continue
		if len(data) != 9:
			print('assumption error')
			pdb.set_trace()
		sequence_class = data[2]
		if sequence_class != 'gene':
			continue
		if data[0] == 'chrX' or data[0] == 'chrY' or data[0] == 'chrM':
			continue

		gene_info = data[8].split(';')

		ensamble_id, gene_id, gene_type, gene_status = extract_gene_name_and_type_from_info(gene_info)

		short_ensamble_id = ensamble_id.split('.')[0]

		if short_ensamble_id not in dicti:
			continue

		if short_ensamble_id in used:
			print('assumption error')
			pdb.set_trace()
		used[short_ensamble_id] = 1

		chrom_num = data[0]
		strand = data[6]

		if strand == '+':
			tss = data[3]
		elif strand == '-':
			tss = data[4]
		else:
			print('assumption error on strands')
			pdb.set_trace()

		dicti[short_ensamble_id] = (chrom_num, tss, strand)
	f.close()

	# More error checking
	for gene_name in gene_names:
		if dicti[gene_name][0] == 'chrom':
			print('assumption eroror')
			pdb.set_trace()


	return dicti





def extract_gene_name_and_type_from_info(gene_info):
	ensamble_id = 'NA'
	gene_id = 'NA'
	gene_type = 'NA'
	gene_status = 'NA'
	for stringer in gene_info:
		if stringer.startswith('ID='):
			ensamble_id = stringer.split('ID=')[1]
		if stringer.startswith('gene_type='):
			gene_type = stringer.split('gene_type=')[1]
		if stringer.startswith('gene_status='):
			gene_status = stringer.split('gene_status=')[1]
		if stringer.startswith('gene_name='):
			gene_id = stringer.split('gene_name=')[1]
	if ensamble_id == 'NA' or gene_id == 'NA' or gene_type == 'NA' or gene_status == 'NA':
		print('assumption erroror')
		pdb.set_trace()
	return ensamble_id, gene_id, gene_type, gene_status

def make_gene_pheno_file(gene_pheno_file, expr_vec, individual_names):
	if len(expr_vec) != len(individual_names):
		print('assumption erooror')
		pdb.set_trace()

	t = open(gene_pheno_file,'w')
	t.write('#IID' + '\tEXPR\n')
	for ii, individual in enumerate(individual_names):
		t.write(individual + '\t' + expr_vec[ii] + '\n')

	t.close()


def generate_gene_info_file(gene_names, gene_name_to_info, gene_name_to_expression, gene_info_file, pseudobulk_expression_dir, study_name, individual_names):
	t2 = open(gene_info_file,'w')
	t2.write('gene_name\tchrom_num\ttss\tstrand\tgene_pheno_file\n')

	for gene_name in gene_names:
		info = gene_name_to_info[gene_name]

		gene_pheno_file = pseudobulk_expression_dir + study_name + '_' + gene_name + '_pheno'
		make_gene_pheno_file(gene_pheno_file, gene_name_to_expression[gene_name], individual_names)


		t2.write(gene_name + '\t' + info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + gene_pheno_file + '\n')
	t2.close()

def generate_plink_individual_file(covariate_file, plink_individual_file):
	raw = np.loadtxt(covariate_file, dtype=str, delimiter='\t')
	individuals = raw[1:, 0]

	t2 = open(plink_individual_file,'w')
	for individual in individuals:
		t2.write(individual + '\n')
	t2.close()


def generate_plink_covariate_file(covariate_file, plink_covariate_file):
	head_count = 0
	f = open(covariate_file)
	t = open(plink_covariate_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		t.write('0' + '\t' + line + '\n')
	f.close()
	t.close()


pseudobulk_expression_dir = sys.argv[1]
gene_annotation_file = sys.argv[2]






# Pseudobulk data set summary file
pseudobulk_data_set_summary_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary.txt'
pseudobulk_data_set_summary_final_file = pseudobulk_expression_dir + 'pseudobulk_data_set_summary_final.txt'

f = open(pseudobulk_data_set_summary_file)
t = open(pseudobulk_data_set_summary_final_file, 'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'gene_info_file\tplink_individual_file\tplink_covariate_file\n')
		continue
	study_name = data[0]
	gene_expression_data_file = data[3]
	covariate_file = data[4]
	num_expression_samples = int(data[5])
	print(study_name)
	if num_expression_samples < 105:
		continue

	gene_names, gene_name_to_expression, individual_names = get_gene_names_from_gene_expression_file(gene_expression_data_file)

	plink_covariate_file = pseudobulk_expression_dir + study_name + '_plink_covariate_file.txt'
	generate_plink_covariate_file(covariate_file, plink_covariate_file)

	
	gene_name_to_info = get_gene_name_to_to_chrom_num_tss_and_strand(gene_names, gene_annotation_file)

	# Generate gene info file
	gene_info_file = pseudobulk_expression_dir + study_name + '_gene_info.txt'
	generate_gene_info_file(gene_names, gene_name_to_info, gene_name_to_expression, gene_info_file, pseudobulk_expression_dir, study_name, individual_names)

	# Generate plink individual file
	plink_individual_file = pseudobulk_expression_dir + study_name + '_plink_individual_names.txt'
	generate_plink_individual_file(covariate_file, plink_individual_file)

	# Print new file names to output file
	t.write(line + '\t' + gene_info_file + '\t' + plink_individual_file + '\t' + plink_covariate_file + '\n')

f.close()
t.close()








# Pseudobulk data set summary file
pseudobulk_data_set_summary_file = pseudobulk_expression_dir + 'cell_type_matched_pseudobulk_data_set_summary.txt'
pseudobulk_data_set_summary_final_file = pseudobulk_expression_dir + 'cell_type_matched_pseudobulk_data_set_summary_final.txt'

f = open(pseudobulk_data_set_summary_file)
t = open(pseudobulk_data_set_summary_final_file, 'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'gene_info_file\tplink_individual_file\tplink_covariate_file\n')
		continue
	study_name = data[0]
	gene_expression_data_file = data[3]
	covariate_file = data[4]
	num_expression_samples = int(data[5])
	print(study_name)
	if num_expression_samples < 105:
		continue

	gene_names, gene_name_to_expression, individual_names = get_gene_names_from_gene_expression_file(gene_expression_data_file)

	plink_covariate_file = pseudobulk_expression_dir + study_name + '_plink_covariate_file.txt'
	generate_plink_covariate_file(covariate_file, plink_covariate_file)

	
	gene_name_to_info = get_gene_name_to_to_chrom_num_tss_and_strand(gene_names, gene_annotation_file)

	# Generate gene info file
	gene_info_file = pseudobulk_expression_dir + study_name + '_gene_info.txt'
	generate_gene_info_file(gene_names, gene_name_to_info, gene_name_to_expression, gene_info_file, pseudobulk_expression_dir, study_name, individual_names)

	# Generate plink individual file
	plink_individual_file = pseudobulk_expression_dir + study_name + '_plink_individual_names.txt'
	generate_plink_individual_file(covariate_file, plink_individual_file)

	# Print new file names to output file
	t.write(line + '\t' + gene_info_file + '\t' + plink_individual_file + '\t' + plink_covariate_file + '\n')

f.close()
t.close()