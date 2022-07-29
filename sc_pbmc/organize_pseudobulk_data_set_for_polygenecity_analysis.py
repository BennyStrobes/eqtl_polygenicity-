import numpy as np 
import os
import sys
import pdb


def extract_h2_info(hsq_file):
	aa = np.loadtxt(hsq_file,dtype=str)
	if len(aa) != 4:
		print('assumption eroror')
		pdb.set_trace()
	return float(aa[1]), float(aa[3])

def filter_to_non_negative_heritable_genes(output_file, output_file2):
	f = open(output_file)
	t = open(output_file2,'w')

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		h2 = float(data[3])
		if h2 > 0:
			t.write(line + '\n')
	f.close()
	t.close()



polygenecity_analysis_input_dir = sys.argv[1] 
data_set_name = sys.argv[2]




gene_names = []
for file_name in os.listdir(polygenecity_analysis_input_dir):
	if file_name.startswith(data_set_name) and file_name.endswith('_genotype.npy'):
		ensamble_id = file_name.split(data_set_name + '_')[1].split('_geno')[0]
		gene_names.append(ensamble_id)

output_file = polygenecity_analysis_input_dir + data_set_name + '_organized_gene_info.txt'
t = open(output_file,'w')

t.write('gene_name\tgene_expression_file\tgenotype_file\tgene_heritability\tgene_heritability_p_value\n')

gene_names = np.sort(gene_names)

for gene_name in gene_names:
	genotype_file = polygenecity_analysis_input_dir + data_set_name + '_' + gene_name + '_genotype.npy'
	expr_file = polygenecity_analysis_input_dir + data_set_name + '_' + gene_name + '_residual_expression.npy'
	hsq_file = polygenecity_analysis_input_dir + data_set_name + '_' + gene_name + '_gene_heritability.hsq'

	if os.path.isfile(expr_file) and os.path.isfile(hsq_file) and os.path.isfile(genotype_file):
		h2, h2_p = extract_h2_info(hsq_file)
		t.write(gene_name + '\t' + expr_file + '\t' + genotype_file + '\t' + str(h2) + '\t' + str(h2_p) + '\n')
	else:
		print('assumption eroror')
		pdb.set_trace()

t.close()



output_file2 = polygenecity_analysis_input_dir + data_set_name + '_organized_gene_info_non_negative_h2.txt'
filter_to_non_negative_heritable_genes(output_file, output_file2)
