import numpy as np 
import os
import sys
import pdb
import gzip
from scipy.stats import chi2





def get_valid_variants(valid_variant_file):
	f = gzip.open(valid_variant_file)
	head_count = 0
	valid_variants = {}
	for line in f:
		line = line.decode('UTF-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_id = 'chr' + data[0] + '_' + data[2]
		rs_id = data[1]
		if snp_id in valid_variants:
			print('assumption eoroooror')
			pdb.set_trace()
		valid_variants[snp_id] = rs_id
	f.close()
	return valid_variants

def p_to_z(P):
    '''Convert P-value and N to standardized beta.'''
    return np.sqrt(chi2.isf(P, 1))


def create_study_summary_statistic_files(tissue_sumstats_file, valid_variants, tissue_name, tissue_sample_size, chrom_num, t, gtex_tissue_sldsc_processed_input_dir):
	chrom_string = 'chr' + str(chrom_num)
	f = open(tissue_sumstats_file)
	head_count = 0
	gene_dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_id = data[1]
		if snp_id.startswith(chrom_string + '_') == False:
			continue
		snp_info = snp_id.split('_')
		short_snp_id = snp_info[0] + '_' + snp_info[1]
		if short_snp_id not in valid_variants:
			continue
		gene_id = data[0]
		snp_id = data[1]
		if gene_id not in gene_dicti:
			gene_dicti[gene_id] = {}
			gene_dicti[gene_id]['snps'] = {}
			gene_dicti[gene_id]['info'] = []
		if snp_id in gene_dicti[gene_id]['snps']:
			print('assumption eroror')
			pdb.set_trace()
		gene_dicti[gene_id]['snps'][snp_id] = 1
		gene_dicti[gene_id]['info'].append((snp_id, float(data[6])))
	f.close()

	sorted_genes = np.sort(np.asarray(list(gene_dicti.keys())))
	for gene in sorted_genes:
		num_snps = len(gene_dicti[gene]['snps'])
		t.write(gene + '\t' + chrom_string + '\t' + tissue_sample_size + '\t' + str(num_snps) + '\n')
		gene_sumstats_file = gtex_tissue_sldsc_processed_input_dir + tissue_name + '_' + gene + '_sumstats'
		g = open(gene_sumstats_file,'w')
		g.write('SNP\tA1\tA2\tZ\tN\n')
		for ele in gene_dicti[gene]['info']:
			snp_id = ele[0]
			pvalue = ele[1]
			snp_info = snp_id.split('_')
			short_snp_id = snp_info[0] + '_' + snp_info[1]
			rs_id = valid_variants[short_snp_id]
			a1 = snp_info[2]
			a2 = snp_info[3]
			z_score = p_to_z(pvalue)
			g.write(rs_id + '\t' + a1 + '\t' + a2 + '\t' + str(z_score) + '\t' + tissue_sample_size + '\n')
		g.close()

	return t

tissue_name = sys.argv[1]
tissue_sumstats_file = sys.argv[2]
tissue_sample_size = sys.argv[3]
sldsc_input_data_dir = sys.argv[4]
gtex_tissue_sldsc_processed_input_dir = sys.argv[5]


# Create file to keep track of all studies (genes)
study_file = gtex_tissue_sldsc_processed_input_dir + tissue_name + '_cis_studies.txt'
t = open(study_file,'w')
t.write('study_name\tchrom_num\tsample_size\tnum_snps\n')


for chrom_num in range(1,2):
	# First create mapping from snp_id to rs_id for valid variants
	valid_variant_file = sldsc_input_data_dir + '1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.' + str(chrom_num) + '.l2.ldscore.gz'
	valid_variants = get_valid_variants(valid_variant_file)
	# Create study file for each gene
	t = create_study_summary_statistic_files(tissue_sumstats_file, valid_variants, tissue_name, tissue_sample_size, chrom_num, t, gtex_tissue_sldsc_processed_input_dir)


t.close()