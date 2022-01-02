import numpy as np 
import os
import sys
import pdb


def snp_converter(input_file, output_file):
	f = open(input_file)
	t = open(output_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		snp_id = data[0] + ':' + data[3] + ':' + data[4] + ':' + data[5]
		snp_id2 = data[0] + ':' + data[3] + ':' + data[5] + ':' + data[4]
		t.write(snp_id + '\n')
		t.write(snp_id2 + '\n')
	f.close()
	t.close()


LDREF_dir = sys.argv[1]
one_KG_snps_dir = sys.argv[2]



for chrom_num in range(1,23):
	input_file = LDREF_dir + '1000G.EUR.' + str(chrom_num) + '.bim'
	output_file = one_KG_snps_dir + '1000G.EUR.' + str(chrom_num) + '.snp_ids'
	snp_converter(input_file, output_file)
