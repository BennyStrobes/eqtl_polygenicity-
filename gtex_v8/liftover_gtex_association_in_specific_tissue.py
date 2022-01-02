import numpy as np 
import os
import sys
import pdb


def association_file_to_bed_file(association_file, bed_file):
	f = open(association_file)
	t = open(bed_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_id = data[1]
		snp_info = snp_id.split('_')
		chrom_num = snp_info[0]
		pos = snp_info[1]
		pos2 = int(pos) + 1
		t.write(chrom_num + '\t' + pos + '\t' + str(pos2) + '\n')
	f.close()
	t.close()

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory):
    stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg38ToHg19.over.chain.gz ' + output_file + ' ' + missing_file
    os.system(stringer)

def recreate_all_association_file_in_hg19_format(hg38_association_file, hg19_association_file, liftover_output_file, liftover_missing_file):
	missing = {}
	f = open(liftover_missing_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		missing[data[0] + '_' + data[1]] = 1
	f.close()
	f = open(hg38_association_file)
	t = open(hg19_association_file,'w')
	g = open(liftover_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		snp_id = data[1]
		snp_info = snp_id.split('_')
		snp_key = snp_info[0] + '_' + snp_info[1]
		if snp_key in missing:
			continue
		hg19_info = next(g).rstrip().split('\t')
		new_snp_id = hg19_info[0] + '_' + hg19_info[1] + '_' + snp_info[2] + '_' + snp_info[3]
		t.write(data[0] + '\t' + new_snp_id + '\t' + '\t'.join(data[2:]) + '\n')
	f.close()
	g.close()
	t.close()


def get_num_lines(file_name):
	counter = 0
	f = open(file_name)
	for line in f:
		counter = counter + 1
	f.close()
	return counter

def error_checking(liftover_output_file, hg19_association_file):
	num_lines_lift = get_num_lines(liftover_output_file)
	num_lines_association = get_num_lines(hg19_association_file) - 1
	if num_lines_association != num_lines_lift:
		print('assumption erororo')


# Command line args
tissue_name=sys.argv[1]
gtex_all_associations_dir= sys.argv[2]
liftover_directory= sys.argv[3]
gtex_hg19_all_associations_dir=sys.argv[4]



# All association file in hg38 format
hg38_association_file = gtex_all_associations_dir + tissue_name + '.allpairs.txt'

# First convert hg38 all association to temporary bed format
temp_hg38_bed = gtex_hg19_all_associations_dir + tissue_name + '_temp_bed_hg38.txt'
association_file_to_bed_file(hg38_association_file, temp_hg38_bed)


# Run liftover
liftover_output_file = gtex_hg19_all_associations_dir + tissue_name + '_liftover_output.txt'
liftover_missing_file = gtex_hg19_all_associations_dir + tissue_name + '_liftover_missing.txt'
run_liftover(temp_hg38_bed, liftover_output_file, liftover_missing_file, liftover_directory)

# All association file in hg19 format
hg19_association_file = gtex_hg19_all_associations_dir + tissue_name + '.hg19.allpairs.txt'
recreate_all_association_file_in_hg19_format(hg38_association_file, hg19_association_file, liftover_output_file, liftover_missing_file)


# Quick error checking
error_checking(liftover_output_file, hg19_association_file)

# Remove temporary files
os.system('rm ' + liftover_output_file)
os.system('rm ' + liftover_missing_file)
os.system('rm ' + temp_hg38_bed)





