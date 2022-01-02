import numpy as np 
import os
import sys
import pdb
import h5py
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io









#######################
# Command line args
#######################
processed_sc_expression_dir = sys.argv[1]
processed_genotype_dir = sys.argv[2]



output_file = processed_genotype_dir + 'sc_rna_seq_ea_individual_list.txt'

# Load in processed-SC Ann-Data file
input_h5py_file = processed_sc_expression_dir + 'scran_normalization_hvg_2000_regress_batch_True_2.h5ad'
adata = sc.read_h5ad(input_h5py_file)


indiz = np.asarray(adata.obs['ind_cov'])
popz = np.asarray(adata.obs['pop_cov'])

unique_individuals_eur = {}
unique_individuals_aa = {}
unique_individuals_asian = {}
unique_individuals_hispanic = {}


for i, indi in enumerate(indiz):
	if popz[i] == 'European':
		unique_individuals_eur[indi] = 1
	elif popz[i] == 'African American':
		unique_individuals_aa[indi] = 1
	elif popz[i] == 'Asian':
		unique_individuals_asian[indi] = 1
	elif popz[i] == 'Hispanic':
		unique_individuals_hispanic[indi] = 1
	else:
		print('assumption eororor')
		pdb.set_trace()

t = open(output_file,'w')
eur_donors = np.sort(np.asarray([*unique_individuals_eur]))

for eur_donor in eur_donors:
	t.write(eur_donor + '\n')
t.close()