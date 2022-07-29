import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData
import h5py

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io







updated_individual_info_file = sys.argv[1]
input_h5py_file = sys.argv[2]
processed_sc_expression_dir = sys.argv[3]


# Load in processed-SC Ann-Data file
adata = sc.read_h5ad(input_h5py_file)

# For development purposes
#adata = adata[:10000,:]


sc.pp.neighbors(adata)


num_cells = adata.shape[0]

resolutions = [.001, .01, .1, 1, 3, 5, 10]

cluster_assignments = []
cluster_names = []
cluster_names.append('leiden_0')

for resolution in resolutions:
	cluster_names.append('leiden_' + str(resolution))
	sc.tl.leiden(adata, resolution=resolution)
	print(len(np.unique(adata.obs['leiden'])))
	cluster_assignments.append(adata.obs['leiden'])

cluster_names = np.asarray(cluster_names)
cluster_assignments =np.transpose(np.asarray(cluster_assignments))

cluster_assignments = np.hstack((np.zeros((num_cells,1)).astype(int), cluster_assignments))

cluster_assignments = np.vstack((cluster_names, cluster_assignments))



output_file = processed_sc_expression_dir + 'cluster_assignments_at_various_resolutions.txt'
np.savetxt(output_file, cluster_assignments, fmt="%s", delimiter='\t', comments='')






