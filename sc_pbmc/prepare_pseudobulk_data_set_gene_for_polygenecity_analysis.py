import numpy as np 
import os
import sys
import pdb





tmp_output_file = sys.argv[1]
gene_pheno_file = sys.argv[2]
covariate_file = sys.argv[3]
gene_output_root = sys.argv[4]


#print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
#print('AM I USING THE CORRECT VERSION OF PYTHON???')
#print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')


genotype_file = tmp_output_file + '_geno_bed'
geno = np.loadtxt(genotype_file,delimiter='\t')


residual_phenotype = tmp_output_file + '.pheno.resid'
raw_resid_phenotype = np.loadtxt(residual_phenotype, dtype=str)
resid_pheno = raw_resid_phenotype[:,2].astype(float)

if geno.shape[0] != len(resid_pheno):
	print('assumption eroror')
	pdb.set_trace()

np.save(gene_output_root + 'genotype.npy',geno)
np.save(gene_output_root + 'residual_expression.npy',resid_pheno)

##########
'''
pheno = np.loadtxt(gene_pheno_file,dtype=str)
pheno = pheno[:,1].astype(float)

cov = np.loadtxt(covariate_file,dtype=str)
cov = cov[:,2:].astype(float)

from sklearn.linear_model import LinearRegression
model = LinearRegression()
aa= model.fit(cov, pheno)
resid2 =pheno - aa.predict(cov)
resid2 = (resid2-np.mean(resid2))/np.std(resid2)
pdb.set_trace()
'''