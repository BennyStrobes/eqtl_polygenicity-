import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import pandas as pd
import scipy.special

def generate_variance_grid(grid_start, grid_end, multiplicative_grid_scaling_factor):
	# Create fixed variance grid
	variance_grid_arr = []
	grid_value =grid_start
	grid_converged = False
	while grid_converged == False:
		variance_grid_arr.append(grid_value)
		new_grid_value = grid_value*multiplicative_grid_scaling_factor
		if new_grid_value > grid_end:
			grid_converged = True 
		else:
			grid_value = new_grid_value
	variance_grid_arr = np.asarray(variance_grid_arr)
	return variance_grid_arr

def load_in_gene_info(gene_info_file):
	return pd.read_csv(gene_info_file,sep='\t')

def initialize_data(gene_df, variance_grid, gene_heritability_normalization):
	# Initialize output variables
	gene_to_beta_mu = {}
	gene_to_beta_var = {}
	gene_to_Z = {}
	gene_to_residual_precision = {}
	gene_to_h2 = {}

	M = len(variance_grid)

	num_genes = gene_df.shape[0]
	# Now loop through num_genes
	for gene_iter in range(num_genes):
		#print(gene_iter)
		# Extract relevent info for this window
		gene_name = gene_df['gene_name'][gene_iter]
		gene_expr_file = gene_df['gene_expression_file'][gene_iter]
		genotype_file = gene_df['genotype_file'][gene_iter]
		# Load in expr data to get num snps
		G = np.load(genotype_file)
		Y = np.load(gene_expr_file)
		num_snps = G.shape[1]

		gene_to_beta_mu[gene_name] = np.zeros(num_snps)
		gene_to_beta_var[gene_name] = np.ones(num_snps)
		gene_to_residual_precision[gene_name] = 1.0

		Z = np.ones((num_snps, M))/M
		gene_to_Z[gene_name] = Z

		if gene_heritability_normalization == 'True_estimate':
			local_residual_variance, local_per_snp_h2, local_gene_h2 = local_heritability_inference(Y, G, max_iters=2000, convergence_thresh=1e-6)
			gene_to_h2[gene_name] = local_per_snp_h2
	
	delta_alpha = np.ones(M)

	return gene_to_beta_mu, gene_to_beta_var, gene_to_Z, gene_to_residual_precision, delta_alpha, gene_to_h2

def update_beta_and_Z(Y, X, beta_mu, Z, variance_grid, expected_ln_pi, residual_precision, per_snp_h2):
	resid = Y - np.dot(X, beta_mu)
	P = X.shape[1]
	beta_squared = np.zeros(P)

	for pp in range(P):
		resid = resid + (X[:, pp]*beta_mu[pp])

		a_terms = -.5*residual_precision*np.sum(np.square(X[:,pp])) - (.5/variance_grid)
		b_term = np.sum(residual_precision*resid*X[:,pp])

		mixture_alpha_var = -1.0/(2.0*a_terms)
		mixture_alpha_mu = b_term*mixture_alpha_var


		un_normalized_lv_weights = expected_ln_pi - (.5*np.log(variance_grid)) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
		#if np.sum(np.isnan(un_normalized_lv_weights)) > 0:
			#pdb.set_trace()
		un_normalized_lv_weights = un_normalized_lv_weights - scipy.special.logsumexp(un_normalized_lv_weights)
		Z[pp,:] = np.exp(un_normalized_lv_weights)


		beta_mu[pp] = np.sum(mixture_alpha_mu*Z[pp,:])
		#beta_var[pp] = np.sum(mixture_alpha_var*Z[pp,:]) + np.sum(np.square(mixture_alpha_mu)*Z[pp,:]) - np.sum(np.square(mixture_alpha_mu*Z[pp,:]))
		beta_squared[pp] = np.sum(Z[pp,:]*(np.square(mixture_alpha_mu) + mixture_alpha_var))

		resid = resid - (X[:, pp]*beta_mu[pp])

	return beta_mu, Z, beta_squared


def update_beta_and_Z_shell(gene_df, gene_to_beta_mu, gene_to_Z, gene_to_residual_precision, gene_to_h2, variance_grid, expected_ln_pi, gene_heritability_normalization):
	num_genes = gene_df.shape[0]

	for gene_iter in range(num_genes):		
		# Extract relevent info for this window
		gene_name = gene_df['gene_name'][gene_iter]
		gene_expr_file = gene_df['gene_expression_file'][gene_iter]
		genotype_file = gene_df['genotype_file'][gene_iter]
		gene_h2 = gene_df['gene_heritability'][gene_iter]

		# Load in data
		Y = np.load(gene_expr_file)
		#num_snps = len(Y)
		#Y = Y.reshape((num_snps,1))
		G = np.load(genotype_file)

		num_snps = G.shape[1]

		# Calculate expected per snp heritability
		per_snp_h2 = gene_h2/num_snps
		if gene_heritability_normalization == 'True':
			beta_mu, Z, beta_squared = update_beta_and_Z(Y, G, gene_to_beta_mu[gene_name], gene_to_Z[gene_name], variance_grid*per_snp_h2, expected_ln_pi, gene_to_residual_precision[gene_name], per_snp_h2)
		elif gene_heritability_normalization == 'False':
			beta_mu, Z, beta_squared = update_beta_and_Z(Y, G, gene_to_beta_mu[gene_name], gene_to_Z[gene_name], variance_grid, expected_ln_pi, gene_to_residual_precision[gene_name], per_snp_h2)
		elif gene_heritability_normalization == 'True_estimate':
			per_snp_h2 = gene_to_h2[gene_name]
			beta_mu, Z, beta_squared = update_beta_and_Z(Y, G, gene_to_beta_mu[gene_name], gene_to_Z[gene_name], variance_grid*per_snp_h2, expected_ln_pi, gene_to_residual_precision[gene_name], per_snp_h2)
		else:
			print('gene heritability method currently not implemented')
			pdb.set_trace()

		residual_precision = update_residual_precision_v3(Y, G, beta_mu, beta_squared)


		gene_to_beta_mu[gene_name] = beta_mu
		gene_to_Z[gene_name] = Z
		gene_to_residual_precision[gene_name] = residual_precision



	return gene_to_beta_mu, gene_to_Z, gene_to_residual_precision

def update_delta_alpha(gene_df, gene_to_Z, M):
	delta_alpha = np.zeros(M) + 1e-30
	num_genes = gene_df.shape[0]

	for gene_iter in range(num_genes):		
		# Extract relevent info for this window
		gene_name = gene_df['gene_name'][gene_iter]

		gene_Z = gene_to_Z[gene_name]

		delta_alpha = delta_alpha + np.sum(gene_Z,axis=0)
	return delta_alpha

def compute_excess_kurtosis_from_mixture_of_gaussians(row_tissue_pi, variance_grid):
	numerator = np.sum(row_tissue_pi*np.square(variance_grid))
	denominator = np.square(np.sum(row_tissue_pi*variance_grid))
	kurt = 3.0*numerator/denominator - 3.0
	return kurt

def print_residual_variance_output_file(gene_df, gene_to_residual_precision, residual_variance_output_file):
	
	t = open(residual_variance_output_file,'w')
	t.write('gene_name\tresidual_variance\n')

	num_genes = gene_df.shape[0]

	for gene_iter in range(num_genes):		
		# Extract relevent info for this window
		gene_name = gene_df['gene_name'][gene_iter]	

		residual_variance = 1.0/gene_to_residual_precision[gene_name]

		t.write(gene_name + '\t' + str(residual_variance) + '\n')

	t.close()

def inference(gene_df, variance_grid, gene_to_beta_mu, gene_to_Z, gene_to_residual_precision, gene_to_h2, delta_alpha, gene_heritability_normalization, output_root, max_iters=1000):

	for itera in range(max_iters):
		expected_ln_pi = scipy.special.psi(delta_alpha) - scipy.special.psi(np.sum(delta_alpha))
		gene_to_beta_mu, gene_to_Z, gene_to_residual_precision = update_beta_and_Z_shell(gene_df, gene_to_beta_mu, gene_to_Z, gene_to_residual_precision, gene_to_h2, variance_grid, expected_ln_pi, gene_heritability_normalization)

		delta_alpha = update_delta_alpha(gene_df, gene_to_Z, len(delta_alpha))

		delta_expected = delta_alpha/np.sum(delta_alpha)
		print('ITERATION ' + str(itera))
		print(delta_expected)
		expected_var = np.sum(variance_grid*delta_expected)
		excess_kurtosis = compute_excess_kurtosis_from_mixture_of_gaussians(delta_expected, variance_grid)

		if np.mod(itera, 5) == 0:
			delta_expected_output_file = output_root + 'temp_' + str(itera) + '_delta_expected.txt'
			np.savetxt(delta_expected_output_file, np.vstack((variance_grid.astype(str), delta_expected.astype(str))), fmt="%s", delimiter='\t')
			kurtosis_output_file = output_root + 'temp_' + str(itera) + '_expected_kurtosis.txt'
			np.savetxt(kurtosis_output_file, np.asarray([expected_var, excess_kurtosis]).astype(str), fmt="%s", delimiter='\t')
			residual_variance_output_file = output_root + 'temp_' + str(itera) + '_residual_variance.txt'
			print_residual_variance_output_file(gene_df, gene_to_residual_precision, residual_variance_output_file)


def univariate_beta_update(X, y, tau, beta_mu, beta_var, residual_precision):
	num_snps = X.shape[1]
	resid_y = y - np.dot(X, beta_mu)
	for kk in np.random.permutation(range(num_snps)):
		resid_y = resid_y + (X[:,kk]*beta_mu[kk])
		beta_var[kk] = 1.0/(tau + residual_precision*np.sum(np.square(X[:,kk])))
		beta_mu[kk] = (beta_var[kk])*np.sum(residual_precision*X[:,kk]*resid_y)
		resid_y = resid_y - (X[:,kk]*beta_mu[kk])
	return beta_mu, beta_var

def local_heritability_inference(Y, G, max_iters=2000, convergence_thresh=1e-6):
	num_snps = G.shape[1]
	beta_mu = np.zeros(num_snps)
	beta_var = np.ones(num_snps)*1e-10
	num_snps = G.shape[1]
	init_gene_h2 = .01
	expected_tau = 1.0/(init_gene_h2/num_snps)
	residual_precision = 1.0 - init_gene_h2
	curr_gene_h2 = num_snps/expected_tau
	converged = False
	for itera in range(max_iters):
		old_gene_h2 = curr_gene_h2

		beta_mu, beta_var = univariate_beta_update(G, Y, expected_tau, beta_mu, beta_var, residual_precision)
		tau_a = num_snps/2.0 
		tau_b = 0.5*np.sum(np.square(beta_mu) + beta_var)
		expected_tau = tau_a/tau_b
		#residual_precision = update_residual_precision(Y, G, beta_mu, beta_var, hyper_param=0.0)
		residual_precision = update_residual_precision_v2(Y, G, beta_mu, beta_var, hyper_param=0.0)

		curr_gene_h2 = num_snps/expected_tau
		delta = np.abs(curr_gene_h2 - old_gene_h2)
		if delta < convergence_thresh:
			converged = True
			break	

	if converged == False:
		print('Convergence issue: Model did not converge after ' + str(max_iters) + ' iterations')

	local_residual_variance = 1.0/residual_precision
	local_per_snp_h2 = 1.0/expected_tau
	local_gene_h2 = (1.0/expected_tau)*num_snps



	return local_residual_variance, local_per_snp_h2, local_gene_h2

def multivariate_beta_update(G, Y, expected_tau, residual_precision):
	S = np.linalg.inv(np.diag(np.ones(G.shape[1])*expected_tau) + residual_precision*np.dot(np.transpose(G), G))
	mu = residual_precision*np.dot(np.dot(S, np.transpose(G)), Y)
	return mu, S



def update_residual_precision(Y, G, beta_mu, beta_var, hyper_param=1e-16):

	beta_squared = np.square(beta_mu) + beta_var
	beta_pred = np.dot(G, beta_mu)
	G_squared = np.square(G)

	resid = np.square(Y) + np.dot(G_squared, beta_squared)
	resid = resid - (2.0*Y*beta_pred)
	# Now add terms with interactions between factors
	resid = resid + (beta_pred*beta_pred - np.dot(G_squared, np.square(beta_mu)))

	new_alpha = hyper_param + (len(Y)/2.0)
	new_beta = hyper_param + (np.sum(resid)/2.0)

	precision = new_alpha/new_beta
	return precision

def update_residual_precision_v2(Y, G, beta_mu, beta_var, hyper_param=1e-16):

	#resid = np.sum(np.square(Y)) - np.dot(np.dot(beta_mu, np.linalg.inv(beta_S)), beta_mu)
	resid = np.sum(np.square(Y)) 
	resid = resid - np.sum(2.0*Y*np.dot(G, beta_mu))
	resid = resid + np.sum((np.dot(np.transpose(G), G))*(np.dot(np.reshape(beta_mu, (len(beta_mu),1)), np.reshape(beta_mu, (1,len(beta_mu)))) + np.diag(beta_var)))

	new_alpha = hyper_param + (len(Y)/2.0)
	new_beta = hyper_param + ((resid)/2.0)

	precision = new_alpha/new_beta
	return precision

def update_residual_precision_v3(Y, G, beta_mu, beta_squared, hyper_param=1e-16):

	#beta_squared = np.square(beta_mu) + beta_var
	beta_pred = np.dot(G, beta_mu)
	G_squared = np.square(G)

	resid = np.square(Y) + np.dot(G_squared, beta_squared)
	resid = resid - (2.0*Y*beta_pred)
	# Now add terms with interactions between factors
	resid = resid + (beta_pred*beta_pred - np.dot(G_squared, np.square(beta_mu)))

	new_alpha = hyper_param + (len(Y)/2.0)
	new_beta = hyper_param + (np.sum(resid)/2.0)

	precision = new_alpha/new_beta

	if precision < 0:
		pdb.set_trace()

	return precision



def update_residual_precision_multivariate(Y, G, beta_mu, beta_S, hyper_param=0.0):
	#resid = (np.sum(np.square(Y - np.dot(G, beta_mu)))/2.0) + (np.trace(np.dot(np.dot(np.transpose(G), G), beta_S))/2.0)
	#resid = np.sum(np.square(Y)) - np.dot(np.dot(beta_mu, np.linalg.inv(beta_S)), beta_mu)
	resid = np.sum(np.square(Y)) 
	resid = resid - np.sum(2.0*Y*np.dot(G, beta_mu))
	resid = resid + np.sum((np.dot(np.transpose(G), G))*(np.dot(np.reshape(beta_mu, (len(beta_mu),1)), np.reshape(beta_mu, (1,len(beta_mu)))) + beta_S))

	new_alpha = hyper_param + (len(Y)/2.0)
	new_beta = hyper_param + ((resid)/2.0)

	precision = new_alpha/new_beta
	return precision

def local_mixture_model_heritability_inference(Y, G,max_iters=10000, convergence_thresh=1e-8):
	num_snps = G.shape[1]
	init_gene_h2 = .01
	expected_tau = 1.0/(init_gene_h2/num_snps)
	residual_precision = 1.0 - init_gene_h2
	residual_precision = 1.0
	# Create variance grid
	grid_start=1e-5
	grid_end=1.0
	multiplicative_grid_scaling_factor=5
	variance_grid = generate_variance_grid(grid_start, grid_end, multiplicative_grid_scaling_factor)

	N = G.shape[0]
	P = G.shape[1]
	M = len(variance_grid)
	beta_mu = np.zeros(P)
	beta_var = np.ones(P)
	Z = np.ones((P, M))/M
	delta_alpha = np.ones(M)
	converged = False
	for itera in range(max_iters):
		old_residual_precision = residual_precision
		expected_ln_pi = scipy.special.psi(delta_alpha) - scipy.special.psi(np.sum(delta_alpha))
		beta_mu, beta_var, Z, beta_squared = update_beta_and_Z(Y, G, beta_mu, beta_var, Z, variance_grid, expected_ln_pi, residual_precision)
		delta_alpha = np.sum(Z,axis=0) 
		residual_precision = update_residual_precision_v3(Y, G, beta_mu, beta_squared, hyper_param=0.0)
		diff = np.abs(residual_precision - old_residual_precision)
		print(delta_alpha/np.sum(delta_alpha))
		if diff < convergence_thresh:
			converged = True
			break
	if converged == False:
		print('Convergence issue: Model did not converge after ' + str(max_iters) + ' iterations')
	pi = delta_alpha/np.sum(delta_alpha)
	expected_var = pi*variance_grid
	residual_var = 1.0/residual_precision
	pdb.set_trace()
	return pi, expected_var*P, residual_var

def local_multivariate_heritability_inference(Y, G, max_iters=10000, convergence_thresh=1e-8):
	num_snps = G.shape[1]
	residual_precision = 1.0 - init_gene_h2



	converged = False
	for itera in range(max_iters):
		old_gene_h2 = curr_gene_h2
		beta_mu, beta_S = multivariate_beta_update(G, Y, expected_tau, residual_precision)
		tau_a = num_snps/2.0 
		tau_b = 0.5*(np.sum(beta_mu*beta_mu)+np.trace(beta_S))
		expected_tau = tau_a/tau_b
		residual_precision = update_residual_precision_multivariate(Y, G, beta_mu, beta_S, hyper_param=0.0)

		curr_gene_h2 = num_snps/expected_tau

		delta = np.abs(curr_gene_h2 - old_gene_h2)
		if delta < convergence_thresh:
			converged = True
			break

	if converged == False:
		print('Convergence issue: Model did not converge after ' + str(max_iters) + ' iterations')

	local_residual_variance = 1.0/residual_precision
	local_per_snp_h2 = 1.0/expected_tau
	local_gene_h2 = (1.0/expected_tau)*num_snps

	return local_residual_variance, local_per_snp_h2, local_gene_h2



		



def local_gene_heritability_estimation(gene_df):
	num_genes = gene_df.shape[0]

	for gene_iter in range(num_genes):
		print(gene_iter)	
		# Extract relevent info for this window
		gene_name = gene_df['gene_name'][gene_iter]
		gene_expr_file = gene_df['gene_expression_file'][gene_iter]
		genotype_file = gene_df['genotype_file'][gene_iter]
		gene_h2 = gene_df['gene_heritability'][gene_iter]

		# Load in data
		Y = np.load(gene_expr_file)
		G = np.load(genotype_file)
		G_stand = np.copy(G)

		Y_stand = (Y - np.mean(Y))/np.std(Y)
		for kk in range(G.shape[1]):
			G_stand[:, kk] = (G[:,kk] - np.mean(G[:,kk]))/np.std(G[:,kk])



		#local_residual_var, local_per_snp_h2, local_h2 = local_heritability_inference(Y, G,max_iters=300)
		#local_mv_residual_var, local_mv_per_snp_h2, local_mv_h2 = local_multivariate_heritability_inference(Y_stand, G_stand,max_iters=10000, convergence_thresh=1e-8)
		pi, expected_var, residual_var = local_mixture_model_heritability_inference(Y_stand, G_stand,max_iters=2000, convergence_thresh=1e-6)
		print(pi)
		print(str(expected_var) + '\t' + str(1.0 - residual_var))
		#print(str(1.0-local_residual_var) + '\t' + str(local_h2) + '\t' + str(gene_h2))
		#print(str(1.0-local_mv_residual_var) + '\t' + str(local_mv_h2) + '\t' + str(gene_h2))


data_set_name = sys.argv[1]
gene_info_file = sys.argv[2]
output_root = sys.argv[3]
gene_heritability_normalization = sys.argv[4]

gene_df = load_in_gene_info(gene_info_file)

#local_gene_heritability_estimation(gene_df)

if gene_heritability_normalization == "True":
	grid_start=1e-5
	grid_end=1e5
	multiplicative_grid_scaling_factor=3.0
elif gene_heritability_normalization == "False":
	grid_start=1e-10
	grid_end=1e0
	multiplicative_grid_scaling_factor=3.0	
elif gene_heritability_normalization == 'True_estimate':
	grid_start=1e-5
	grid_end=1e5
	multiplicative_grid_scaling_factor=3.0	
else:
	print('gene heritability method currently not implemented')
	pdb.set_trace()
variance_grid = generate_variance_grid(grid_start, grid_end, multiplicative_grid_scaling_factor)

gene_to_beta_mu, gene_to_beta_var, gene_to_Z, gene_to_residual_precision, delta_alpha, gene_to_h2 = initialize_data(gene_df, variance_grid, gene_heritability_normalization)

inference(gene_df, variance_grid, gene_to_beta_mu, gene_to_Z, gene_to_residual_precision, gene_to_h2, delta_alpha, gene_heritability_normalization, output_root)

