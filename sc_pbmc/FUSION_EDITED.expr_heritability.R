# ==== TODO
# * Make sure BLUP/BSLMM weights are being scaled properly based on MAF

suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('methods'))

option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gemma", action="store", default="gemma", type='character',
              help="Path to plink executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--rn", action="store_true", default=FALSE,
              help="Rank-normalize the phenotype after all QC: [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"),			  
  make_option("--models", action="store", default="blup,lasso,top1,enet", type='character',
              help="Comma-separated list of prediction models [default: %default]\n
					Available models:\n
					top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
					blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
					bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
					lasso:\t LASSO regression (with heritability used as lambda)\n
					enet:\t Elastic-net regression (with mixing parameter of 0.5)\n")			  
		  
)

opt = parse_args(OptionParser(option_list=option_list))
models = unique( c(unlist(strsplit(opt$models,','))) )
M = length(models)

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}


# Perform i/o checks here:
files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$covar) ) files = c(files,opt$covar)

for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gcta could not be executed, set with --PATH_gcta\n" , sep='', file=stderr() )
	cleanup()
	q()
}


# ---

fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

# Make/fetch the phenotype file
if ( !is.na(opt$pheno) ) {
	pheno.file = opt$pheno
	pheno = read.table(pheno.file,as.is=T)
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	m = m[m.keep]
	pheno = pheno[m,]
} else {
	pheno.file = paste(opt$tmp,".pheno",sep='')
	pheno = fam[,c(1,2,6)]
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}



# Load in the covariates if needed
if ( !is.na(opt$covar) ) {
	covar = ( read.table(opt$covar,as.is=T,head=F) )
	#print(covar)
	#print(pheno)
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,2]) , paste(covar[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,]
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
	if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid)
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# --- HERITABILITY ANALYSIS
if ( is.na(opt$hsq_set) ) {
if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")

# 1. generate GRM
arg = paste( opt$PATH_plink," --allow-no-sex --bfile ",opt$tmp," --make-grm-bin --out ",opt$tmp,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 2. estimate heritability
if ( !is.na(opt$covar) ) {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",raw.pheno.file," --qcovar ",opt$covar," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
} else {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
}
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 3. evaluate LRT and V(G)/Vp
if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
	cat(opt$tmp,"does not exist, likely GCTA could not converge, skipping gene\n",file=stderr())
	q()
}

hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))

if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq,"LRT P-value:",hsq.pv,'\n')
if ( opt$save_hsq ) cat( opt$out , hsq , hsq.pv , '\n' , file=paste(opt$out,".hsq",sep='') )

} 

# read in genotypes
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2
sds = apply(genos$bed,2,sd)
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)

write.table(genos$bed, file=paste0(opt$tmp, "_geno_bed"), sep='\t', row.names=FALSE, col.names=FALSE)
