library(ashr)
library(mashr)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
set.seed(2020)

current = getwd()
setwd(current)
folder = "5_mashr"

CORR_TYPE <- "withCorr"
MODEL_FOLDER <- "eQTL_QN"
MODEL_TYPE <- c("eQTL_QN_combined")
COV_TYPE <- c("canonical_and_dataDriven")

## folders
INPUTS_dir <- paste0("inputs/",folder,"/",MODEL_FOLDER,"/")
OUTOUTS_dir <- paste0("outputs/",folder,"/",MODEL_FOLDER,"/",CORR_TYPE,"/")
system(paste0("mkdir -p ",OUTOUTS_dir))

for(i in 1:length(MODEL_TYPE)){
	MODEL <- MODEL_TYPE[i]

	## read in 
	beta_SE_matrix_ALL <- readRDS(paste0(INPUTS_dir,"eQTL_ALL_snps_tested.rds"))
	beta_SE_matrix_STRONG <- readRDS(paste0(INPUTS_dir,MODEL,".rds"))

	## get random subset of all tests performed
	random.subset = sample(1:nrow(beta_SE_matrix_ALL$Bhat), 200000)

	## estimate null correlation structure
	data.temp = mash_set_data(as.matrix(beta_SE_matrix_ALL$Bhat[random.subset,]), as.matrix(beta_SE_matrix_ALL$Shat[random.subset,]))
	Vhat = estimate_null_correlation_simple(data.temp)
	rm(data.temp)

	## transform to mashr format
	data_RANDOM = mash_set_data(as.matrix(beta_SE_matrix_ALL$Bhat[random.subset, ]), as.matrix(beta_SE_matrix_ALL$Shat[random.subset, ]), V = Vhat)
	data_STRONG = mash_set_data(as.matrix(beta_SE_matrix_STRONG$Bhat), as.matrix(beta_SE_matrix_STRONG$Shat), V = Vhat)

	for(j in 1:length(COV_TYPE)){
		COV <- COV_TYPE[j]
		OUTS_BY_COV_dir <- paste0("outputs/",folder,"/",MODEL_FOLDER,"/",CORR_TYPE,"/",COV,"/")
		system(paste0("mkdir -p ",OUTS_BY_COV_dir))

		## write covariance matrix 
		write.table(Vhat, paste0(OUTS_BY_COV_dir,MODEL,"_correlationMatrix_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		
		## set up the data-driven covariance matrix
		U.pca = cov_pca(data_STRONG, 5)
		## extreme deconvolution
		U.ed = cov_ed(data_STRONG, U.pca)
		saveRDS(U.pca, paste0(OUTS_BY_COV_dir,MODEL,"_covariancePCAs_",CORR_TYPE,"_",COV,"COV.rds"))
		saveRDS(U.ed, paste0(OUTS_BY_COV_dir,MODEL,"_U_ed_covarianceMatrices_",CORR_TYPE,"_",COV,"COV.rds"))

		## set up the canonical covariance matrix
		U.c = cov_canonical(data_RANDOM)
		saveRDS(U.c, paste0(OUTS_BY_COV_dir,MODEL,"_U_c_covarianceMatrices_",CORR_TYPE,"_",COV,"COV.rds"))

		## create infection covariance matrices
		U.NI = matrix(nrow = 10, ncol = 10)
		U.NI[,1:5] <- c(rep(1,5),rep(0,5))
		U.NI[,6:10] <- c(rep(0,10))
		U.flu = matrix(nrow = 10, ncol = 10)
		U.flu[,1:5] <- c(rep(0,10))
		U.flu[,6:10] <- c(rep(0,5),rep(1,5))
		U.infection <- list(U.NI, U.flu)
		names(U.infection) <- c("NI","flu")

		## run mash with both canonical and data-driven cov matrices -- fit the model 
		## fit mash to the random tests using both data-driven and canonical covariances -- we have to fit using a random set of tests, and not a dataset that is enriched for strong tests 
		## the outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time)
		m = mash(data_RANDOM, Ulist = c(U.ed, U.c, U.infection), outputlevel = 1)

		## use the fit from the previous run of mash by specifying g=get_fitted_g(m), fixg=TRUE to compute posterior summaries for any subset of tests
		m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)

		## assess sharing of significant signals among each pair of conditions by posterior means
		m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh = 0.1, factor = 0.5)
		write.table(m.pairwise_PM, paste0(OUTS_BY_COV_dir,MODEL,"_pairwise_sharing_PM_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)

		pdf(paste0(OUTS_BY_COV_dir,MODEL,"_pairwiseSharing_sigbyPosteriorMeans_",CORR_TYPE,"_",COV,"COV.pdf"), width=7, height=7)
		corrplot(m.pairwise_PM, method= 'color', cl.lim=c(0,1), type='full', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'pairwise sharing by magnitude \nof estimated effect sizes (PM)', mar = c(4,0,4,0), order = "hclust")
		dev.off()

		## extract the estimates of the mixture proportions for different types of covariance matrix
		pdf(paste0(OUTS_BY_COV_dir,MODEL,"_mixtureProportions_piEstimates_",CORR_TYPE,"_",COV,"COV.pdf"), width=8, height=8)
		par(mar = c(10,4,8,4))
		barplot(get_estimated_pi(m2), las = 2)
		dev.off()

		## write posterior outs
		write.table(get_lfsr(m2), paste0(OUTS_BY_COV_dir,MODEL,"_lfsr_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_lfdr(m2), paste0(OUTS_BY_COV_dir,MODEL,"_lfdr_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_pm(m2), paste0(OUTS_BY_COV_dir,MODEL,"_posteriorMeans_PM_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_psd(m2), paste0(OUTS_BY_COV_dir,MODEL,"_posteriorStandardDevs_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
	}
	print(MODEL)
}

