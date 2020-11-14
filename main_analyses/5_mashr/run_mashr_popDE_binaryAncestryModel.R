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
MODEL_FOLDER <- "popDE_binaryAncestryModel_QN"
MODEL_TYPE <- c("popDE_binaryAncestryModel_QN_combined")
COV_TYPE <- c("canonical_and_dataDriven")

## folders
INPUTS_dir <- paste0("inputs/",folder,"/",MODEL_FOLDER,"/")
OUTOUTS_dir <- paste0("outputs/",folder,"/",MODEL_FOLDER,"/",CORR_TYPE,"/")
system(paste0("mkdir -p ",OUTOUTS_dir))

for(i in 1:length(MODEL_TYPE)){
	MODEL <- MODEL_TYPE[i]

	## read in
	beta_SE_matrix <- readRDS(paste0(INPUTS_dir,MODEL,".rds"))

	## transform to mashr format
	data = mash_set_data(as.matrix(beta_SE_matrix$Bhat), as.matrix(beta_SE_matrix$Shat))

	## estimate null correlations
	V = estimate_null_correlation_simple(data)
	data_V = mash_update_data(data, V = V)

	for(j in 1:length(COV_TYPE)){
		COV <- COV_TYPE[j]
		OUTS_BY_COV_dir <- paste0("outputs/",folder,"/",MODEL_FOLDER,"/",CORR_TYPE,"/",COV,"/")
		system(paste0("mkdir -p ",OUTS_BY_COV_dir))

		## write covariance matrix 
		write.table(V, paste0(OUTS_BY_COV_dir,MODEL,"_correlationMatrix_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)

		## set up the canonical covariance matrix
		U.c = cov_canonical(data_V)
		saveRDS(U.c, paste0(OUTS_BY_COV_dir,MODEL,"_U_c_covarianceMatrices_",CORR_TYPE,"_",COV,"COV.rds"))
		
		## set up the data-driven covariance matrix
		m.1by1 = mash_1by1(data_V)
		strong = get_significant_results(m.1by1, 0.05)
		betas <- as.matrix(beta_SE_matrix$Bhat)
		strong_betas <- betas[strong,]
			
		U.pca = cov_pca(data_V, 5, subset=strong)
		## extreme deconvolution
		U.ed = cov_ed(data_V, U.pca, subset=strong)
		saveRDS(U.pca, paste0(OUTS_BY_COV_dir,MODEL,"_covariancePCAs_",CORR_TYPE,"_",COV,"COV.rds"))
		saveRDS(U.ed, paste0(OUTS_BY_COV_dir,MODEL,"_U_ed_covarianceMatrices_",CORR_TYPE,"_",COV,"COV.rds"))

		## run mash with both canonical and data-driven cov matrices -- fit the model 
		## sample from mash object (sample from posteriors)
		m = mash(data_V, c(U.c, U.ed), algorithm.version = 'R', posterior_samples = 100)

		## assess sharing of significant signals among each pair of conditions by posterior means
		m.pairwise_PM <- get_pairwise_sharing(m, lfsr_thresh = 0.1, factor = 0.5)
		write.table(m.pairwise_PM, paste0(OUTS_BY_COV_dir,MODEL,"_pairwise_sharing_PM_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)

		pdf(paste0(OUTS_BY_COV_dir,MODEL,"_pairwiseSharing_sigbyPosteriorMeans_",CORR_TYPE,"_",COV,"COV.pdf"), width=7, height=7)
		corrplot(m.pairwise_PM, method= 'color', cl.lim=c(0,1), type='full', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'pairwise sharing by magnitude \nof estimated effect sizes (PM)', mar = c(4,0,4,0), order = "hclust")
		dev.off()

		## extract the estimates of the mixture proportions for different types of covariance matrix
		pdf(paste0(OUTS_BY_COV_dir,MODEL,"_mixtureProportions_piEstimates_",CORR_TYPE,"_",COV,"COV.pdf"), width=8, height=8)
		par(mar = c(10,4,8,4))
		barplot(get_estimated_pi(m), las = 2)
		dev.off()

		## write posterior outs
		write.table(get_lfsr(m), paste0(OUTS_BY_COV_dir,MODEL,"_lfsr_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_lfdr(m), paste0(OUTS_BY_COV_dir,MODEL,"_lfdr_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_pm(m), paste0(OUTS_BY_COV_dir,MODEL,"_posteriorMeans_PM_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
		write.table(get_psd(m), paste0(OUTS_BY_COV_dir,MODEL,"_posteriorStandardDevs_",CORR_TYPE,"_",COV,"COV.txt"), quote = FALSE)
	}
	print(MODEL)
}

