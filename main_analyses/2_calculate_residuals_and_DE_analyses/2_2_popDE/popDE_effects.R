library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)

current = getwd()
folder = "2_calculate_residuals_and_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/2_2_popDE_results/"))

## inputs
meta_data_file <- paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt")
residuals_dir <- paste0("inputs/",folder,"/corrected_expression_QN/")
celltypes <- c("B","CD4_T","CD8_T","monocytes_combined","NK_combined")

## outputs
results_dir <- paste0("outputs/",folder,"/2_2_popDE_results/")

## read in meta data from individuals and remove individuals with no quantitative admixture value
meta_data <- read.table(paste0(meta_data_file), header = TRUE, sep = ",")
meta_data$capture <- as.factor(meta_data$capture)
meta_data$indiv_ID <- as.factor(meta_data$indiv_ID)
meta_data$infection_status = factor(meta_data$infection_status, levels=c("NI","flu"))
meta_data <- meta_data[complete.cases(meta_data$YRI_GRCh38), ]

for (i in 1:length(celltypes)){
	cell_type_i <- celltypes[i]
	meta_data_i <- meta_data
	reorder_names <- rownames(meta_data_i)

	## read in corrected expression
	residuals <- read.table(paste0(residuals_dir,cell_type_i,"_corrected_expression.txt"), header = TRUE, sep = ",", row.names = NULL)
	residuals[which(is.na(residuals)),][1] <- "NA"
	rownames(residuals) <- residuals$row.names; residuals$row.names <- NULL
	residuals <- residuals[colnames(residuals) %in% meta_data_i$infection_ID]

	## read in weights
	weights <- read.table(paste0(residuals_dir,cell_type_i,"_weights.txt"), header = TRUE, sep = ",", row.names = NULL)
	weights[which(is.na(weights)),][1] <- "NA"
	rownames(weights) <- weights$row.names; weights$row.names <- NULL
	weights <- weights[colnames(weights) %in% meta_data_i$infection_ID]

	if(length(reorder_names) == dim(residuals)[2]){
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}else{
		correct_names <- colnames(residuals)
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% correct_names,]
		weights <- weights[,colnames(weights) %in% correct_names]
		reorder_names <- rownames(meta_data_i)
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}

	## subset correct meta_data bimodal proportion
	counts <- subset(meta_data_i, select = c(paste0(cell_type_i,"_geneProp")))
	colnames(counts)[1] <- "counts"
	meta_data_i <- meta_data_i[,c(1:22)]
	meta_data_i <- cbind(meta_data_i, counts)

	length(which(colnames(residuals)!=rownames(meta_data_i)))
	length(which(colnames(weights)!=rownames(meta_data_i)))
	length(which(colnames(residuals)!=colnames(weights)))

	## model popDE effects using a nested model
	design = model.matrix(~ age_Scale + counts + infection_status + YRI_Scale:infection_status, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, weights = weights, design)
	vfit <- eBayes(vfit)

	## collect outputs
	infection_betas_sign = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusflu"))]); colnames(infection_betas_sign)[1] <- "infection_betas_sign"

	betas_NI = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:YRI_Scale"))]); colnames(betas_NI)[1] <- "betas_NI"
	p_values_NI = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:YRI_Scale"))]); colnames(p_values_NI)[1] <- "pvalues_NI"
	fdrs_NI = as.data.frame(p.adjust(p_values_NI[,1], method = "BH")); colnames(fdrs_NI)[1] <- "fdrs_NI"

	betas_flu = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:YRI_Scale"))]); colnames(betas_flu)[1] <- "betas_flu"
	p_values_flu = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:YRI_Scale"))]); colnames(p_values_flu)[1] <- "pvalues_flu"
	fdrs_flu = as.data.frame(p.adjust(p_values_flu[,1], method = "BH")); colnames(fdrs_flu)[1] <- "fdrs_flu"

	results <- cbind(infection_betas_sign, betas_NI, p_values_NI, fdrs_NI, betas_flu, p_values_flu, fdrs_flu)
  	assign(paste0(cell_type_i,"_results_popDE"), results)
  	write.table(results, paste0(results_dir,cell_type_i,"_resultsALL.txt"), quote = FALSE, sep = ",")

  	## calculate standard errors
  	SEs_NI = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:YRI_Scale"))]); colnames(SEs_NI) <- "SEs_NI"

  	SEs_flu = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:YRI_Scale"))]); colnames(SEs_flu) <- "SEs_flu"

  	results_SE <- cbind(SEs_NI, SEs_flu)
  	write.table(results_SE, paste0(results_dir,cell_type_i,"_results_SEs.txt"), quote = FALSE, sep = ",")

  	print(cell_type_i)
}
