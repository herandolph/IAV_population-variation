library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(mgsub)

current = getwd()
folder = "4_cisRegression_model"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/binary_ancestry/"))

## inputs
expression_dir <- paste0("inputs/",folder,"/corrected_expression_QN/")
celltypes <- c("B","CD4_T","CD8_T","monocytes_combined","NK_combined")

## outputs
results_dir <- paste0("outputs/",folder,"/binary_ancestry/")

## read in meta data from individuals
meta_data <- read.table(paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt"), header = TRUE, sep = ",")
meta_data$capture <- as.factor(meta_data$capture)
meta_data$infection_status = factor(meta_data$infection_status, levels=c("NI","flu"))
meta_data$ethnicity = factor(meta_data$ethnicity, levels=c("EUR","AFR"))
meta_data <- meta_data[complete.cases(meta_data$YRI_GRCh38), ]

for (i in 1:length(celltypes)){
	cell_type_i <- celltypes[i]
	meta_data_i <- meta_data
	reorder_names <- rownames(meta_data_i)

	## read in corrected expression
	residuals <- read.table(paste0(expression_dir,cell_type_i,"_corrected_expression.txt"), header = TRUE, sep = ",", row.names = NULL)
	residuals[which(is.na(residuals)),][1] <- "NA"
	rownames(residuals) <- residuals$row.names; residuals$row.names <- NULL
	residuals <- residuals[colnames(residuals) %in% meta_data_i$infection_ID]

	## read in weights
	weights <- read.table(paste0(expression_dir,cell_type_i,"_weights.txt"), header = TRUE, sep = ",", row.names = NULL)
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

	## model popDE effects but with binary ancestry label (to compare with predicted population difference values downstream)
	design = model.matrix(~ age_Scale + counts + infection_status + ethnicity:infection_status, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, weights = weights, design)
	vfit <- eBayes(vfit)

	betas_NI = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:ethnicityAFR"))]); colnames(betas_NI)[1] <- "betas_NI"
	
	betas_flu = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:ethnicityAFR"))]); colnames(betas_flu)[1] <- "betas_flu"

	results <- cbind(betas_NI, betas_flu)
  	write.table(results, paste0(results_dir,cell_type_i,"_resultsALL.txt"), quote = FALSE, sep = ",")

	SEs_NI = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:ethnicityAFR"))]); colnames(SEs_NI) <- "SEs_NI"
  	SEs_flu = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:ethnicityAFR"))]); colnames(SEs_flu) <- "SEs_flu"
  	results_SE <- cbind(SEs_NI, SEs_flu)

  	write.table(results_SE, paste0(results_dir,cell_type_i,"_results_SEs.txt"), quote = FALSE, sep = ",")

  	print(cell_type_i)
}

