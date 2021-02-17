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
system(paste0("mkdir -p inputs/",folder,"/corrected_expression_QN/"))

## inputs
meta_data_file <- paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt")
pseudobulk_dir <- paste0("inputs/",folder,"/2_1_correctedExp_and_infection_DE/pseudobulk/")
celltypes <- c("B","CD4_T","CD8_T","monocytes_combined","NK_combined")

## outputs
residuals_dir <- paste0("inputs/",folder,"/corrected_expression_QN/")
results_dir <- paste0("outputs/",folder,"/results/")

## read in meta data from individuals
meta_data <- read.table(paste0(meta_data_file), header = TRUE, sep = ",")
meta_data$capture <- as.factor(meta_data$capture)
meta_data$indiv_ID <- as.factor(meta_data$indiv_ID)
meta_data$infection_status = factor(meta_data$infection_status, levels=c("NI","flu"))

for (i in 1:length(celltypes)){
	cell_type_i <- celltypes[i]
	meta_data_i <- meta_data
	reorder_names <- rownames(meta_data_i)

	## read in pseudobulk
	reads <- readRDS(paste0(pseudobulk_dir,cell_type_i,"_pseudobulk.rds"))

	if(length(reorder_names) == dim(reads)[2]){
		reads <- reads[reorder_names]
	}else{
		correct_names <- colnames(reads)
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% correct_names,]
		reorder_names <- rownames(meta_data_i)
		reads <- reads[reorder_names]
	}

	## remove lowly-expressed genes
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design = model.matrix(~ 0 + capture, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]

	v <- voom(dge, design, plot = FALSE)

	## median logCPM filtering
	if(cell_type_i == "CD4_T" || cell_type_i == "monocytes" || cell_type_i == "monocytes_combined"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 1.5), ]
	}else if(cell_type_i == "B" || cell_type_i == "CD8_T"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 2.5), ]
	}else{
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 4.0), ]
		}

	## voom after removal of lowly expressed genes
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design <- model.matrix(~ 0 + capture, data = meta_data_i)
	v <- voom(dge, design, plot = FALSE)

	## fit
	fit <- lmFit(v, design)
	fit <- eBayes(fit)

	## get residuals (intercept == 0 right now) to regress out capture effect
	residuals <- residuals.MArrayLM(object = fit, v)

	length(which(colnames(residuals)!=rownames(meta_data_i)))

	## get average capture effect
	avg_capture_effect <- rowMeans(fit$coefficients)

	## add average capture effect back into residuals 
	corrected_expression <- apply(residuals,2,function(x){x + avg_capture_effect}) 
	weights <- v$weights
	colnames(weights) <- colnames(corrected_expression)
	rownames(weights) <- rownames(corrected_expression)

	## split corrected_expression matrix into NI and flu conditions
	list_NI <- rownames(meta_data_i[meta_data_i$infection_status %in% "NI",])
	list_flu <- rownames(meta_data_i[meta_data_i$infection_status %in% "flu",])
	corrected_expression_NI <- corrected_expression[,colnames(corrected_expression) %in% list_NI]
	corrected_expression_flu <- corrected_expression[,colnames(corrected_expression) %in% list_flu]

	quantile_norm <- function(df){

		## QUANTILE NORMALIZE WITHIN CONDITION
		quantile_expression <- matrix(, nrow = nrow(df), ncol = ncol(df))

		for (j in 1:nrow(quantile_expression)){
			exp <- df[j,]
			exp_QN <- qqnorm(exp, plot = FALSE)
			exp_QN <- exp_QN$x
			quantile_expression[j,] <- exp_QN
		}

		rownames(quantile_expression) <- rownames(df)
		colnames(quantile_expression) <- colnames(df)
		return(quantile_expression)
	}

	quantile_expression_NI <- quantile_norm(corrected_expression_NI)
	quantile_expression_flu <- quantile_norm(corrected_expression_flu)
	
	length(which(rownames(quantile_expression_NI)!=rownames(quantile_expression_flu)))
	quantile_expression <- cbind(quantile_expression_NI, quantile_expression_flu)
	
	## write capture-corrected expression and weights to model popDE and popDR later
	write.table(quantile_expression, paste0(residuals_dir,cell_type_i,"_corrected_expression.txt"), quote = FALSE, sep = ",")

	## write weights
	write.table(weights, paste0(residuals_dir,cell_type_i,"_weights.txt"), quote = FALSE, sep = ",")
  	print(cell_type_i)
}


