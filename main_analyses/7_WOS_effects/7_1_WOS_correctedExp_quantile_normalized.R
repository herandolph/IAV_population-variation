library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)

current = getwd()
folder = "7_WOS_effects"
setwd(current)
system(paste0("mkdir -p inputs/",folder,"/pseudobulk/"))
system(paste0("mkdir -p outputs/",folder,"/7_1_corrected_expression_outs/"))

## read in inputs
meta_data <- as.data.frame(read.csv(paste0("inputs/",folder,"/pseudobulk/meta_data_with_batches.txt"), header = TRUE, sep = ","))
rownames(meta_data) <- meta_data$Sample.ID
meta_data$age_Scale <- scale(meta_data$Age)
meta_data$BMI_Scale <- scale(meta_data$BMI)
meta_data$batch <- factor(meta_data$batch)
celltypes <- c("CD4_T","CD8_T","B","NK","monocytes_CD14_CD16")

for (i in 1:length(celltypes)){
	cell_type_i <- celltypes[i]
	meta_data_i <- meta_data
	reorder_names <- rownames(meta_data_i)

	## read in pseudobulk
	reads <- readRDS(paste0("inputs/",folder,"/pseudobulk/",cell_type_i,"_pseudobulk.rds"))

	if(length(reorder_names) == dim(reads)[2]){
		reads <- reads[reorder_names]
	}else{
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% colnames(reads), ]
		reads <- reads[, colnames(reads) %in% meta_data_i$Sample.ID]
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% colnames(reads), ]
		reorder_names <- rownames(meta_data_i)
		reads <- reads[reorder_names]
	}
	length(which(colnames(reads)!=rownames(meta_data_i)))

	## remove lowly-expressed genes and plot the effect
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design = model.matrix(~ 0 + batch, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]
	v <- voom(dge, design, plot = FALSE)

	## filter lowly-expressed genes
	if(cell_type_i == "CD4_T" || cell_type_i == "monocytes_CD14_CD16"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 1.5), ]
	}else if(cell_type_i == "CD8_T" || cell_type_i == "NK"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 2), ]
	}else{
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 3), ]
	}

	## voom after removal of lowly-expressed genes
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design <- model.matrix(~ 0 + batch, data = meta_data_i)
	v <- voom(dge, design, plot = FALSE)

	## fit
	fit <- lmFit(v, design)
	fit <- eBayes(fit)

	## get residuals (intercept == 0 right now) to regress out batch effect
	residuals <- residuals.MArrayLM(object = fit, v)
	length(which(colnames(residuals)!=rownames(meta_data_i)))

	## get average batch effect
	avg_batch_effect <- rowMeans(fit$coefficients)

	## add average batch effect back into residuals 
	corrected_expression <- apply(residuals,2,function(x){x + avg_batch_effect}) 
	weights <- v$weights
	colnames(weights) <- colnames(corrected_expression)
	rownames(weights) <- rownames(corrected_expression)

	## split corrected_expression matrix into T1 and T2 conditions
	list_T1 <- rownames(meta_data_i[meta_data_i$Blood.draw.time.point %in% "T1",])
	list_T2 <- rownames(meta_data_i[meta_data_i$Blood.draw.time.point %in% "T2",])
	corrected_expression_T1 <- corrected_expression[,colnames(corrected_expression) %in% list_T1]
	corrected_expression_T2 <- corrected_expression[,colnames(corrected_expression) %in% list_T2]

	quantile_norm <- function(df){

		## QUANTILE NORMALIZE WITHIN TIME POINT
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

	quantile_expression_T1 <- quantile_norm(corrected_expression_T1)
	quantile_expression_T2 <- quantile_norm(corrected_expression_T2)

	length(which(rownames(quantile_expression_T1)!=rownames(quantile_expression_T2)))
	quantile_expression <- cbind(quantile_expression_T1, quantile_expression_T2)
	order <- colnames(corrected_expression)
	quantile_expression <- as.data.frame(quantile_expression)
	quantile_expression <- quantile_expression[order]

	## write batch-corrected expression and weights
	write.table(quantile_expression, paste0("outputs/",folder,"/7_1_corrected_expression_outs/",cell_type_i,"_corrected_expression_QN.txt"), quote = FALSE, sep = ",")
	write.table(weights, paste0("outputs/",folder,"/7_1_corrected_expression_outs/",cell_type_i,"_weights.txt"), quote = FALSE, sep = ",")
	print(cell_type_i)
}

