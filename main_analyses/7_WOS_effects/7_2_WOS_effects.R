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
system(paste0("mkdir -p inputs/",folder,"/corrected_expression_QN/"))
system(paste0("mkdir -p outputs/",folder,"/7_2_WOS_effects/"))

## source function
setwd("common_functions")
source("permFDR.R")

setwd(current)

## read in inputs
meta_data <- as.data.frame(read.csv(paste0("inputs/",folder,"/corrected_expression_QN/meta_data_with_propBimodal.txt"), header = TRUE, sep = ","))
rownames(meta_data) <- meta_data$Sample.ID
meta_data$age_Scale <- scale(meta_data$Age)
meta_data$BMI_Scale <- scale(meta_data$BMI)
meta_data$batch <- factor(meta_data$batch)
meta_data$Sex <- factor(meta_data$Sex, levels = c("Female","Male"))

## change NA values in race/ethnicity to most common values
meta_data$Ethnicity <- mapvalues(meta_data$Ethnicity, from = "Unknown / Not Reported", to = "NOT Hispanic or Latino")
meta_data$Race <- mapvalues(meta_data$Race, from = "Unknown / Not Reported", to = "White")
meta_data$Race <- factor(meta_data$Race, levels = c("White","Asian","Black or African American","Native Hawaiian or Other Pacific Islander","American Indian/Alaska Native","More Than One Race"))
meta_data$Ethnicity <- factor(meta_data$Ethnicity, levels = c("NOT Hispanic or Latino","Hispanic or Latino"))
celltypes <- c("CD4_T","CD8_T","B","NK","monocytes_CD14_CD16")

## remove NA -- one individual is missing WOS score
meta_data <- meta_data[complete.cases(meta_data$Who.Ordinal.Scale), ]
meta_data <- meta_data[!(meta_data$Study.Subject.ID %in% c("INCOV132","INCOV145")),]

for (i in 1:length(celltypes)){

	cell_type_i <- celltypes[i]
	meta_data_i <- meta_data

	if(cell_type_i == "B"){
		meta_data_i <- meta_data_i[!(meta_data_i$Study.Subject.ID %in% c("INCOV086","INCOV109","INCOV121","INCOV144")),]
	}

	reorder_names <- rownames(meta_data_i)

	## read in permutations
	permutations_T1 <- read.table(paste0("inputs/",folder,"/corrected_expression_QN/permutations/permutations_",cell_type_i,"_T1.txt"), header = TRUE, sep = ","); rownames(permutations_T1) <- permutations_T1$genes
	permutations_T2 <- read.table(paste0("inputs/",folder,"/corrected_expression_QN/permutations/permutations_",cell_type_i,"_T2.txt"), header = TRUE, sep = ","); rownames(permutations_T2) <- permutations_T2$genes

	## read in quantile normalized corrected expression
	residuals <- read.table(paste0("inputs/",folder,"/corrected_expression_QN/",cell_type_i,"_corrected_expression_QN.txt"), header = TRUE, sep = ",", check.names = FALSE)
	residuals <- residuals[colnames(residuals) %in% meta_data_i$Sample.ID]

	## read in weights
	weights <- read.table(paste0("inputs/",folder,"/corrected_expression_QN/",cell_type_i,"_weights.txt"), header = TRUE, sep = ",", check.names = FALSE)
	weights <- weights[colnames(weights) %in% meta_data_i$Sample.ID]

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
	
	length(which(colnames(residuals)!=rownames(meta_data_i)))
	length(which(colnames(weights)!=rownames(meta_data_i)))
	length(which(colnames(residuals)!=colnames(weights)))

	## subset correct meta_data counts
	counts <- subset(meta_data_i, select = c(paste0(cell_type_i,"_geneProp")))
	colnames(counts)[1] <- "counts"
	meta_data_i <- meta_data_i[,!(grepl("_geneProp", colnames(meta_data_i)))]
	meta_data_i <- cbind(meta_data_i, counts)

	## model WOS effects
	design = model.matrix(~ age_Scale + counts + Sex + Ethnicity + Race + Blood.draw.time.point + Who.Ordinal.Scale:Blood.draw.time.point, data = meta_data_i)
	## remove any columns that are all 0s
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, weights = weights, design)
	vfit <- eBayes(vfit)

	get_geno_res <- function(beta_col, name){
		betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% beta_col)]); colnames(betas)[1] <- paste0("betas_",name)
		p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% beta_col)]); colnames(p_values)[1] <- paste0("pvalues_",name)
		fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- paste0("fdrs_",name)
		out <- cbind(betas, p_values, fdrs)
		return(out)
	}

	T1_WOS <- get_geno_res("Blood.draw.time.pointT1:Who.Ordinal.Scale", "T1")
	T2_WOS <- get_geno_res("Blood.draw.time.pointT2:Who.Ordinal.Scale", "T2")

	## bind all results
	results <- cbind(T1_WOS, T2_WOS)
	results$genes <- rownames(results)

  	## CALCULATE FDRS USING PERMUTED DATA
  	permutations_T1 <- permutations_T1[order(match(rownames(permutations_T1), rownames(results))), , drop = FALSE]
  	length(which(rownames(results)!=rownames(permutations_T1)))

	perm_fdrs = permFDR(full_data = results, full_column_id = "pvalues_T1", perm_data = permutations_T1, perm_column_ids = "all", output_name = paste0("outputs/",folder,"/7_2_WOS_effects/"))
	fdrs_intermediate <- perm_fdrs$fdrs
	colnames(fdrs_intermediate)[c(8,9)] <- paste0(colnames(fdrs_intermediate)[c(8,9)],"_T1")

	## now do T2
	permutations_T2 <- permutations_T2[order(match(rownames(permutations_T2), rownames(fdrs_intermediate))), , drop = FALSE]
  	length(which(rownames(fdrs_intermediate)!=rownames(permutations_T2)))

	perm_fdrs = permFDR(full_data = fdrs_intermediate, full_column_id = "pvalues_T2", perm_data = permutations_T2, perm_column_ids = "all", output_name = paste0("outputs/",folder,"/7_2_WOS_effects/"))
	fdrs_final <- perm_fdrs$fdrs
	colnames(fdrs_final)[c(10,11)] <- paste0(colnames(fdrs_final)[c(10,11)],"_T2")

	## order
	fdrs_final <- fdrs_final[order(match(rownames(fdrs_final), rownames(results))), , drop = FALSE]

	## write results
	write.table(fdrs_final, file = paste0("outputs/",folder,"/7_2_WOS_effects/",cell_type_i,"_resultsALL_with_qvalues.txt"), quote = FALSE, sep = ",")

  	## CALCULATE STANDARD ERRORS
  	SEs_T1 = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("Blood.draw.time.pointT1:Who.Ordinal.Scale"))]); colnames(SEs_T1) <- "SEs_T1"

  	SEs_T2 = as.data.frame(sqrt(vfit$s2.post) * vfit$stdev.unscaled[, which(colnames(vfit$coefficients) %in% c("Blood.draw.time.pointT2:Who.Ordinal.Scale"))]); colnames(SEs_T2) <- "SEs_T2"

  	results_SE <- cbind(SEs_T1, SEs_T2)
  	length(which(rownames(results_SE)!=rownames(fdrs_final)))
  	write.table(results_SE, paste0("outputs/",folder,"/7_2_WOS_effects/",cell_type_i,"_results_SEs.txt"), quote = FALSE, sep = ",")

  	print(cell_type_i)
}

