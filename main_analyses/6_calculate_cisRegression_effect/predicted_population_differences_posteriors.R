library(dplyr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(limma)

current = getwd()
setwd(current)
folder = "6_calculate_cisRegression_effect/predicted_population_differences"

## outputs
system(paste0("mkdir -p outputs/",folder,"/"))
results_dir <- paste0("outputs/",folder,"/")

## read in genotypes
genotypes <- read.table(paste0("inputs/",folder,"/genotypes.txt"), header = TRUE)

## read in significant popDE genes
popDE_genes <- readRDS(paste0("inputs/",folder,"/significant_popDE_genes_continuousGA_QN_lfsr0.1.rds"))

## read in meta data
meta_data <- read.table(paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt"), header = TRUE, sep = ",")
meta_data$indiv_ID <- as.factor(meta_data$indiv_ID)
meta_data$infection_status = factor(meta_data$infection_status, levels=c("NI","flu"))

## read in mashr eQTL results
eQTL_Results <- read.table(paste0("inputs/",folder,"/eQTL_combined_posteriorMeans.txt"), header = TRUE)

## read in top SNPs
topSNPs <- as.character(read.table(paste0("inputs/",folder,"/topSNP_perCondition_perGene.txt"), header = FALSE)$V1)

calculatePredicted <- function(infection, celltype){

	## subset matrix eQTL results on celltype, infection
	eQTL <- subset(eQTL_Results, select = paste0("beta_",celltype,"_",infection))
	colnames(eQTL)[1] <- "beta"
	eQTL$geneSNPs <- rownames(eQTL)
	eQTL <- separate(eQTL, geneSNPs, into = c("gene","snps"), sep = "_", extra = "merge", remove = FALSE)

	## get celltype popDE genes
	CT_popDE <- as.character(unlist(popDE_genes[grepl(paste0(celltype,"_",infection), names(popDE_genes))]))

	## subset eQTL results on top SNPs and popDE genes
	eQTL <- eQTL[eQTL$geneSNPs %in% topSNPs,]
	eQTL <- eQTL[eQTL$gene %in% CT_popDE,]

	## get effect sizes
	eQTL <- subset(eQTL, select = c(snps, gene, beta))

	## get best eQTL SNPgene
	SNPgenes <- as.character(eQTL$snps)

	## pull the genotypes for these SNPgenes across individuals
	genotypes_subset <- genotypes[rownames(genotypes) %in% SNPgenes, ]
	eQTL <- eQTL[!duplicated(eQTL$snps),]

	## change -9 to NA
	genotypes_subset[genotypes_subset == -9] <- NA

	## put genotypes and eQTL results in order
	order_vec <- as.character(eQTL$snps)
	genotypes_subset <- genotypes_subset[match(order_vec, rownames(genotypes_subset)),]
	length(which(eQTL$snps!=rownames(genotypes_subset)))

	## grab AFR/EUR binary call from meta data
	AFR_ids <- as.character(meta_data[meta_data$ethnicity %in% "AFR",]$indiv_ID)
	EUR_ids <- as.character(meta_data[meta_data$ethnicity %in% "EUR",]$indiv_ID)

	AFR_genotypes <- genotypes_subset[, colnames(genotypes_subset) %in% AFR_ids]
	EUR_genotypes <- genotypes_subset[, colnames(genotypes_subset) %in% EUR_ids]

	## multiply genotypes by effect sizes
	AFR_popValues <- sweep(AFR_genotypes, MARGIN = 1, eQTL$beta, `*`)
	EUR_popValues <- sweep(EUR_genotypes, MARGIN = 1, eQTL$beta, `*`)

	## create predicted exp matrix
	length(which(rownames(AFR_popValues)!=rownames(EUR_popValues)))
	predicted_expression_matrix <- cbind(AFR_popValues, EUR_popValues)
	predicted_expression_matrix$snps <- rownames(predicted_expression_matrix)
	predicted_expression_matrix <- join(predicted_expression_matrix, eQTL, by = "snps")
	rownames(predicted_expression_matrix) <- predicted_expression_matrix$gene
	predicted_expression_matrix$gene <- NULL
	predicted_expression_matrix$snps <- NULL
	predicted_expression_matrix$beta <- NULL
	colnames(predicted_expression_matrix) <- paste0(colnames(predicted_expression_matrix),"_",infection)

	## get correct meta data
	meta_data_i <- meta_data[rownames(meta_data) %in% colnames(predicted_expression_matrix),]
	reorder_names <- rownames(meta_data_i)
	residuals <- predicted_expression_matrix[reorder_names]

	length(which(colnames(residuals)!=rownames(meta_data_i)))

	## run model
	design = model.matrix(~ YRI_Scale, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, design)
	vfit <- eBayes(vfit)

	betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("YRI_Scale"))]); colnames(betas)[1] <- "betas"

	## write outs
	write.table(betas, paste0(results_dir,celltype,"_",infection,".txt"), quote = FALSE, row.names = TRUE)
}

calculatePredicted("NI", "B")
calculatePredicted("flu", "B")
calculatePredicted("NI", "CD4T")
calculatePredicted("flu", "CD4T")
calculatePredicted("NI", "CD8T")
calculatePredicted("flu", "CD8T")
calculatePredicted("NI", "monocytes")
calculatePredicted("flu", "monocytes")
calculatePredicted("NI", "NK")
calculatePredicted("flu", "NK")
