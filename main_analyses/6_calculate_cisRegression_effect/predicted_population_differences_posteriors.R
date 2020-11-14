library(dplyr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)

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

	## write these as the predicted expression matrices
	write.table(predicted_expression_matrix, paste0(results_dir,celltype,"_",infection,"_predicted_expression.txt"), quote = FALSE)

	## take the rowMean of these values across individuals for each gene 
	get_popAvg <- function(df, name){
		popAvg <- as.data.frame(rowMeans(df, na.rm = TRUE))
		colnames(popAvg)[1] <- name
		popAvg$snps <- rownames(popAvg)
		popAvg <- popAvg[, c(2,1)]
	}

	AFR_popAvg <- get_popAvg(AFR_popValues, "AFR_popAvg")
	EUR_popAvg <- get_popAvg(EUR_popValues, "EUR_popAvg")

	## join
	popAvg <- join(AFR_popAvg, EUR_popAvg, by = "snps")

	## calculate difference (binary call)
	popAvg$difference <- popAvg$AFR_popAvg - popAvg$EUR_popAvg

	## add back in gene names
	popAvg <- join(popAvg, eQTL, by = "snps")

	## write outs
	write.table(popAvg, paste0(results_dir,celltype,"_",infection,".txt"), quote = FALSE, row.names = FALSE)
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
