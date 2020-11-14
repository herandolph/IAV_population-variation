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
system(paste0("mkdir -p outputs/",folder,"/cis_regressed/"))

## inputs
expression_dir <- paste0("inputs/",folder,"/corrected_expression_QN/")
celltypes <- c("B","CD4_T","CD8_T","monocytes_combined","NK_combined")

## outputs
results_dir <- paste0("outputs/",folder,"/cis_regressed/")

## read in meta data from individuals
meta_data <- read.table(paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt"), header = TRUE, sep = ",")
meta_data$capture <- as.factor(meta_data$capture)
meta_data$infection_status = factor(meta_data$infection_status, levels=c("NI","flu"))
meta_data <- meta_data[complete.cases(meta_data$YRI_GRCh38), ]

## read in genotypes
genotypes <- read.table(paste0("inputs/",folder,"/genotypes.txt"), header = TRUE)

## read in top SNPs
topSNPs <- read.table(paste0("inputs/",folder,"/topSNP_perCondition_perGene.txt"), header = FALSE)
topSNPs <- separate(topSNPs, V1, into = c("gene","snp"), extra = "merge", sep = "_")
## subset genotypes to pertinent ones
genotypes$snp <- rownames(genotypes)
genotypes_subset <- merge(genotypes, topSNPs, by = "snp")
rownames(genotypes_subset) <- genotypes_subset$gene; genotypes_subset$gene <- NULL; genotypes_subset$snp <- NULL
genotypes_subset[genotypes_subset == -9] <- NA

## get list of common genes to subset cell type expression data on 
common_genes <- rownames(genotypes_subset)

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

	## subset residuals and weights on common genes
	residuals <- residuals[rownames(residuals) %in% common_genes,]
	weights <- weights[rownames(weights) %in% common_genes,]

	length(which(colnames(residuals)!=rownames(meta_data_i)))
	length(which(colnames(weights)!=rownames(meta_data_i)))
	length(which(colnames(residuals)!=colnames(weights)))
	length(which(rownames(residuals)!=rownames(weights)))

	for(j in 1:length(common_genes)){
		if(j%%1000 == 0) print(j)
		gene <- common_genes[j]

		## add genotype into meta data
		gene_genotype <- as.data.frame(t(genotypes_subset[rownames(genotypes_subset) %in% gene,]))
		colnames(gene_genotype)[1] <- "genotype"
		gene_genotype$indiv_ID <- rownames(gene_genotype)
		meta_data_loop <- join(meta_data_i, gene_genotype, by = "indiv_ID")
		meta_data_loop$genotype <- as.numeric(meta_data_loop$genotype)
		rownames(meta_data_loop) <- meta_data_loop$infection_ID

		## model the continuous genetic ancestry effect while taking into account the genotype of the top cis SNP
		design = model.matrix(~ age_Scale + counts + genotype + infection_status + YRI_Scale:infection_status, data = meta_data_loop)

		## remove individuals with an NA genotype
		residuals_loop <- residuals[,colnames(residuals) %in% rownames(design)]
		weights_loop <- residuals[,colnames(weights) %in% rownames(design)]

		length(which(colnames(residuals_loop)!=rownames(design)))
		length(which(colnames(weights_loop)!=rownames(design)))
		length(which(colnames(residuals_loop)!=colnames(weights_loop)))

		vfit <- lmFit(residuals_loop, weights = weights_loop, design)
		vfit <- eBayes(vfit)

		## for each gene, store each output separately
		betas_NI = as.data.frame(vfit$coefficients[rownames(vfit$coefficients) %in% gene, which(colnames(vfit$coefficients) %in% c("infection_statusNI:YRI_Scale"))]); colnames(betas_NI)[1] <- "betas_NI"
		betas_flu = as.data.frame(vfit$coefficients[rownames(vfit$coefficients) %in% gene, which(colnames(vfit$coefficients) %in% c("infection_statusflu:YRI_Scale"))]); colnames(betas_flu)[1] <- "betas_flu"
		results_j <- cbind(betas_NI, betas_flu)

		if(j == 1){
			beta_outs <- results_j
		}else{
			beta_outs <- rbind(beta_outs, results_j)
		}

		s2_post <- as.data.frame(vfit$s2.post); colnames(s2_post)[1] <- "s2_post"
		rownames(s2_post) <- rownames(vfit$stdev.unscaled)
		SEs_NI = as.data.frame(sqrt(s2_post[rownames(s2_post) %in% gene,]) * vfit$stdev.unscaled[rownames(vfit$coefficients) %in% gene, which(colnames(vfit$coefficients) %in% c("infection_statusNI:YRI_Scale"))]); colnames(SEs_NI) <- "SEs_NI"
		SEs_flu = as.data.frame(sqrt(s2_post[rownames(s2_post) %in% gene,]) * vfit$stdev.unscaled[rownames(vfit$coefficients) %in% gene, which(colnames(vfit$coefficients) %in% c("infection_statusflu:YRI_Scale"))]); colnames(SEs_flu) <- "SEs_flu"
		results_SE_j <- cbind(SEs_NI, SEs_flu)

		if(j == 1){
			SE_outs <- results_SE_j
		}else{
			SE_outs <- rbind(SE_outs, results_SE_j)
		}

		pvalues_NI = as.data.frame(vfit$p.value[rownames(vfit$p.value) %in% gene, which(colnames(vfit$p.value) %in% c("infection_statusNI:YRI_Scale"))]); colnames(pvalues_NI)[1] <- "pvalues_NI"
		pvalues_flu = as.data.frame(vfit$p.value[rownames(vfit$p.value) %in% gene, which(colnames(vfit$p.value) %in% c("infection_statusflu:YRI_Scale"))]); colnames(pvalues_flu)[1] <- "pvalues_flu"
		results_pvals_j <- cbind(pvalues_NI, pvalues_flu)

		if(j == 1){
			pvals_outs <- results_pvals_j
		}else{
			pvals_outs <- rbind(pvals_outs, results_pvals_j)
		}
	}

	rownames(beta_outs) <- common_genes
	rownames(SE_outs) <- common_genes
	rownames(pvals_outs) <- common_genes

	write.table(beta_outs, paste0(results_dir,cell_type_i,"_resultsALL.txt"), quote = FALSE, sep = ",")
	write.table(SE_outs, paste0(results_dir,cell_type_i,"_results_SEs.txt"), quote = FALSE, sep = ",")
	write.table(pvals_outs, paste0(results_dir,cell_type_i,"_results_pvals.txt"), quote = FALSE, sep = ",")
	print(cell_type_i)
}

