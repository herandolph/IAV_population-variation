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
system(paste0("mkdir -p inputs/",folder,"/corrected_expression/"))
system(paste0("mkdir -p outputs/",folder,"/2_1_DE_infection_results/"))

## source function
setwd("common_functions")
source("permFDR.R")

setwd(current)
## inputs
meta_data_file <- paste0("inputs/",folder,"/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt")
pseudobulk_dir <- paste0("inputs/",folder,"/2_1_correctedExp_and_infection_DE/pseudobulk/")
celltypes <- c("B","CD4_T","CD8_T","monocytes_combined","NK_combined")

## outputs
residuals_dir <- paste0("inputs/",folder,"/corrected_expression/")
results_dir <- paste0("outputs/",folder,"/2_1_DE_infection_results/")

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

	## get average capture effect
	avg_capture_effect <- rowMeans(fit$coefficients)

	## add average capture effect back into residuals 
	corrected_expression <- apply(residuals,2,function(x){x + avg_capture_effect}) 
	weights <- v$weights
	colnames(weights) <- colnames(corrected_expression)
	rownames(weights) <- rownames(corrected_expression)

	## write capture-corrected expression and weights to model popDE and popDR later
	write.table(corrected_expression, paste0(residuals_dir,cell_type_i,"_corrected_expression.txt"), quote = FALSE, sep = ",")
	write.table(weights, paste0(residuals_dir,cell_type_i,"_weights.txt"), quote = FALSE, sep = ",")

	length(which(colnames(corrected_expression)!=rownames(meta_data_i)))

	## subset correct meta_data bimodal proportion
	counts <- subset(meta_data_i, select = c(paste0(cell_type_i,"_geneProp")))
	colnames(counts)[1] <- "counts"
	meta_data_i <- meta_data_i[,c(1:22)]
	meta_data_i <- cbind(meta_data_i, counts)

	# read in permutations
	permutations <- read.table(paste0("inputs/",folder,"/2_1_correctedExp_and_infection_DE/permutations_",cell_type_i,".txt"), header = TRUE, sep = ",", row.names = NULL)
	permutations$genes <- as.character(permutations$genes)
	permutations[which(is.na(permutations)),][1] <- "NA"

	## model infection differential expression
	design = model.matrix(~ indiv_ID + counts + infection_status*age_Scale, data = meta_data_i)
	design <- design[, colSums(design != 0) > 0]

	## remove global age effects
	design <- subset(design, select = -c(age_Scale))

	vfit <-lmFit(corrected_expression, weights = weights, design)
	vfit <- eBayes(vfit)

	## collect outputs
	betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusflu"))]); colnames(betas)[1] <- "betas"
	p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("infection_statusflu"))]); colnames(p_values)[1] <- "pvalues"
	fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- "fdrs"
	t_stats = as.data.frame(cbind(rownames(topTable(vfit, coef = "infection_statusflu", number = Inf)), topTable(vfit, coef = "infection_statusflu", number = Inf)$t))
	colnames(t_stats) <- c("genes","t_stat")

	results <- cbind(betas, p_values, fdrs)
	results$genes <- rownames(results)
	results <- join(results, t_stats, by = "genes")
	
	## match row names for permutations
	perms <- permutations[match(results$genes, permutations$genes),]
	length(which(results$genes!=perms$genes))
	rownames(results) <- results$genes; results$genes <- NULL
	rownames(perms) <- perms$genes; perms$genes <- NULL

	## calculate qvalues using a permutation-based null
	perm_fdrs = permFDR(full_data = results, full_column_id = "pvalues", perm_data = perms, perm_column_ids = "all", output_name = paste0(results_dir))
	write.table(perm_fdrs$fdrs, file = paste0(results_dir,cell_type_i,"_results_with_qvalues.txt"), quote = FALSE, sep = ",")

	corr_fdrs <- subset(perm_fdrs$fdrs)

  	if(i == 1){
  		numGenes <- dim(dge)[1]
  		celltype <- cell_type_i
  	}else{
  		numGenes[i] <- dim(dge)[1]
  		celltype[i] <- cell_type_i
  	}

  	if(i == 1){
  		numDEgenes_05 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_01 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.01 & abs(corr_fdrs[,"betas"]) > 0.5))
  		celltype <- cell_type_i
  	}else{
  		numDEgenes_05[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_01[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.01 & abs(corr_fdrs[,"betas"]) > 0.5))
  		celltype[i] <- cell_type_i
  	}

  	print(cell_type_i)
}

numDEgenes_per_CT <- as.data.frame(cbind(celltype, numDEgenes_05, numDEgenes_01))
write.table(numDEgenes_per_CT, paste0(results_dir,"number_of_DEG_by_CT_qvalues_logFCcutoff.txt"), quote = FALSE, sep = ",")
