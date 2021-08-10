library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(resample)

current = getwd()
folder = "2_calculate_residuals_and_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/2_1_DE_infection_results/"))

celltypes <- c("B","CD4_T","CD8_T","NK_combined","monocytes_combined")
inputs_dir <- paste0("inputs/",folder,"/infection_effects/")
results_dir <- paste0("outputs/",folder,"/2_1_DE_infection_results/")

## find genes DE in at least one celltype
sigList <- list()

for(i in (1:length(celltypes))){

	celltype_i <- celltypes[i]

	## read in data
	DE_results <- read.table(paste0(inputs_dir,celltype_i,"_results_with_qvalues.txt"), header = TRUE, sep = ",", row.names = NULL)
	DE_results[which(is.na(DE_results)),][1] <- "NA"
	rownames(DE_results) <- DE_results$row.names; colnames(DE_results)[1] <- "genes"

	## subset on genes that are DE + have a |beta| > 0.5
	betas <- filter(DE_results, Fdr_SAB < 0.05 & abs(betas) > 0.5)
	sigList[[i]] <- as.character(betas$genes)
}

## get union of this list
sigDEgenes <- Reduce(union, sigList)

## get betas in a format that can used to calculate specificity score
for(i in (1:length(celltypes))){

	celltype_i <- celltypes[i]
	celltype_i_name <- gsub("_combined","",celltype_i)
	celltype_i_name <- gsub("_","",celltype_i_name)

	## read in data
	DE_results <- read.table(paste0(inputs_dir,celltype_i,"_results_with_qvalues.txt"), header = TRUE, sep = ",", row.names = NULL)
	DE_results[which(is.na(DE_results)),][1] <- "NA"
	rownames(DE_results) <- DE_results$row.names; colnames(DE_results)[1] <- "genes"

	## get results for cell type
	betas <- subset(DE_results, select = c(betas))

	## change colname for merge later
	colnames(betas)[1] <- celltype_i_name

	## assign
	betas$genes <- rownames(betas)
	assign(paste0(celltype_i_name), betas)
}

log2fcMerged <- join(B, CD4T, by = "genes", type = "inner")
log2fcMerged <- join(log2fcMerged, CD8T, by = "genes", type = "inner")
log2fcMerged <- join(log2fcMerged, NK, by = "genes", type = "inner")
log2fcMerged <- join(log2fcMerged, monocytes, by = "genes", type = "inner")
rownames(log2fcMerged) <- log2fcMerged$genes
log2fcMerged$genes <- NULL

## subset on genes that are DE with |logFC| > 0.5 in at least one CT
log2fcMerged <- log2fcMerged[rownames(log2fcMerged) %in% sigDEgenes,]

## calculate CoV by row (across cell types per gene)
specificity_score <- matrix(ncol = 2, nrow = nrow(log2fcMerged))

for(i in 1:nrow(log2fcMerged)){

	row <- as.numeric(as.character(log2fcMerged[i,]))
	gene <- rownames(log2fcMerged[i,])
	CoV <- sd(row)/mean(row)

	specificity_score[i,] <- cbind(gene, CoV)
}

specificity_score <- as.data.frame(specificity_score)
colnames(specificity_score)[2] <- "specificityScore"
colnames(specificity_score)[1] <- "genes"
specificity_score <- specificity_score[,c(2,1)]
specificity_score$specificityScore <- as.numeric(as.character(specificity_score$specificityScore))

## take the absolute value
specificity_score$specificityScore <- abs(specificity_score$specificityScore)

## rank
high_to_low <- specificity_score[order(-specificity_score$specificityScore),]
low_to_high <- specificity_score[order(specificity_score$specificityScore),]

write.table(high_to_low, paste0(results_dir,"ranked_specificityScore_CoV.txt"), row.names = FALSE, quote = FALSE)
write.table(as.character(high_to_low$genes), paste0(results_dir,"ranked_specificityScore_CoV_highToLow_list.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(as.character(low_to_high$genes), paste0(results_dir,"ranked_specificityScore_CoV_lowToHigh_list.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)

