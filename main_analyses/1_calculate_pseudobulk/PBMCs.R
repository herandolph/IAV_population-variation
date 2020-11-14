library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(textTinyR)
library(pbapply)

current = getwd()
folder = "1_calculate_pseudobulk"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## read in single cell data (raw UMI counts stored in RNA slot)
allCells_integrated <- readRDS(file = paste0("inputs/",folder,"/mergedAllCells_withCellTypeIdents_CLEAN.rds"))
DefaultAssay(allCells_integrated) <- "RNA"

## do not consider cells labeled as "highly infected"
dat <- subset(allCells_integrated, idents = c("monocytes","infected_monocytes","CD4_T","CD8_T","B","NK","NK_high_response","DC","neutrophils","NKT"))

## get raw UMI counts
raw_data_sparse <- GetAssayData(dat, assay = "RNA", slot = "counts")
meta_data <- dat@meta.data
sample_colname <- "sample_condition"

IDs <- as.data.frame(meta_data)[, sample_colname]
unique_ID_list <- as.list(unique(IDs))

## calculate pseudobulk by summing across all cells identified in each individual,cdt pair
pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(raw_data_sparse[,IDs == x, drop = FALSE], rowSums = TRUE)}))
cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
rownames(pseudobulk) <- rownames(x = dat)

saveRDS(pseudobulk, file = paste0("outputs/",folder,"/allCellTypes_pseudobulk.rds"))