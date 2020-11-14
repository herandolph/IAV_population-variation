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

celltypes <- c("CD4_T","CD8_T","B")

## calculate pseudobulk sums per sample for the main cell types
for (i in 1:length(celltypes)){
  name <- celltypes[i]

  ## read in single cell data (raw UMI counts stored in RNA slot) for each cluster
  dat <- readRDS(paste0("inputs/",folder,"/",name,"_cluster_singlets.rds"))
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

  saveRDS(pseudobulk, file = paste0("outputs/",folder,"/",name,"_pseudobulk.rds"))

  print(i)
}

## calculate pseudobulk sums per sample for the cell types in which there are multiple clusters (NK and monocytes)
monocytes <- readRDS(paste0("inputs/",folder,"/monocytes_cluster_singlets.rds"))
mono_infected <- readRDS(paste0("inputs/",folder,"/infected_monocytes_cluster_singlets.rds"))
NK <- readRDS(paste0("inputs/",folder,"/NK_cluster_singlets.rds"))
NK_high_response <- readRDS(paste0("inputs/",folder,"/NK_high_response_cluster_singlets.rds"))

get_merged <- function(df1, df2, name){
    dat <- merge(df1, y = df2)
    raw_data_sparse <- GetAssayData(dat, assay = "RNA", slot = "counts")
    meta_data <- dat@meta.data
    sample_colname <- "sample_condition"

    IDs <- as.data.frame(meta_data)[, sample_colname]
    unique_ID_list <- as.list(unique(IDs))
    pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(raw_data_sparse[,IDs == x, drop = FALSE], rowSums = TRUE)}))
    cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
    colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
    rownames(pseudobulk) <- rownames(x = dat)

    saveRDS(pseudobulk, file = paste0("outputs/",folder,"/",name,"_pseudobulk.rds"))
}

get_merged(monocytes, mono_infected, "monocytes_combined")
get_merged(NK, NK_high_response, "NK_combined")

