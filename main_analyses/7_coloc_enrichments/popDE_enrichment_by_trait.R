library(dplyr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(data.table)
library(mgsub)

current = getwd()
setwd(current)
folder = "7_coloc_enrichments/popDE_genes"

## inputs
trait_group <- c("autoimmune")
coloc_filter <- "sum0.5_ratio0.8"
TYPE_SNP <- "gwas_loci"

## outputs 
system(paste0("mkdir -p outputs/",folder,"/by_trait/"))
results_dir <- paste0("outputs/",folder,"/by_trait/")

## read in significant popDE genes
popDE_genes <- readRDS(paste0("inputs/",folder,"/significant_popDE_genes_continuousGA_QN_lfsr0.1.rds"))
popDE_genes <- Reduce(union, popDE_genes)

## read in background genes
background_geneset <- as.character(read.table(paste0("inputs/",folder,"/ALL_UNIQUE_genes_tested.txt"), header = FALSE)$V1)

## read in coloc hits from publicly-available data
broader_coloc <- read.table(paste0("inputs/",folder,"/",trait_group,"_gene_colocLoci.txt"), header = TRUE, sep = "\t")
broader_coloc <- subset(broader_coloc, select = c("trait","gene",paste0(TYPE_SNP)))
broader_coloc <- broader_coloc[broader_coloc$gwas_loci != "",]
broader_coloc <- broader_coloc[broader_coloc$gene %in% background_geneset,]

stat_distribution <- function(snp_type, coloc_filter){

	set.seed(2020)
	for(j in 1:length(trait_groups)){
		coloc <- broader_coloc
		colnames(coloc)[2] <- "pid"

		## read in data
		scFLU <- as.data.frame(fread(paste0("inputs/",folder,"/",coloc_filter,"_",trait_group,"_hits.txt")))
		scFLU$V1 <- NULL
		scFLU <- subset(scFLU, select = c("trait","gene",paste0(TYPE_SNP)))
		scFLU <- scFLU[scFLU$gwas_loci != "",]
		colnames(scFLU)[2] <- "pid"

		## merge with larger dataset
		coloc <- subset(coloc, select = c(trait, pid))
		scFLU <- subset(scFLU, select = c(trait, pid))
		all_traits_datasets <- rbind(coloc, scFLU)

		plots_list <- list()
		for(z in 1:length(unique(all_traits_datasets$trait))){
			trait <- as.character(unique(all_traits_datasets$trait)[z])

			all_traits_subset <- all_traits_datasets[all_traits_datasets$trait %in% trait,]
			## calculate the number of unique genes that are associated with coloc hits
			obs_genes <- as.character(unique(all_traits_subset$pid))
			obs_num <- length(unique(all_traits_subset$pid))

			## popDE genes within this list
			popDE_num <- length(obs_genes[obs_genes %in% popDE_genes])

			## calculate percentage
			obs_percentage <- popDE_num/obs_num

			## calculate a null distribution
			npermutations <- 1000
			for(i in 1:npermutations){

				## sample number rows equal to the number of unique genes
				null_vec <- sample(background_geneset, obs_num)

				## calculate the number these genes that are popDE
				null_popDE <- length(null_vec[null_vec %in% popDE_genes])

				## null percent observed
				null_percentage <- null_popDE/obs_num

				if(i == 1){
					permutation_outs <- null_percentage
				}else{
					permutation_outs <- rbind(permutation_outs, null_percentage)
				}
			}

			permutation_outs <- as.data.frame(permutation_outs)
			colnames(permutation_outs) <- "null"
			## calculate a p-value
			permutation_outs$eval <- ifelse(permutation_outs$null >= obs_percentage, "TRUE", "FALSE")
			equal_or_greater <- as.numeric(table(permutation_outs$null >= obs_percentage)["TRUE"])

			write.table(permutation_outs, paste0(results_dir,"/",trait_group,"_",trait,"_",snp_type,"_",coloc_filter,"_permutation_outs.txt"), quote = FALSE, row.names = FALSE)

			## plot 
			plots_list[[z]] <- ggplot(permutation_outs, aes(null)) +
						geom_bar() +
						theme_bw() +
						geom_vline(xintercept = obs_percentage, color = "red", linetype = "dashed") + 
						ggtitle(paste0(trait,", ","n unique genes = ",obs_num, "\nperc null tests >= obs perc popDE genes = ", equal_or_greater/npermutations)) +
						theme(plot.title = element_text(size = 6),
							  axis.title.x = element_text(size = 8),
							  axis.title.y = element_text(size = 8),
							  axis.text.x = element_text(size = 7),
							  axis.text.y = element_text(size = 7))
		}

		plots <- lapply(plots_list, ggplotGrob)
		ggsave(paste0(results_dir,"/",trait_group,"_",snp_type,"_",coloc_filter,"_distribution.pdf"), marrangeGrob(plots, nrow = 4, ncol = 3), width = 7, height = 9, units = "in")
	}
}

stat_distribution(TYPE_SNP, "sum0.5_ratio0.8")

