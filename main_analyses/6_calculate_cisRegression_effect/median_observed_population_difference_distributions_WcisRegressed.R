library(dplyr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(mgsub)
library(stringr)

current = getwd()
setwd(current)
folder = "6_calculate_cisRegression_effect/median_obs_population_differences_effect_cisRegression"

## inputs
celltypes <- c("B","CD4T","CD8T","monocytes","NK")
color_list <- c("#525564","#74828F","#96C0CE","#C25B56","#cc8f8c")

## outputs
system(paste0("mkdir -p outputs/",folder,"/"))
results_dir <- paste0("outputs/",folder,"/")
plots_dir <- paste0(results_dir,"distributions/"); dir.create(plots_dir)
outputs_dir <- paste0(results_dir,"tables/"); dir.create(outputs_dir)

## read in lists of popDE genes for which there is an eQTL
popDE_WITH_eqtl_NI <- readRDS(paste0("inputs/",folder,"/popDE_genes_with_eQTLeffect_NI.rds"))
popDE_WITH_eqtl_flu <- readRDS(paste0("inputs/",folder,"/popDE_genes_with_eQTLeffect_flu.rds"))
popDE_WITH_eqtl <- c(popDE_WITH_eqtl_NI, popDE_WITH_eqtl_flu)

## read in observed popDE results
## read in mashr outputs from continuous genetic ancestry model
real_Results <- read.table(paste0("inputs/",folder,"/popDE_continuousGA_QN_combined_posteriorMeans.txt"), header = TRUE, row.names = NULL)
real_Results[which(is.na(real_Results)),][1] <- "NA"
rownames(real_Results) <- real_Results$row.names; real_Results$row.names <- NULL
colnames(real_Results) <- mgsub(colnames(real_Results), c("_combined","_T"), c("","T"))

## read in observed popDE results with the cis-regression effect taken into account
cisRegress_Results <- read.table(paste0("inputs/",folder,"/popDE_continuousGA_QN_cisRegressed_combined_posteriorMeans.txt"), header = TRUE)
colnames(cisRegress_Results) <- mgsub(colnames(cisRegress_Results), c("_combined","_T"), c("","T"))

## se function
se <- function(x) sd(x)/sqrt(length(x))

plotDistribution <- function(infection){
	
	set.seed(2020)
	## read in significant GO terms
	GO_terms = read.table(paste0("inputs/",folder,"/GO_terms_",infection,".txt"), sep = "\t", header = TRUE)
	terms <- as.character(GO_terms$GOID)
	names <- as.character(GO_terms$GOTerm)

	## for each GO term, write and plot the median observed population differences (and p-values) before and after cis-regression
	for(z in 1:length(terms)){
		if(z %% 10 == 0){print(z)}
		term <- terms[z]
		name <- names[z]

		gene_list <- as.character(GO_terms[z,]$All.Associated.Genes)
		gene_list <- gsub("\\[|\\]", "", gene_list)
		gene_list <- unlist(str_split(gene_list, ", "))

		get_pvalue <- function(type, df_Results, df_Pvalue){

			for(i in 1:length(celltypes)){

				celltype_GE <- celltypes[i]
				celltype <- mgsub(celltype_GE,c("_combined","_T"),c("","T"))

				## get popDE eQTL
				with_eqtl <- as.character(unlist(popDE_WITH_eqtl[grepl(paste0(celltype,"_",infection), names(popDE_WITH_eqtl))]))

				## get genes that are genetic, popDE + in the GO term list
				genes_of_interest <- intersect(gene_list, with_eqtl)

				if(length(genes_of_interest) > 0){
					## read in observed expresion data
					## subset observed on specific cell type
					df <- subset(df_Results, select = paste0(celltype,"_",infection))
					df$gene <- rownames(df)
					colnames(df) <- c("observed_difference","gene")

					## combine dfs to plot correlation with inner join
					popDifferences <- df[df$gene %in% genes_of_interest,]
					median_diff <- median(popDifferences$observed_difference)
					se_diff <- se(popDifferences$observed_difference)

					## to get pvalue, sample from popDE genes with an eQTL as the null (with_eqtl_df)
					obs_num <- length(genes_of_interest)
					with_eqtl_df <- df[df$gene %in% with_eqtl,]
					background <- with_eqtl_df$observed_difference

					## calculate a null distribution
					npermutations <- 1000
					for(j in 1:npermutations){

						## sample number rows equal to the number of unique genes
						null_vec <- sample(background, obs_num)

						## take the median
						null <- median(null_vec)

						if(j == 1){
							permutation_outs <- null
						}else{
							permutation_outs <- rbind(permutation_outs, null)
						}
					}

					permutation_outs <- as.data.frame(permutation_outs, stringsAsFactors = FALSE)
					colnames(permutation_outs) <- "null"

					## directional pvalue but consider the sign of the effect size in the real, observed data only
					if(type == "real"){
						median_diff_pval <- median_diff
					}else if(type == "cis_regress"){
						df_pval <- subset(df_Pvalue, select = paste0(celltype,"_",infection))
						df_pval$gene <- rownames(df_pval)
						colnames(df_pval) <- c("observed_difference","gene")
						popDifferences_pval <- df_pval[df_pval$gene %in% genes_of_interest,]
						median_diff_pval <- median(popDifferences_pval$observed_difference)
					}

					if(median_diff_pval < 0){
						pvalue <- as.numeric(table(permutation_outs$null <= median_diff)["TRUE"])/npermutations
					}else if(median_diff_pval > 0){
						pvalue <- as.numeric(table(permutation_outs$null >= median_diff)["TRUE"])/npermutations
					}

					if(i == 1){
						outs <- c(celltype, median_diff, se_diff, pvalue, length(genes_of_interest))
					}else{
						outs_i <- c(celltype, median_diff, se_diff, pvalue, length(genes_of_interest))
						outs <- rbind(outs, outs_i)
					}
				}else{print("no genes")}
			}
			if(length(genes_of_interest) > 0){return(outs)}
		}

		## get results from get_pvalue function
		real_df <- as.data.frame((get_pvalue("real", real_Results)), stringsAsFactors = FALSE)
		if(nrow(real_df) != 0){real_df$model_type <- "with_top_cis_effect"}
		cisRegress_df <-  as.data.frame((get_pvalue("cis_regress", cisRegress_Results, real_Results)), stringsAsFactors = FALSE)
		if(nrow(real_df) != 0){cisRegress_df$model_type <- "cis_regressed"}
		outs <- rbind(real_df, cisRegress_df)

		if(nrow(outs) == 10){
			outs <- as.data.frame(outs, stringsAsFactors = FALSE)
			colnames(outs) <- c("celltype","obsDiff","se","pvalue","n_genes","model_type")
			outs$obsDiff <- as.numeric(as.character(outs$obsDiff))
			outs$se <- as.numeric(as.character(outs$se))
			outs$pvalue <- as.numeric(as.character(outs$pvalue))

			## write results
			write.table(outs, paste0(outputs_dir,gsub(" ","_",name),"_",infection,"_median.txt"), quote = FALSE, row.names = FALSE)

			## plot results
			outs$model_type <- gsub("_"," ",outs$model_type)
			outs$model_type <- factor(outs$model_type, levels = c("with top cis effect","cis regressed"))
			outs$celltype_nGene <- paste0(outs$celltype,"\nn = ",outs$n_genes)

			plot <- ggplot(outs, aes(obsDiff, celltype_nGene, color = celltype)) +
						facet_grid(. ~ model_type) +
						geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
						geom_errorbar(aes(xmin = obsDiff-se, xmax = obsDiff+se), width = 0.2, position = position_dodge(0.1)) +
						geom_text(aes(label = paste0("p = ",pvalue)), hjust = 0.5, vjust = 1.75, color = "black", size = 2) +
						geom_point(size = 1) +
						scale_color_manual(values = color_list) +
						theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
							  plot.title = element_text(size = 8),
							  axis.title.x = element_text(size = 8),
							  axis.title.y = element_text(size = 8),
							  axis.text.x = element_text(size = 7),
							  axis.text.y = element_text(size = 7),
							  legend.title = element_blank(),
							  panel.grid.major = element_blank(),
							  panel.grid.minor = element_blank(),
							  legend.position = "none",
							  panel.background = element_rect(fill = "white"),
							  strip.background = element_rect(fill = NA),
							  strip.text = element_text(size = 7)) +
						xlab("median observed expression difference") +
						ylab("") +
						ggtitle(paste0(name,", ",infection))

			pdf(paste0(plots_dir,gsub(" ","_",name),"_",infection,"_median.pdf"), width = 5.5, height = 2.5)
			print(plot)
			dev.off()
		}else{print("one or more cell types with no genes")}
	}
}

plotDistribution("NI")
plotDistribution("flu")
