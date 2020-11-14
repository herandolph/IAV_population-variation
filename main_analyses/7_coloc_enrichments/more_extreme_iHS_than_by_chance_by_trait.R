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
folder = "7_coloc_enrichments/iHS"

## inputs
## NOTE: declare pop equal to "CEU" or "YRI" to obtain enrichments for high iHS values in each population
pop <- "CEU"
trait_group <- "autoimmune"
## NOTE: declare trait equal to one of the below options to run the enrichment for the colocalization hits for that trait
# options: ae ms ra allergy asthma lange_cd lange_ibd lange_uc eur_cd eur_ibd euc_uc tagc_asthma vyse_sle eagle
trait <- "ae"
coloc_filter <- "sum0.5_ratio0.8"
TYPE_SNP <- "gwas_loci"

## outputs 
system(paste0("mkdir -p outputs/",folder,"/by_trait/"))
results_dir <- paste0("outputs/",folder,"/by_trait/")

## read in ihs values
iHS <- read.table(paste0("inputs/",folder,"/iHS_",pop,".txt"), header = TRUE)

## subset on all snps tested in the eQTL analysis
tested_SNPs <- as.character(read.table(paste0("inputs/",folder,"/ALL_UNIQUE_snps_tested.txt"))$V1)
## subset on tested snps
tested_SNPs_iHS <- iHS[iHS$chr_pos %in% tested_SNPs, ]
tested_SNPs_iHS$abs_Std_iHS <- abs(tested_SNPs_iHS$Std.iHS)

## read in allele frequencies calculated from CEU and YRI populations
allele_freqs <- as.data.frame(fread(paste0("inputs/",folder,"/",pop,"_AF_GRCh37.txt")))

add_AFs <- function(af_df, tested_snps_df){
	af_df <- subset(af_df, select = c(chr_pos,freq_ALT))
	af_df$MAF <- ifelse(af_df$freq_ALT > 0.5, 1 - af_df$freq_ALT, af_df$freq_ALT)
	tested_snps_df <- join(tested_snps_df, af_df, by = "chr_pos", type = "inner")
	tested_snps_df <- tested_snps_df[!duplicated(tested_snps_df$chr_pos),]
	tested_snps_df$bin <- findInterval(tested_snps_df$MAF, seq(0, 0.5, by = 0.05), rightmost.closed = TRUE)
	return(tested_snps_df)
}

tested_SNPs_iHS <- add_AFs(allele_freqs, tested_SNPs_iHS)

## read in coloc hits from publicly-available data
broader_coloc <- read.table(paste0("inputs/",folder,"/",trait_group,"_gene_colocLoci.txt"), header = TRUE, sep = "\t")
broader_coloc <- subset(broader_coloc, select = c("trait","gene",paste0(TYPE_SNP)))
broader_coloc <- broader_coloc[broader_coloc$gwas_loci != "",]
colnames(broader_coloc)[3] <- "chr_pos"

## calculate the 95th percentile for iHS
iHS_95 <- as.numeric(quantile(tested_SNPs_iHS$abs_Std_iHS, probs = 0.95, na.rm = TRUE))

npermutations <- 1000
stat_distribution <- function(snp_type){

		set.seed(2020)
		## read in data
		scFLU <- as.data.frame(fread(paste0("inputs/",folder,"/",coloc_filter,"_",trait_group,"_hits.txt")))
		scFLU$V1 <- NULL
		scFLU <- subset(scFLU, select = c("trait","gene",paste0(TYPE_SNP)))
		scFLU <- scFLU[scFLU$gwas_loci != "",]
		colnames(scFLU)[3] <- "chr_pos"

		## subset on only those test in single-cell flu data
		broader_coloc <- broader_coloc[broader_coloc$chr_pos %in% tested_SNPs_iHS$chr_pos,]

		## merge with larger dataset
		broader_coloc <- subset(broader_coloc, select = c(trait, chr_pos))
		scFLU <- subset(scFLU, select = c(trait, chr_pos))
		all_traits_datasets <- rbind(broader_coloc, scFLU)

		## merge in ihs data 
		with_ihs <- join(all_traits_datasets, tested_SNPs_iHS, by = "chr_pos")
		coloc_hits_df <- with_ihs[complete.cases(with_ihs), ]
		coloc_hits_df <- subset(coloc_hits_df, select = c(trait, chr_pos, abs_Std_iHS, MAF, bin))
		coloc_hits_df$trait <- factor(coloc_hits_df$trait)

		## pick trait
		if(trait %in% coloc_hits_df$trait){

			coloc_hits_subset <- coloc_hits_df[coloc_hits_df$trait %in% trait,]
			coloc_hits_subset <- coloc_hits_subset[!duplicated(coloc_hits_subset$chr_pos),]
			coloc_hits_num <- length(unique(coloc_hits_subset$chr_pos))

			## calculate the number of unique SNPs > 95th percentile
			obs <- subset(coloc_hits_subset, abs_Std_iHS > iHS_95)
			obs <- obs[!duplicated(obs$chr_pos),]
			obs_95 <- length(unique(obs$chr_pos))

			## calculate percentage
			obs_percentage <- obs_95/coloc_hits_num

			tested_SNPs_iHS <- subset(tested_SNPs_iHS, select = c(chr_pos, abs_Std_iHS, MAF, bin))
			## calculate a null distribution
			for(i in 1:npermutations){
				if(i %% 100 == 0){print(i)}

				## sample snps by allele frequency bins 
				## get the number of snps per bin and sample that
				for(z in 1:length(table(coloc_hits_subset$bin))){

					snp_bin <- names(table(coloc_hits_subset$bin)[z])

					table_snps <- table(coloc_hits_subset$bin)
					bin_sample <- subset(tested_SNPs_iHS, bin == snp_bin)[sample(nrow(tested_SNPs_iHS[tested_SNPs_iHS$bin == snp_bin,]), as.numeric(table_snps[names(table(coloc_hits_subset$bin)) %in% snp_bin])),]

					if(z == 1){
						null_df <- bin_sample
					}else{
						null_df <- rbind(null_df, bin_sample)
					}
				}

				## calculate the number of snps in the null that are > 95th percentile
				null_hits <- subset(null_df, abs_Std_iHS > iHS_95)
				null_hits <- length(unique(null_hits$chr_pos))

				## null percent observed
				null_percentage <- null_hits/coloc_hits_num

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

			write.table(permutation_outs, paste0(results_dir,trait_group,"_",trait,"_",snp_type,"_",coloc_filter,"_",pop,"_permutation_outs.txt"), quote = FALSE, row.names = FALSE)

			## plot 
			plot <- ggplot(permutation_outs, aes(null)) +
						geom_bar() +
						theme_bw() +
						geom_vline(xintercept = obs_percentage, color = "red", linetype = "dashed") + 
						ggtitle(paste0("n unique SNPs = ",coloc_hits_num, ", perc null tests >= obs 95th perc = ", equal_or_greater/npermutations)) +
						theme(plot.title = element_text(size = 8),
							  axis.title.x = element_text(size = 8),
							  axis.title.y = element_text(size = 8),
							  axis.text.x = element_text(size = 7),
							  axis.text.y = element_text(size = 7))

			pdf(paste0(results_dir,trait_group,"_",trait,"_",snp_type,"_",coloc_filter,"_",pop,"_distribution.pdf"), width = 4, height = 3)
			print(plot)
			dev.off()
	}else{print("trait not in coloc hit list")}
}

stat_distribution(TYPE_SNP)
