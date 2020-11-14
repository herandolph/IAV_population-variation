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
coloc_filter <- "sum0.5_ratio0.8"
TYPE_SNP <- "gwas_loci"

## outputs 
system(paste0("mkdir -p outputs/",folder,"/"))
results_dir <- paste0("outputs/",folder,"/")

## read in ihs values
iHS <- read.table(paste0("inputs/",folder,"/iHS_",pop,".txt"), header = TRUE)

## subset on all snps tested in the eQTL analysis
tested_SNPs <- as.character(read.table(paste0("inputs/",folder,"/ALL_UNIQUE_snps_tested.txt"))$V1)
## subset on tested snps
tested_SNPs_iHS <- iHS[iHS$chr_pos %in% tested_SNPs, ]
tested_SNPs_iHS$abs_Std_iHS <- abs(tested_SNPs_iHS$Std.iHS)

## read in allele frequencies calculated from CEU and YRI populations
allele_freqs <- as.data.frame(fread(paste0("inputs/",folder,"/",pop,"_AF_GRCh37.txt")))

## read in list of snps that are not independent by pairwise LD
nonIndependent <- read.table(paste0("inputs/",folder,"/nonindependent_",trait_group,"_",TYPE_SNP,"_",coloc_filter,".txt"), header = TRUE)

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
broader_coloc <- subset(broader_coloc, select = c("gene",paste0(TYPE_SNP)))
broader_coloc <- broader_coloc[broader_coloc$gwas_loci != "",]
colnames(broader_coloc)[2] <- "chr_pos"

## calculate the 95th percentile for iHS
iHS_95 <- as.numeric(quantile(tested_SNPs_iHS$abs_Std_iHS, probs = 0.95, na.rm = TRUE))

## read in chr1 data for null sampling
chr1 <- read.table(paste0("inputs/",folder,"/formatted_for_sampling_1000G_chr1.txt"), header = TRUE)
chr1 <- merge(chr1, tested_SNPs_iHS, by.x = "chr_pos_A", by.y = "chr_pos")

npermutations <- 1000
stat_distribution <- function(snp_type){

		set.seed(2020)
		## read in data
		scFLU <- as.data.frame(fread(paste0("inputs/",folder,"/",coloc_filter,"_",trait_group,"_hits.txt")))
		scFLU$V1 <- NULL
		scFLU <- subset(scFLU, select = c("gene",paste0(TYPE_SNP)))
		scFLU <- scFLU[scFLU$gwas_loci != "",]
		colnames(scFLU)[2] <- "chr_pos"

		## subset on only those test in single-cell flu data
		broader_coloc <- broader_coloc[broader_coloc$chr_pos %in% tested_SNPs_iHS$chr_pos,]

		## merge with larger dataset
		broader_coloc <- subset(broader_coloc, select = c(gene, chr_pos))
		scFLU <- subset(scFLU, select = c(gene, chr_pos))
		all_traits_datasets <- rbind(broader_coloc, scFLU)

		## merge in ihs data 
		with_ihs <- join(all_traits_datasets, tested_SNPs_iHS, by = "chr_pos")
		coloc_hits_df <- with_ihs[complete.cases(with_ihs), ]
		coloc_hits_df <- coloc_hits_df[!duplicated(coloc_hits_df$chr_pos),]
		coloc_hits_df <- subset(coloc_hits_df, select = c(gene, chr_pos, abs_Std_iHS, MAF, bin))

		## subset only on independent coloc hits -- ie remove those snps that have an r2 > 0.8 for each other
		independent_hits <- coloc_hits_df[!coloc_hits_df$chr_pos %in% nonIndependent$chr_pos,]
		non_independent_hits <- coloc_hits_df[coloc_hits_df$chr_pos %in% nonIndependent$chr_pos,]

		## for each unique gene, pick the snp with the highest selection statistic
		list_sample <- list()
		for(i in 1:length(unique(non_independent_hits$gene))){
			gene_i <- as.character(unique(non_independent_hits$gene)[i])

			df_i <- non_independent_hits[non_independent_hits$gene %in% gene_i, ]
			df_i <- df_i[which(df_i$abs_Std_iHS == max(df_i$abs_Std_iHS)), ]
			df_i <- df_i[sample(nrow(df_i), 1), ]
			list_sample[[i]] <- as.character(df_i$chr_pos)
		}

		## collapse list and add back into independent loci data
		list_sample <- Reduce(union, list_sample)
		non_independent_keep <- non_independent_hits[non_independent_hits$chr_pos %in% list_sample,]
		coloc_hits_df <- rbind(independent_hits, non_independent_keep)

		coloc_hits <- length(unique(coloc_hits_df$chr_pos))

		## calculate the number of unique SNPs > 95th percentile
		obs <- subset(coloc_hits_df, abs_Std_iHS > iHS_95)
		obs <- obs[!duplicated(obs$chr_pos),]
		obs_95 <- length(unique(obs$chr_pos))

		## calculate percentage
		obs_percentage <- obs_95/coloc_hits

		tested_SNPs_iHS <- subset(tested_SNPs_iHS, select = c(chr_pos, abs_Std_iHS, MAF, bin))
		## calculate a null distribution
		for(i in 1:npermutations){
			if(i %% 100 == 0){print(i)}

			## sample snps by allele frequency bins 
			## get the number of snps per bin and sample that
			## do this for independent hits first
			for(z in 1:length(table(independent_hits$bin))){

				snp_bin <- names(table(independent_hits$bin)[z])

				table_snps <- table(independent_hits$bin)
				bin_sample <- subset(tested_SNPs_iHS, bin == snp_bin)[sample(nrow(tested_SNPs_iHS[tested_SNPs_iHS$bin == snp_bin,]), as.numeric(table_snps[names(table(independent_hits$bin)) %in% snp_bin])),]

				if(z == 1){
					null_df_indep <- bin_sample
				}else{
					null_df_indep <- rbind(null_df_indep, bin_sample)
				}
			}

			## now for linked hits
			for(z in 1:length(unique(non_independent_hits$gene))){

				gene_z <- as.character(unique(non_independent_hits$gene)[z])
				df_z <- non_independent_hits[non_independent_hits$gene %in% gene_z, ]
				df_z_hit <- df_z[which(df_z$abs_Std_iHS == max(df_z$abs_Std_iHS)), ]
				snp_bin <- df_z_hit[sample(nrow(df_z_hit), 1), ]$bin

				subset_chr1 <- chr1[chr1$bin == snp_bin,]
				## subset on instances in chr1 where the number of snps in LD is the same as the number observed
				names_numSnp_match <- names(table(subset_chr1$chr_pos_A)[table(subset_chr1$chr_pos_A) == nrow(df_z)])
				subset_chr1 <- subset_chr1[subset_chr1$chr_pos_A %in% names_numSnp_match,]

				## now, from these, sample one snp set randomly
				sample_snp_chr1 <- as.character(subset_chr1[sample(nrow(subset_chr1), 1), ]$chr_pos_A)
				sample_df_chr1 <- subset_chr1[subset_chr1$chr_pos_A %in% sample_snp_chr1,]

				## pick the highest selection stat from these
				sample_df_chr1 <- sample_df_chr1[which(sample_df_chr1$abs_Std_iHS == max(sample_df_chr1$abs_Std_iHS)), ]
				sample_df_chr1 <- sample_df_chr1[sample(nrow(sample_df_chr1), 1), ]
				snp_to_pick <- as.character(sample_df_chr1$chr_pos_A)

				## grab from tested snps df
				bin_sample <- tested_SNPs_iHS[tested_SNPs_iHS$chr_pos %in% snp_to_pick,]

				if(z == 1){
					null_df_linked <- bin_sample
				}else{
					null_df_linked <- rbind(null_df_linked, bin_sample)
				}
			}

			## combine the linked and non-linked hits data
			null_df <- rbind(null_df_indep, null_df_linked)

			## calculate the number of snps in the null that are > 95th percentile
			null_hits <- subset(null_df, abs_Std_iHS > iHS_95)
			null_hits <- length(unique(null_hits$chr_pos))

			## null percent observed
			null_percentage <- null_hits/coloc_hits

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

		write.table(permutation_outs, paste0(results_dir,trait_group,"_",snp_type,"_",coloc_filter,"_",pop,"_permutation_outs.txt"), quote = FALSE, row.names = FALSE)

		## plot 
		plot <- ggplot(permutation_outs, aes(null)) +
					geom_bar() +
					theme_bw() +
					geom_vline(xintercept = obs_percentage, color = "red", linetype = "dashed") + 
					ggtitle(paste0("n unique SNPs = ",coloc_hits, ", perc null tests >= obs 95th perc = ", equal_or_greater/npermutations)) +
					theme(plot.title = element_text(size = 8),
						  axis.title.x = element_text(size = 8),
						  axis.title.y = element_text(size = 8),
						  axis.text.x = element_text(size = 7),
						  axis.text.y = element_text(size = 7))

		pdf(paste0(results_dir,trait_group,"_",snp_type,"_",coloc_filter,"_",pop,"_distribution.pdf"), width = 4, height = 3)
		print(plot)
		dev.off()
		print(trait_group)
}

stat_distribution(TYPE_SNP)
