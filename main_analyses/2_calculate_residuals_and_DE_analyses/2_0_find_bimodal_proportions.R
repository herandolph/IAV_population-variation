library(mixtools)
library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrepel)
library(stringr)

current = getwd()
folder = "2_calculate_residuals_and_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/2_0_bimodal_proportions/"))

## residuals inputs
residuals_dir <- paste0("inputs/",folder,"/corrected_expression/")
results_dir <- paste0("outputs/",folder,"/2_0_bimodal_proportions/")
celltypes <- c("monocytes_combined","CD8_T","CD4_T","B","NK_combined")
thresholds <- list(c(-0.3,5), c(0.2,5), c(-0.3,2), c(0.5,5), c(2,6))

for(i in 1:length(celltypes)){
	celltype_i <- celltypes[i]
	threshold_i <- thresholds[[i]]

	## expression estimates
	residuals <- read.table(paste0(residuals_dir,celltype_i,"_corrected_expression.txt"), header = TRUE, sep = ",", row.names = NULL)
	residuals[which(is.na(residuals)),][1] <- "NA"
	rownames(residuals) <- residuals$row.names; residuals$row.names <- NULL

	## convert to long format and cbind
	residuals_m <- melt(residuals)
	expression_ests <- separate(residuals_m, variable, into = c("indivID","infection"), sep = "_") 

	## separate for NI and FLU
	expression_NI <- expression_ests[expression_ests$infection %in% "NI",]
	expression_FLU <- expression_ests[expression_ests$infection %in% "flu",]

	individuals <- unique(expression_ests$indivID)

	for(j in 1:length(individuals)){
		indiv_j <- individuals[j]
		NI_j <- expression_NI[which(expression_NI$indivID %in% indiv_j),]
		FLU_j <- expression_FLU[which(expression_FLU$indivID %in% indiv_j),]

		## subset
		NI_j_subs <- NI_j[which(NI_j$value > threshold_i[1] & NI_j$value < threshold_i[2]),]
		FLU_j_subs <- FLU_j[which(FLU_j$value > threshold_i[1] & FLU_j$value < threshold_i[2]),]

		## calculate density distributions
		d_NI <- density(NI_j_subs$value)
		d_FLU <- density(FLU_j_subs$value)

		## get the min peak (calculate discrete analog of second derivative -- positive at local minima + take the first of these values)
		min_NI <- d_NI$x[which(diff(sign(diff(d_NI$y)))==2)+1][1]
		min_FLU <- d_FLU$x[which(diff(sign(diff(d_FLU$y)))==2)+1][1]

		## collect how many genes fall below this threshold for each sample
		genes_prop_NI <- round(length(which(NI_j_subs$value < min_NI))/dim(NI_j_subs)[1], 5)
		genes_prop_FLU <- round(length(which(FLU_j_subs$value < min_FLU))/dim(FLU_j_subs)[1], 5)

		## put this in a table
		if(j == 1){
			out_j <- c(indiv_j, genes_prop_NI, genes_prop_FLU)
			out <- out_j
		}else{
			out_j <- c(indiv_j, genes_prop_NI, genes_prop_FLU)
			out <- rbind(out, out_j)
		}
	}

	out <- as.data.frame(out)
	colnames(out) <- c("indivID","prop_NI","prop_FLU")
	rownames(out) <- out$indivID
	out$prop_NI <- as.numeric(as.character(out$prop_NI))
	out$prop_FLU <- as.numeric(as.character(out$prop_FLU))

	if(celltype_i == "CD8_T"){
		out[out$indivID == "HMN52534",]$prop_NI <- 0.5
		out[out$indivID == "HMN52534",]$prop_FLU <- 0.6
	}else if(celltype_i == "B"){
		out[out$indivID == "HMN52533",]$prop_FLU <- 0.01
		out[out$indivID == "HMN52548",]$prop_NI <- 0
		out[out$indivID == "HMN171215",]$prop_FLU <- 0.01
		out[out$indivID == "HMN83557",]$prop_NI <- 0.4
	}else if(celltype_i == "CD4_T"){
		out[out$indivID == "HMN83573",]$prop_NI <- 0.01
		out[out$indivID == "HMN83580",]$prop_FLU <- 0.01
		out[out$indivID == "HMN52537",]$prop_FLU <- 0.03
		out[out$indivID == "HMN52547",]$prop_FLU <- 0.03
		out[out$indivID == "HMN52540",]$prop_FLU <- 0.03
		out[out$indivID == "HMN171219",]$prop_FLU <- 0.01
		out[out$indivID == "HMN52553",]$prop_NI <- 0.01
		out[out$indivID == "HMN52553",]$prop_FLU <- 0.01
	}else if(celltype_i == "monocytes_combined"){
		out[out$indivID == "HMN83553",]$prop_NI <- 0.52
		out[out$indivID == "HMN83556",]$prop_NI <- 0.01
		out[out$indivID == "HMN83561",]$prop_NI <- 0.01
		out[out$indivID == "HMN83575",]$prop_FLU <- 0.30
		out[out$indivID == "HMN83577",]$prop_NI <- 0.01
		out[out$indivID == "HMN52536",]$prop_NI <- 0.01
		out[out$indivID == "HMN52544",]$prop_NI <- 0.005
		out[out$indivID == "HMN52544",]$prop_FLU <- 0.005
		out[out$indivID == "HMN52561",]$prop_NI <- 0.01
		out[out$indivID == "HMN171218",]$prop_NI <- 0.01
		out[out$indivID == "HMN171218",]$prop_FLU <- 0.03
		out[out$indivID == "HMN171232",]$prop_NI <- 0.75
		out[out$indivID == "HMN171232",]$prop_FLU <- 0.65
		out[out$indivID == "HMN171221",]$prop_NI <- 0.02
		out[out$indivID == "HMN171238",]$prop_NI <- 0.65
		out[out$indivID == "HMN171238",]$prop_FLU <- 0.85
		out[out$indivID == "HMN171223",]$prop_NI <- 0.025
		out[out$indivID == "HMN171223",]$prop_FLU <- 0.025
	}else if(celltype_i == "NK_combined"){
		out[out$indivID == "HMN83558",]$prop_NI <- 0.70
		out[out$indivID == "HMN83562",]$prop_NI <- 0.65
		out[out$indivID == "HMN83575",]$prop_NI <- 0.80
		out[out$indivID == "HMN83575",]$prop_FLU <- 0.55
		out[out$indivID == "HMN83573",]$prop_NI <- 0.60
		out[out$indivID == "HMN52543",]$prop_NI <- 0.85
		out[out$indivID == "HMN52560",]$prop_FLU <- 0.005
		out[out$indivID == "HMN52549",]$prop_FLU <- 0.75
		out[out$indivID == "HMN171218",]$prop_NI <- 0.005
		out[out$indivID == "HMN171218",]$prop_FLU <- 0.005
	}

	write.table(out, paste0(results_dir,celltype_i,"_bimodal_proportions.txt"), quote = FALSE, row.names = FALSE)
	print(celltype_i)
}
