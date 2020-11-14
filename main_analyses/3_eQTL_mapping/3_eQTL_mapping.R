# the following Rscript is intended to be submitted to a cluster system
# the accompanying .sbatch file and submission scripts (.sh) files are included and should be modified as necessary (currently, they are written for the slurm scheduler)
# RUN WITH 90G RAM

##########################################################
#### 1. load args & dependencies #########################
args <- commandArgs(TRUE)

## declare condition equal to "NI or "flu" to obtain cis-eQTLs for each condition
condition <- args[1]
## declare pseudobulk type
pseudobulk <- args[2]
## declare number of expression PCs to regress out
expPCs_reg <- args[3]
## declare directory for temp files
temp_dir <- args[4]
## declare output directory
out_dir <- args[5]
## name of the metadata file to use
metadata_file_name <- args[6]
## name of the counts matrix to use
counts_file_name <- args[7]
## gene positions file to use
genepos_file_name <- args[8]
## declare which cell type to use
celltype <- args[9]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

## number of iteractions for permutations that will be used for FDR correction 
iterations = 10

## load required libraries
library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(dplyr)

## source permFDR function, then re-set working directory to where you are working
current = getwd()
setwd(current)
setwd("common_functions")
source("permFDR.R")
setwd(current)

##########################################################
#### 2. load input files & clean/declare output files ####

genepos = read.table(paste0("inputs/3_eQTL_mapping/",genepos_file_name),header = TRUE, stringsAsFactors = FALSE)
colnames(genepos)[1] <- "Gene_ID"
colnames(genepos)[2] <- "chromosome"
snpspos = read.table(paste0("inputs/3_eQTL_mapping/SNP_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
gtypes = read.table(paste0("inputs/3_eQTL_mapping/genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)

## create directory structure to save outputs
system(paste0("mkdir -p outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/"))
system(paste0("mkdir -p outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/raw_results"))
## create directory for qqplots
system(paste0("mkdir -p outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/qqplots"))

## erase any previous files from and re-create the temp file
system(paste0("rm -rf outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition))
system(paste0("mkdir -p outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition))

############################################################################################################
#### 3. obtain clean input files:                                                                       ####
####   - expression tables output from pseudobulk, where first n principal components are regressed out ####
####   - covariate tables for matrixEQTL, including:                                                    ####
####     2 first genotypes PCs, age, bimodal proportion                                                 ####

## input the phenotype (for ex. the expression data) along with the corresponding meta data
metadata_whole = read.table(paste0("inputs/3_eQTL_mapping/",metadata_file_name), header = TRUE, sep = ",")
reads_whole = read.table(paste0("inputs/3_eQTL_mapping/PSEUDOCOUNTS_MATRICES/",pseudobulk,"/",celltype,"_",counts_file_name), header = TRUE, sep = ",", row.names = NULL)
reads_whole[which(is.na(reads_whole)),][1] <- "NA"
rownames(reads_whole) <- reads_whole$row.names; reads_whole$row.names <- NULL

### subset genes position file to only genes that are present in the RNA-seq data
genepos <- genepos[genepos$Gene_ID %in% rownames(reads_whole),]

## remove duplicate gene positions
genepos <- genepos %>% distinct(Gene_ID, .keep_all = TRUE)

## remove the 10 flu genes from the RNA-seq data
reads_whole <- reads_whole[rownames(reads_whole) %in% genepos$Gene_ID,]

## select only samples for which genotype data is available in the global metadata table and order samples
metadata_whole <- metadata_whole[complete.cases(metadata_whole$YRI_GRCh38), ]
rownames(metadata_whole) <- metadata_whole$infection_ID
## remove samples from meta data that we do not want from RNA-seq data
metadata_whole <- metadata_whole[which(rownames(metadata_whole) %in% colnames(reads_whole)),]
metadata_whole = metadata_whole[order(rownames(metadata_whole)),]

## subset the corresponding columns in the reads matrix and order samples and genes
# this removes samples that are in the reads_whole matrix that are not in the metadata_whole matrix
reads_whole = reads_whole[,which(colnames(reads_whole) %in% rownames(metadata_whole))]
# make sure that the samples are ordered
reads_whole = reads_whole[order(rownames(reads_whole)),order(colnames(reads_whole))]

## this is a final check to make sure the order of elements in the metadata_whole and reads_whole
length(which(rownames(metadata_whole)!=colnames(reads_whole)))

## the input data matrices are already either voomed-transformed or logCPM, so do not need to redo this
voomed_reads_whole = reads_whole

## filter metadata based on condition
metadata = metadata_whole[which(metadata_whole$infection_status==condition),]

## factor variables
metadata$infection_status = factor(metadata$infection_status)
metadata$indiv_ID = factor(metadata$indiv_ID)
metadata$ethnicity = factor(metadata$ethnicity)

## filter voomed expression per condition
voomed_reads = voomed_reads_whole[,which(colnames(voomed_reads_whole) %in% rownames(metadata))]

## check again the coherence of samples order (should be 0)
length(which(colnames(voomed_reads)!=rownames(metadata)))

## shift from sampleIDs to genotyping_IDs (there only will be one sample per genotype in every analysis)
## this removes the condition from the col names of voomed_reads and row names of metadata
colnames(voomed_reads)=metadata$indiv_ID
rownames(metadata)=metadata$indiv_ID

## recover alphabetical order with the new IDs
voomed_reads = voomed_reads[,order(colnames(voomed_reads))]
metadata = metadata[order(metadata$indiv_ID),]
## check to make sure length = 0
length(which(colnames(voomed_reads)!=rownames(metadata)))


### build matrixEQTL input -- expression tables: regressing out first n PCs from voomed_reads
## empirically remove PCs from the phenotype data
pc_set = c(1:expPCs_reg)

## regress those out
pca_rm <- function(input_data, pc_set) {
    pca = prcomp(t(input_data), na.action = na.omit)
    new = input_data
    new = apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
    new = t(new)
    colnames(new) = colnames(input_data)
    rownames(new) = rownames(input_data)
    return(new)
}
expression = pca_rm(voomed_reads, pc_set)

## perform genotype PC analysis and clean genotype data
## set the column names of gtypes to be the samples
samples=colnames(gtypes)
gtypes_pca=data.frame(snp_id=rownames(gtypes),gtypes)

## NOTE: IF YOU RUN INTO AN ERROR: RUN snpgdsClose(genofile) TO RESET FILE OPEN/CLOSE
# snpgdsClose(genofile)

## this creates a SNP genotype dataset from the gtypes_pca matrix
snpgdsCreateGeno(paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/GDS_genotypes.gds"),
                 genmat = as.matrix(gtypes_pca[, samples]),
                 sample.id = unique(samples),
                 snp.id = gtypes_pca$snp_id,
                 snpfirstdim=TRUE)

## this command tells you the total number of samples and SNPs in the .gds file
snpgdsSummary(paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/GDS_genotypes.gds"))

## load .gds file in as genofile
genofile <- snpgdsOpen(paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/GDS_genotypes.gds"))

## perform a PCA on genotype information
pca <- snpgdsPCA(genofile)
## subset the first 3 PCs
tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],
                  PC3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)

## create covariates table
## make sure that the pcs_genotypes file is properly labeled and ordered
pcs_genotypes = tab[which(tab$sample.id %in% rownames(metadata)),]
pcs_genotypes = pcs_genotypes[order(pcs_genotypes$sample.id),]
length(which(rownames(metadata)!=pcs_genotypes$sample.id))
metadata$PC1 = pcs_genotypes$PC1
metadata$PC2 = pcs_genotypes$PC2

## subset correct meta_data bimodal proportion (by cell type)
geneProp <- subset(metadata, select = c(paste0(celltype,"_geneProp")))
colnames(geneProp)[1] <- "geneProp"
metadata <- cbind(metadata, geneProp)

covariates = t(model.matrix(~ PC1 + PC2 + age_Scale + geneProp, data = metadata))
covariates = covariates[2:nrow(covariates),]

## check to make sure that rownames are correct and match gene_pos
expression = expression[which(rownames(expression) %in% genepos$Gene_ID),]
genepos = genepos[order(genepos$Gene_ID),]

## subset the genotypes file with the individuals present in each condition
genotypes = gtypes[,which(colnames(gtypes) %in% colnames(covariates))]
genotypes = genotypes[which(rownames(genotypes) %in% snpspos$snp),]
genotypes = genotypes[,order(colnames(genotypes))]
length(which(rownames(genotypes)!=snpspos$snp))

#########################################
#### 5. check that input files match ####

## check that all inputs are in the correct order (should all be 0)
length(which(rownames(metadata)!=colnames(covariates)))
length(which(rownames(metadata)!=colnames(expression)))
length(which(rownames(metadata)!=colnames(genotypes)))
length(which(rownames(genotypes)!=snpspos$snp))

## gene-wise check -- this step removes any genes for which there is no expression data (should all be 0)
genepos_trimmed <- genepos[(genepos$Gene_ID %in% rownames(expression)),]
genepos <- genepos_trimmed
length(which(rownames(expression)!=genepos_trimmed$Gene_ID))
length(which(rownames(expression)!=genepos$Gene_ID))

## this output will tell you how many genes for which there is expression data (they should all have the same length)
length(genepos_trimmed$Gene_ID)
length(genepos$Gene_ID)
length(rownames(expression))

##############################################
#### 6. save matrixEQTL temp input files. ####

expression_file_name = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/expression.txt")
covariates_file_name = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/covariates.txt")
SNP_file_name = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",temp_dir,"/",condition,"/genotypes.txt")

## write files
write.table(genotypes,SNP_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(expression,expression_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(covariates,covariates_file_name, quote=F, sep="\t", row.names=TRUE)

## in this loop iter = 0 runs the actual analyses, iters 1 to iters 10 run permutations for FDR correction
permuted_pvalues_folder = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/raw_results/")
for(iteration in 0:iterations){
    
  ###############################################################
  #### 7. permute genotype data (only for iterations > 0) #######
  
  if(iteration > 0){
    
    cols <- colnames(genotypes)
    cols.perm <- sample(cols)
    if(iteration==1){
      random_individuals_df = data.frame(cols.perm)
    }else{
      random_individuals_df = cbind(random_individuals_df,cols.perm)
    }
    genotypes <- genotypes[,cols.perm]
    colnames(genotypes) <- cols
    write.table(genotypes, SNP_file_name, sep="\t", quote = FALSE)
  }

    ##############################################################
    #### 8. prepare & run matrix EQTL ############################
    
    ## load phenotype data
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";     # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## load covariates
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    
    ## load genotype data
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA";
    snps$fileOmitCharacters = "-9" # denote missing values
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    ## set MatrixeQTL options
    useModel = modelLINEAR
    output_file_name_cis = tempfile()
    pvOutputThreshold_cis = 1
    pvOutputThreshold = 0;
    errorCovariance = numeric()
    cisDist = 1e5
    output_file_name = tempfile()
    output_file_name_cis = tempfile()
    
    ## begin MatrixeQTL
    me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold = pvOutputThreshold,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE);

    ## create and save qqplot for each iteration
    png(filename =  paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/qqplots/qqplot",iteration,".png"))
    plot(me, pch = 16, cex = 0.7)
    dev.off()

##############################################################
#### 9. write temporal output files ##########################
    
    unlink(output_file_name_cis);
    
    if(iteration == 0){
        write.table(me$cis$eqtls, file = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/raw_results/result_original.txt"))}else{
            write.table(me$cis$eqtls, file = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/raw_results/result_permuted_",iteration,".txt"))}
    
}

################################################################################################
#### 10. write best SNP-gene associations for downstream & prepare files for FDR corrections ###

## select the top cis-SNP for each gene in true and permuted files.
for(iteration in 0:iterations)
{ 
    if(iteration == 0){
        event = read.table(paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/raw_results/result_original.txt"),header=TRUE)
        event.sort <- event[order(event[,2],event[,5]),]
        event.bestQTL <- event.sort[!duplicated(event.sort$gene),]
        event.bestQTL <- event.bestQTL[order(event.bestQTL[,5]),]
    }else{
        event = read.table(paste0(permuted_pvalues_folder,"/result_permuted_",iteration,".txt"),header=TRUE)
        event.sort <- event[order(event[,2],event[,5]),]
        event.bestQTL <- event.sort[!duplicated(event.sort$gene),]
        event.bestQTL <- event.bestQTL[order(event.bestQTL[,5]),]
    }
    if(iteration==0){
        original_best_EQTL = event.bestQTL
    }else{
        if(iteration==1)
        {   
            permuted1_best_EQTL = event.bestQTL
            # only pull out the p-value
            Permutation_Input = event.bestQTL[4]
        }else{
            Permutation_Input = cbind(Permutation_Input,event.bestQTL[4])}
    }
}

##################################################
#### 11. run FDR corrections and write outputs ###

cis_eQTL_fdrs = permFDR(full_data = original_best_EQTL, full_column_id = "pvalue", perm_data = Permutation_Input, perm_column_ids = "all", output_name = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir))
write.table(cis_eQTL_fdrs$fdrs,file = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/results_best_SNPs.txt"), quote = FALSE)

## if you want to save the permuted SNPs
permuted_cis_eQTL_fdrs = permFDR(full_data = permuted1_best_EQTL, full_column_id = "pvalue", perm_data = Permutation_Input, perm_column_ids = "all", output_name = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/permuted1"))
write.table(permuted_cis_eQTL_fdrs$fdrs, file = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/permuted1_best_SNPs.txt"), quote = FALSE)

############################################################
#### 12. create qqplot using p-values of best SNPs only ####

results_best_SNPs <- read.table(paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/results_best_SNPs.txt"))

## make and save best SNPs qqplot
png(filename = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/best_SNPs_qqplot.png"))
qqplot(-log10(permuted1_best_EQTL[,4]), -log10(results_best_SNPs[,4]), ylim= c(0, 40), xlim= c(0, 40), main ="best SNPs qqplot", xlab = "-log10(theoretical p-values)", ylab = "-log10(observed p-values)")
abline(c(0,1),col="red")
dev.off()

#####################################################################################
#### 13. can also use qvalue to calc FDRs (although not used in this manuscript) ####

## calculate qvalues using the best_SNPs files
emp_pvalues <- empPvals(stat = -log10(results_best_SNPs[,4]), stat0 = -log10(as.matrix(Permutation_Input)), pool = T)
qvalues <- qvalue(p = emp_pvalues)$qvalue

## bind qvalue output with results_best_SNPs file (assumes order is the same in both)
results_best_SNPs_with_qval <- cbind(results_best_SNPs, qvalues)
## write table with this column added
write.table(x = results_best_SNPs_with_qval, file = paste0("outputs/3_eQTL_mapping/",pseudobulk,"/",celltype,"/",condition,"/",out_dir,"/results_best_SNPs_with_qval.txt"))
