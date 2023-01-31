setwd("/coursedata/users/leey17")

library(edgeR)
library(dplyr)


####### NO TREATMENT vs 1 HOUR ##########

files <- data.frame(files = c(list.files("/coursedata/project/RNA-seq-data/htseq_count_results", pattern = c("noTreatment_0hour"), full.names = F), 
                              list.files("/coursedata/project/RNA-seq-data/htseq_count_results", pattern = c("dexamethasone_1hour"), full.names = F)))




# Read in the htseq-count files using 'readDGE()' function and add the group info simultaneously   
y <- readDGE(files, path = "/coursedata/project/RNA-seq-data/htseq_count_results",
             group = factor(c(1,1,1,2,2,2,2)))


# You dont need to do anything here!
# You will get this warning:"Meta tags detected: __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique" from previous run
# This is because the special counters from htseq-count are also included in the result files.  
# Remove the Meta tags (i.e. special counters) using the following command
meta_tags <- grep("__", rownames(y))
y <- y[-meta_tags,]




# We will now perform RPKM filtering. 
# For that you need to first calculate the RPKM values for each gene in each replicate of a sample (i.e. each element of the y$counts matrix).
# Use the equation in the lecture 4 (Part 1) slides to calculate the RPKM value. 
# This calculation will take a while.

geneInfo <- read.delim("/coursedata/project/RNA-seq-data/EnsemblGeneInfo.txt", header = T) 
# use the 'gene end base' and 'gene start base' in this file to calculate the gene lengths: gene_end_base-gene_start_base+1

ensemblIds <- rownames(y)
rpkm_values <- matrix(nrow = nrow(y), ncol = ncol(y))
rownames(rpkm_values) <- ensemblIds
colnames(rpkm_values) <- y$samples$file
l <- geneInfo$Gene_End_Base - geneInfo$Gene_Start_Base + 1 

for(i in 1:nrow(rpkm_values)){
  for(j in 1:ncol(rpkm_values)){
    rpkm_values[i,j] = 10^9*y$counts[i,j]/(y[["samples"]][["lib.size"]][j]*l[i])
  }
 # Calculate RPKM values
}

# Identify the genes that have max(mean(RPKMs of 0 hour samples), mean(RPKMs of 1 hour samples)) > 0.5
# You don't need to code anything here! The following code will filter out the genes for you. 
# In the 3rd column "TRUE" indicates retaining the gene and "FALSE" indicates filtering out.
# Make sure you check the mean_rpkm object and see if the "TRUE" and "FALSE" are correct based on the second column of mean_rpkm 
mean_rpkm <- c()
for(i in 1:nrow(rpkm_values)){
  max_mean_rpkm <- max(mean(as.numeric(as.character(rpkm_values[i,1:3]))), mean(as.numeric(as.character(rpkm_values[i,4:7]))))
  if(max_mean_rpkm > 0.5){
    mean_rpkm <- rbind(mean_rpkm, cbind(rownames(rpkm_values)[i], max_mean_rpkm, TRUE))
  }else{
    mean_rpkm <- rbind(mean_rpkm, cbind(rownames(rpkm_values)[i], max_mean_rpkm, FALSE))
  }
}
keep <- as.logical(mean_rpkm[,3])
names(keep) <- mean_rpkm[,1]
y <- y[keep, , keep.lib.sizes=FALSE]




# Normalize the data using TMM normalization 
y <- calcNormFactors(y)




# You dont need to do anything here!
# Define design matrix
design <- model.matrix(~y$samples$group)
print(design) ### Check the design




# Perform trended dispersion estimation
y <- estimateGLMTrendedDisp(y,design)




# We will be performing likelihood ratio test (NOT quasi-likelihood F-test). So, choose the following glm* functions carefully!
# Fit a negative binomial generalized log-linear model to the read counts for each gene
fit <- glmQLFit(y,design)




# Perform likelihood ratio tests 
lrt <-glmLRT(fit)




# You dont need to do anything here!
# Extract the LRT results of all genes 
results <- data.frame(topTags(lrt, n = nrow(y)))
results <- data.frame(Ensembl_Gene_IDs = rownames(results), results)
results$Ensembl_Gene_IDs <- as.character(results$Ensembl_Gene_IDs)




# You dont need to do anything here!
# Annotate Ensembl gene IDs with gene names from EnsemblGeneInfo.txt file
geneInfo$Ensembl_Gene_IDs <- as.character(geneInfo$Ensembl_Gene_IDs)
results_new <- left_join(results, geneInfo, by = "Ensembl_Gene_IDs") 




# Write out the result file (Make the file name). You will need this file in Part 3
filename = "//coursedata/users/leey17/result_part1.txt"
write.table(results_new, file = filename, sep = "\t", append = T, quote = F, col.names = T, row.names = F)





# Plot a MA-plot using plotSmear() function and add this figure in your report.Highlight genes with FDR<0.05.  
pdf("MAplot.pdf")
de = decideTestsDGE(lrt,p=0.05,adjust="BH")
detags = rownames(lrt)[as.logical(de)]
plotSmear(lrt, de.tags=detags,
          main = 'MA Plot for Count Data',
          smearWidth = 0.5)
dev.off()








