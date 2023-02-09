#!usr/bin/env Rscript
#Skrypt stworzony na potrzeby normalizacji danych do projektu Bioniformatyka roślin przez Julie Denkewicz (119644)

#loading the library (if not installed, run the command > BiocManager::install("DESeq2"))
library(DESeq2)

#loading data
command_arg = commandArgs(trailingOnly = TRUE)
Data = read.csv(command_arg[1], sep = "\t", skip = 1, row.names = 1)

#checking the correctness of the implemented data (not required)
#head(Data)
sample = names(Data)

cond_1 = rep("cond1",3)
cond_2 = rep("cond2", 3)
condition = factor(c(cond_1, cond_2))

colData = data.frame(samples = sample, condition = condition)
dds_matrix = DESeqDataSetFromMatrix(countData = Data, colData = colData, design = ~condition)

#Differential gene expression analysis based on the negative binomial distribution
dds = DESeq(dds_matrix)
results = results(dds)

#filtering results (not required but recommended)
result_mean = results[results$baseMean != 0,] #without baseMean = 0
result_fold = result_mean[result_mean$log2FoldChange >1 |result_mean$log2FoldChange < -1,] #change twice as large/small
result_na = result_fold[!is.na(result_fold$padj),] #without empty padj
results = result_na[result_na$padj <0.05,]

result_frame = as.data.frame(results)
#saving the file
write.csv2(result_frame, "119644_differential_gene_expression_filtered.csv")
