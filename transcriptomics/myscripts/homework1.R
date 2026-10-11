setwd("transcriptomics")
setwd("mydata")

## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  


####################################################

### Import our data

####################################################


# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

conds_subset <- subset(conds, treatment %in% c("AM", "OWA"))
conds_subset

counts_subset <- countsTableRound[, rownames(conds_subset)]

dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData=conds_subset, 
                              design= ~ treatment)
dim(dds)
# [1] 130580     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 15,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)

dds_F4 <- subset(dds, select = generation == 'F4')
dim(dds_F4)

dds_F11 <- subset(dds, select = generation == 'F11')
dim(dds_F11)

# Perform DESeq2 analysis on each generational subset
dds_F0 <- DESeq(dds_F0)

dds_F4 <- DESeq(dds_F4)

dds_F11 <- DESeq(dds_F11)

####################################################

### Check on the DE results from each DESeq 

####################################################

resultsNames(dds_F0)


res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
head(res_OWAvsAM_F0)  
summary(res_OWAvsAM_F0)

resultsNames(dds_F4)

res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
head(res_OWAvsAM_F4)  
summary(res_OWAvsAM_F4)

resultsNames(dds_F11)

res_OWAvsAM_F11 <- results(dds_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
head(res_OWAvsAM_F11)  
summary(res_OWAvsAM_F11)

