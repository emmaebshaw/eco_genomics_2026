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

# Filter out OA and OW data from conds and counts table
conds_subset <- subset(conds, treatment %in% c("AM", "OWA"))
conds_subset

rownames(conds_subset)

counts_subset <- countsTableRound[,rownames(conds_subset)]


# make sure all sample names in conds_subset and counts_subset are the same
all(colnames(counts_subset) == rownames(conds_subset)) 

dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData=conds_subset, 
                              design= ~ treatment)
dim(dds)


# Filter out low counts (lose genes with fewer than 15 counts in 75% of samples)
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

#save these results files 
# for F0
res_OWAvsAM_F0.df <- as.data.frame(res_OWAvsAM_F0)
res_OWAvsAM_F0.df$fullID <- rownames(res_OWAvsAM_F0.df)

parts <- strsplit(res_OWAvsAM_F0.df$fullID, "::")

res_OWAvsAM_F0.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F0.df,
  "myresults/F0_OWAvsAM_hw_results.csv",
  row.names = FALSE
)

# For F4
res_OWAvsAM_F4.df <- as.data.frame(res_OWAvsAM_F4)
res_OWAvsAM_F4.df$fullID <- rownames(res_OWAvsAM_F4.df)

parts <- strsplit(res_OWAvsAM_F4.df$fullID, "::")

res_OWAvsAM_F4.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F4.df,
  "myresults/F4_OWAvsAM_hw_results.csv",
  row.names = FALSE
)
# For OA vs AM
res_OWAvsAM_F11.df <- as.data.frame(res_OWAvsAM_F11)
res_OWAvsAM_F11.df$fullID <- rownames(res_OWAvsAM_F11.df)

parts <- strsplit(res_OWAvsAM_F11.df$fullID, "::")

res_OWAvsAM_F11.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F11.df,
  "myresults/F11_OWAvsAM_hw_results.csv",
  row.names = FALSE
)

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_F0, gene="TRINITY_DN1306_c0_g1::TRINITY_DN1306_c0_g1_i3::g.7453::m.7453 ", intgroup = (c("treatment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

d <-plotCounts(dds_F4, gene="TRINITY_DN1219_c0_g1::TRINITY_DN1219_c0_g1_i6::g.6988::m.6988", intgroup = (c("treatment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


d <-plotCounts(dds_F11, gene="TRINITY_DN2155_c0_g1::TRINITY_DN2155_c0_g1_i3::g.10600::m.10600", intgroup = (c("treatment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# make Eueler plot to visualize differences in expression across generations between OWA and AM

# make lists of DEGs in OWA vs AM for each generation
# For F0
res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
res_OWAvsAM_F0 <- res_OWAvsAM_F0[!is.na(res_OWAvsAM_F0$padj),]
degs_OWAvsAM_F0 <- row.names(res_OWAvsAM_F0[res_OWAvsAM_F0$padj < 0.05,])

# For F4
res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
res_OWAvsAM_F4 <- res_OWAvsAM_F4[!is.na(res_OWAvsAM_F4$padj),]
degs_OWAvsAM_F4 <- row.names(res_OWAvsAM_F4[res_OWAvsAM_F4$padj < 0.05,])

# For F11
res_OWAvsAM_F11 <- results(dds_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
res_OWAvsAM_F11 <- res_OWAvsAM_F11[!is.na(res_OWAvsAM_F11$padj),]
degs_OWAvsAM_F11 <- row.names(res_OWAvsAM_F11[res_OWAvsAM_F11$padj < 0.05,])

library(eulerr)

length(degs_OWAvsAM_F0)  #3907
length(degs_OWAvsAM_F4) #85
length(degs_OWAvsAM_F11) #1439
# Intersections
length(intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F4))  # 36
length(intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F11))  # 233
length(intersect(degs_OWAvsAM_F4,degs_OWAvsAM_F11))  # 7

# Shared across all
intF0_F4 <- intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F4)
length(intersect(degs_OWAvsAM_F11,intF0_F4)) # 2

fit1 <- euler(c("F0" = 3907, "F4" = 85, "F11" = 1439, "F0&F4" = 36, "F0&F11" = 233, "F4&F11" = 7, "F0&F4&F11" = 2))
plot(fit1,  lty = 1:3, quantities = TRUE)

######################################
# Prepare to run GO analysis####
######################################
library(topGO)

mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

# Filter GO terms by size####
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]

# Create a function to run a topGO analysis ####
run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

# Run TopGO function for 3 different contrasts####
GO_F0 <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_hw_results.csv",
  "myresults/GO_F0_OWAvsAM_hw_BP.csv"
)

GO_F4 <- run_topGO_contrast(
  "myresults/F4_OWAvsAM_hw_results.csv",
  "myresults/GO_F4_OWAvsAM_hw_BP.csv"
)

GO_F11 <- run_topGO_contrast(
  "myresults/F11_OWAvsAM_hw_results.csv",
  "myresults/GO_F11_OWAvsAM_hw_BP.csv"
)

# Make a bubble plot for topgo categories####

# Read TopGO results
# F0
go <- read.csv(
  "myresults/GO_F0_OWAvsAM_hw_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: F0 OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
# F4
go <- read.csv(
  "myresults/GO_F4_OWAvsAM_hw_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: F4 OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
# F11
go <- read.csv(
  "myresults/GO_F11_OWAvsAM_hw_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: F11 OA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )

