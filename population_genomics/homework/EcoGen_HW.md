# Population Genomics Homework

**By: Emma Shaw**

------------------------------------------------------------------------

### Overview:

------------------------------------------------------------------------

### PCA Generating & LD Pruning:

This code was run in bash in the terminal. the module plink v1.9.0-b.7.7 was used for Linkage Disequilibrium pruning on the cleaned VCF file. An LD threshold of 0.5 was used. After LD pruning,983923 SNPs were kept in our VCF, and 837827 SNPs were removed. A PCA was generated for the remaining SNPs after pruning, and the eigenvalues for the first 20 PCs were recorded.

```         
# Define directories 
DATA_DIR=/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/mydata
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz
OUT_DIR=/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework

# check the variable holds what you expect, and that the file is really there
echo $VCF
ls -lh $VCF

cd $DATA_DIR
pwd ##check you are in the right place

module load plink
plink --version   # using plink v1.9.0-b.7.7

plink --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.5 \
    --out $OUT_DIR/Lg_LD_0.5
    
# how many SNPs did we keep and remove?
wc -l $OUT_DIR/Lg_LD_0.5.prune.in $OUT_DIR/Lg_LD_0.5.prune.out
# 983923 SNPs were kept
# 837827 SNPs were removed

# make pca with the SNPs that were kept 
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract $OUT_DIR/Lg_LD_0.5.prune.in \
    --make-bed --pca --out $OUT_DIR/Lg_LD_0.5
    
# first 3 individuals, first 4 PCs (cut -f keeps only the columns we ask for)
head -n 3 $OUT_DIR/Lg_LD_0.5.eigenvec | cut -d " " -f 1-6
cat $OUT_DIR/Lg_LD_0.5.eigenval

# list of variance eigenvalues for first 20 PCs
#6.24863
#4.15376
#2.05948
#1.55805
#1.37952
#1.33352
#1.30756
#1.28463
#1.27861
#1.23942
#1.22941
#1.20565
#1.19107
#1.16842
#1.16235
#1.15438
#1.14881
#1.13823
#1.13256
#1.12893


```

------------------------------------------------------------------------

### Checking Variance Explained by PCs:

This R code was run in Rstudio. Percent variance explained was calculated for each PC (1-20) using its eigenvalue. After looking that the graph generated below, it was determined that PCs 1,2,and 3 were of interest using the broken elbow test.

```         
#load libraries
library(tidyverse)

setwd("/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework") 

# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.5.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.5.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so drop column 2
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)


pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

OUTDIR="/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.5.png"), width = 6, height = 4)
```

**Image generated:**

![](variance_0.5.png){width="430"}

------------------------------------------------------------------------

### Plotting PCA:

This R code was run in Rstudio. PCA graphs for PC1 vs PC2 and PC2 vs PC3 were generated. SItes were ordered by latituded, and coloring was done by site.

```         
# Load in metadata and join with pca results
meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")
pc3_lab <- paste0("PC3 (", round(pve$pve[3], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()


ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.5_pca1v2.png"), width = 8, height = 6)

pca_plot <- ggplot(pca, aes(PC2, PC3, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc2_lab, y = pc3_lab, colour = "Site (north to south)") +
  theme_bw()


ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.5_pca2v3.png"), width = 8, height = 6)

```

------------------------------------------------------------------------

### Admixture Prep: 

This bash code was run in terminal. A copy was made of the original .bim file before reformatting the .bim file to be compatible with ADMIXTURE.

```         
# keep a copy of the original .bim file (only run this line once! see below)
cp Lg_LD_0.5.bim Lg_LD_0.5.bim.orig

# read the original, replace column 1 (the chromosome) with 0, and write the new .bim
# OFS="\t" keeps the columns tab-separated
awk 'BEGIN{OFS="\t"} {$1="0"; print}' Lg_LD_0.5.bim.orig > Lg_LD_0.5.bim

# compare before and after
head -n 3 Lg_LD_0.5.bim.orig
head -n 3 Lg_LD_0.5.bim
```

------------------------------------------------------------------------

### ADMIXTURE sbatch job:

ADMIXTURE was run for K=1-10. Although, only tests for K=1-6 were complete after the 24 hours allotted for the job.

```         
#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=eshaw7@uvm.edu

cd /users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework
module load admixture

INPUT="Lg_LD_0.5.bed"
THREADS=${SLURM_CPUS_PER_TASK}   # use as many threads as CPUs we asked for

rm -f cv_errors.tsv
for K in $(seq 1 10); do
    admixture --cv -j${THREADS} "${INPUT}" ${K} | tee "admixture_K${K}.out"
    CV=$(grep "CV error" "admixture_K${K}.out" | awk '{print $NF}')
    echo -e "${K}\t${CV}" >> cv_errors.tsv
done
```

------------------------------------------------------------------------

------------------------------------------------------------------------

### Checking CV Error:

```         
library(tidyverse)

setwd("/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework")   # the folder with your plink output

cv <- read_tsv("cv_errors.tsv", col_names = c("K", "CV_error"))

cv_plot <- ggplot(cv, aes(K, CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K (number of ancestral populations)", y = "Cross-validation error") +
  theme_bw()

OUTDIR="/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework"

ggsave(cv_plot, file = paste0(OUTDIR, "/cv_error.png"), width = 6, height = 4)

```

------------------------------------------------------------------------

### Plotting ADMIXTURE:

```         

#Plotting ancestry Proportions
fam  <- read_table("Lg_LD_0.5.fam", col_names = FALSE)
meta <- read_csv("sample_metadata.csv")

# a function that reads the .Q file for one value of K and reshapes it
# into one row per individual per cluster
read_q <- function(K) {
  q <- read_table(paste0("Lg_LD_0.5.", K, ".Q"), col_names = FALSE)
  names(q) <- paste0("Cluster", 1:K)
  q %>%
    mutate(ID = fam$X2, K = K) %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry")
}

# map_dfr() runs read_q() for each K and binds the results together row-wise
# (change 2:5 to look at other values of K)
q_all <- map_dfr(1:5, read_q) %>%
  left_join(meta, by = "ID") %>%
  # order sites north to south, as in the PCA plots
  mutate(Site = fct_reorder(Site, -Latitude),
         K = paste0("K = ", K))

# check every individual matched the metadata: this should print 0
sum(is.na(q_all$Site))

admix_plot <- ggplot(q_all, aes(ID, Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ Site, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x     = element_blank(),         # too many sample names to read
        panel.grid      = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),      # thin gap between sites
        panel.spacing.y = unit(0.8, "lines"),      # gap between values of K
        strip.text.x    = element_text(angle = 90, hjust = 0),
        legend.position = "none")

OUTDIR="/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/homework"

ggsave(admix_plot, file = paste0(OUTDIR, "/admixture_K1-5.png"), width = 12, height = 7)
```
