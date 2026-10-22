# Population Genomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Emma Shaw

------------------------------------------------------------------------

## 10.13.2026 - Popgen Day 2

-   Intro to study system Lottia gigantea
-   learned basic bash commands for navigating vcf files

**Working Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics`

**Input Files:**

`none`

**Output Files:**

`none`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

```         
# Load many modules with gcc dependencies
module load gcc
module load module_name

# view bfctools resions and record the versions you use in your methods
bcftools --version   

# access bcftools manual
bcftools -h        

# grep (global regular expression print) grabs and prints lines that #match pattern
# -v inverts the command and prints the opposite 
bcftools head file_name.vcf.gz | grep -v "##contig="

# see how the file was made and what each field means
bcftools view -h file_name.vcf.gz | grep "^##INFO"
bcftools view -h file_name.vcf.gz | grep "^##FORMAT"


# bcftools also writes every command it runs into the header (in lines containing
# Command=), so you can often reconstruct how a file was made and filtered:
bcftools view file_name.vcf.gz | grep "Command="

# The VCF also stores our sample names, which we can list with bcftools query -l:
bcftools query -l Lottia_gigantea_pl_flt3.renamed.vcf.gz

# count the samples by piping into wc -l (counts lines)
bcftools query -l file_name.vcf.gz | wc -l
```

------------------------------------------------------------------------

## 10.15.2026 - Popgen Day 2

-   learned how to submit sbatch job

**Working Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics`

**Input Files:**

`Lottia_gigantea_pl_flt3.renamed.vcf.gz`

**Output Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/missing_data`

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/depth_data`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`vcf_stats.sh`

**Code:**

```         
#!/bin/bash
#SBATCH --job-name=vcf_stats
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=30:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eshaw7@uvm.edu

module load gcc
module load vcftools

##set your path to where the data is, change to your own path!!
INDIR="/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/mydata"
##set your path to where the data will be stored, change to your own path!!
OUTDIR="/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults"

##change to directory where data is
cd $INDIR
##change out directory to where the results will go
vcftools --gzvcf Lottia_gigantea_pl_flt3.renamed.vcf.gz --missing-indv --out $OUTDIR/missing_data
vcftools --gzvcf Lottia_gigantea_pl_flt3.renamed.vcf.gz --depth --out $OUTDIR/depth_data
```

------------------------------------------------------------------------

## 10.20.2026 - Popgen Day 3

-   

**Working Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics`

**Input Files:**

`Lottia_gigantea_pl_flt3.renamed.vcf.gz`

**Output Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/missing_data`

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/depth_data`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`vcf_stats.sh`

**Code:**

```         
DATA_DIR=/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/mydata
```

------------------------------------------------------------------------

## 10.22.2026 - Popgen Day 4

-   **Genetic Diversity Can be Measured in:**

    -   heterozygosity

    -   private alleles

    -   Fst

    -   Nucleotide Diversity

    -   Inbreeding coefficient

-   **Hardy Weinberg Equilibrium**

    -   p^2^ +2pq +q^2^ =1

-   **Indbreeding Ceefficient**

    -   Ho= heterozygous individuals

    -   He= total individuals

    -   F=1-Ho/He

-   **Tajima's D**

    -   used to detect selection

    -   greater than 0 is a signature of balancing selection

    -   less than 0 shows a selective sweep

**Working Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics`

**Input Files:**

`Lottia_gigantea_pl_flt3.renamed.vcf.gz`

**Output Directory:**

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/missing_data`

`/gpfs1/home/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/myresults/depth_data`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`vcf_stats.sh`

**Code:**

```         
DATA_DIR=/users/e/s/eshaw7/projects/eco_genomics__2026/population_genomics/mydata
```
