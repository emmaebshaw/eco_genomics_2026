# keep a copy of the original .bim file (only run this line once! see below)
cp Lg_LD_0.5.bim Lg_LD_0.5.bim.orig

# read the original, replace column 1 (the chromosome) with 0, and write the new .bim
# OFS="\t" keeps the columns tab-separated
awk 'BEGIN{OFS="\t"} {$1="0"; print}' Lg_LD_0.5.bim.orig > Lg_LD_0.5.bim

# compare before and after
head -n 3 Lg_LD_0.5.bim.orig
head -n 3 Lg_LD_0.5.bim

module load admixture

INPUT="Lg_LD_0.5.bed"
K_MIN=1
K_MAX=10
THREADS=1   # we can increase this when running as a batch job

