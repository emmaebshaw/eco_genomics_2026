
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


