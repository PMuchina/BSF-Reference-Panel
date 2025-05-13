#!/bin/bash

#SNP filtering tutorial: https://ddocent.com/filtering/

bcftools index --threads 8 output.vcf.gz

#step 1: filters on high missignesss (>50%), minor allele count (mac3), and Quality (Q30)
vcftools --gzvcf output.vcf.gz  \
         --max-missing 0.5 \
         --mac 3 \
         --minQ 30 \
         --recode \
         --recode-INFO-all \
         --stdout | bgzip -c --threads 8 > output.mac3.Q30.vcf.gz

#step 2: Remove individuals with high missing data (> 50%)
vcftools --gzvcf output.mac3.Q30.vcf.gz --missing-indv
mawk '!/IN/' out.imiss | cut -f5 > totalmissing
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

vcftools --gzvcf output.mac3.Q30.vcf.gz \
         --remove lowDP.indv \
         --recode \
         --recode-INFO-all \
         --stdout | bgzip -c --threads 8 > output.mac3.Q30.lDP.vcf.gz

#Step 3: Decompose biallelic block substitutions
# More info on vt here: https://genome.sph.umich.edu/wiki/Vt#Decompose_biallelic_block_substitutions
bcftools index --threads 8 output.mac3.Q30.lDP.vcf.gz

vt decompose_blocksub output.mac3.Q30.lDP.vcf.gz | \
vt uniq - | \
bgzip -c --threads 8 > output.mac3.Q30.lDP.dec.dup.vcf.gz

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.vcf.gz

#Step 4: Remove Indels and multiallelic (--min-alleles/--max-alleles 2),Kept only biallelic SNPs
bcftools view -m2 -M2 -v output.mac3.Q30.lDP.dec.dup.vcf.gz \
                     -Oz -o output.mac3.Q30.lDP.dec.dup.snps.vcf.gz --threads 8

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.snps.vcf.gz

#Step 5: Filter on meanDP,adjust based on your data.
vcftools --gzvcf output.mac3.Q30.lDP.dec.dup.snps.vcf.gz \
         --min-meanDP 5 \
         --max-meanDP 25 \
         --recode \
         --recode-INFO-all \
         --stdout | bgzip -c --threads 8 > output.mac3.Q30.lDP.dec.dup.snps.mnDP.vcf.gz

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.snps.mnDP.vcf.gz

#Step 6: Missigness, adjust based on your data and downstream analysis intended
vcftools --gzvcf output.mac3.Q30.lDP.dec.dup.snps.mnDP.vcf.gz \
         --max-missing 0.95 \
         --recode \
         --recode-INFO-all \
         --stdout | bgzip -c --threads 8 > output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.vcf.gz

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.vcf.gz

#Step 7: Pop-specific missigness (Based on the number of populations you have) (optional)
# Downdload the script here: https://www.ddocent.com/filtering/

bash pop_missing_filter.sh output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.vcf.gz popmap 0.1 3 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop

bgzip --threads 8 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.recode.vcf

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.recode.vcf.gz

#Step 8: Allelic balance, adjust based on your data
bcftools filter -i '(INFO/AB > 0.25 && INFO/AB < 0.75) || INFO/AB < 0.01' -Oz -o output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.AB.vcf.gz \
                                                                                 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.recode.vcf.gz --threads 8

bcftools index --threads 8 output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.AB.vcf.gz

#Step 9: Minor Allele Frequency (MAF) threshold to exclude monomorphic sites (which have a MAF of 0).
vcftools --gzvcf output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.AB.vcf.gz \
         --maf 0.001 \
         --recode \
         --recode-INFO-all \
         --stdout | bgzip -c --threads 8 > output.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.AB.lmaf.vcf.gz





