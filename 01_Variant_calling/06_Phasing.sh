#!/bin/bash

#Phasing one chromosome at a time
VCF='Chr1.mac3.Q30.lDP.dec.dup.snps.mnDP.miss0.95.pop.AB.lmaf.vcf.gz'  # Change this to your input VCF file path
Chr='Chr1'         # Change this to your desired chromosome prefix (e.g., Chr1)

#Step 1: Phase the VCF using Beagle
java -Xmx100g -jar beagle/beagle5.jar gt=${VCF} \
                   out=${Chr}_phased gp=true burnin=10 iterations=40 impute=false nthreads=40 chrom=${Chr}

#Step 2: Index the phased VCF output
bcftools index --threads 8 ${Chr}_phased.vcf.gz