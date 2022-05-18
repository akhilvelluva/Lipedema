#!/bin/bash
input="vcf"
ref="references/GRCh38.p13.genome.fa"
output1="snp/variants"
output2="indel/variants"
mkdir -p $output1
mkdir -p $output2
vcfs=$(ls $input | grep .vcf.gz$)

for vc in $vcfs; do
    echo "__COMMAND__"
    echo "gatk SelectVariants -R $ref -V $input/$vc -selectType SNP -o $output1/${vc}.vcf.gz"
    gatk SelectVariants -R $ref -V $input/$vc --select-type-to-include SNP -O $output1/${vc}_snp.vcf.gz
    echo "gatk SelectVariants -R $ref -V $input/$vc -selectType INDEL -o $output2/${vc}.vcf.gz"
    gatk SelectVariants -R $ref -V $input/$vc --select-type-to-include INDEL -O $output2/${vc}_indel.vcf.gz
done
     



