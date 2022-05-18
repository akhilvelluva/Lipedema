#!/bin/bash
input="snp/variants"
java_util="/usr/bin/java"
snpeff_util="/work/users/pz192nijo/SnpEff/snpEff/snpEff.jar"
ref="GRCh38.p13"
output="snp/annotation"
mkdir -p $output

snps=$(ls $input | grep .vcf.gz$)

for snp in $snps; do
    echo "__COMMAND__"
    $java_util -Xmx8g -jar $snpeff_util -v $ref $input/$snp > $output/annotated_${snp}.vcf
done
     



