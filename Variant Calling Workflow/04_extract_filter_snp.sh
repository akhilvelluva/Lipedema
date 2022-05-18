#!/bin/bash
input="snp/annotation"
java_util="/usr/bin/java"
snpeff_util="/work/users/pz192nijo/SnpEff/snpEff/SnpSift.jar"
perl_util="/work/users/pz192nijo/SnpEff/snpEff/scripts/vcfEffOnePerLine.pl"
ref="GRCh38.p13"
output="snp/extract"
mkdir -p $output

snps=$(ls $input | grep .vcf$)

for snp in $snps; do
    echo "__COMMAND__"
    cat  $input/$snp | $perl_util | $java_util -Xmx8g -jar $snpeff_util extractFields - "CHROM" "POS" "ID" "AF" "REF" "ALT" "ANN[*].GENE" "ANN[*].BIOTYPE" "ANN[*].IMPACT" "ANN[*].HGVS_P" "GEN[*].GT" "GEN[*].AD" "GEN[*].DP" "GEN[*].GQ" "GEN[*].PL" "DP" "QUAL" | awk '$8~/(HIGH|MODERATE)/ && $7=="protein_coding" && $9!="" && $15>9 && $16>29' - |  awk '!seen[$1,$2]++' - > $output/${snp}.txt
	
done
