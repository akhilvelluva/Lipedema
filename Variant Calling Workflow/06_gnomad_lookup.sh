#!/bin/bash
input="snp/extract"
java_util="/usr/bin/java"
ref="GRCh38.p13"
output="snp/gnomad"
mkdir -p $output

snps=$(ls $input | grep .txt$)

for snp in $snps; do
    echo "__COMMAND__"
    python intersect.py $input/$snp gnomad.table	
done
