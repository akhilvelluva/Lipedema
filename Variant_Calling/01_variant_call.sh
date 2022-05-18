#!/bin/bash
input="bam"
ref="references/GRCh38.p13.genome.fa"
output="vcf"
mkdir -p $output
samples=$(ls $input | grep .bam$)
#echo $samples
for i in $samples; do
    echo "___COMMAND___"
    echo "gatk -T HaplotypeCaller -R $ref -I $input/$i -O $output/${i}.vcf.gz"
    gatk HaplotypeCaller -R $ref -I $input/$i -O $output/${i}.vcf.gz
done



