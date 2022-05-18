#!/bin/bash
###################################
dir="/raw/fastq"                  # The name of directory where you have all the raw fastq files 
samples=$(ls $dir|sed 's/_R1_001.fastq.gz//g'|sed 's/_R2_001.fastq.gz//g'|uniq) # This is only for paired end fastq files 
for i in $samples; do 
echo $i
bwa mem -M -t 46 Refrence.fna $dir$i"_R1_001.fastq.gz" $dir$i"_R2_001.fastq.gz" | samtools sort -@46 -o $dir$i"_BWA_Sort.bam"
done
