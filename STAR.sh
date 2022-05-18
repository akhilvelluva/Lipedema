#!/bin/bash

dir="/raw/fastq"                # The name of directory where you have all the raw fastq files 
starout= "/home/star-out/"      # The name of directory where you need to store the mapped files

samples=$(ls $dir|sed 's/_R1_001.fastq.gz//g'|sed 's/_R2_001.fastq.gz//g'|uniq) # This is only for paired end fastq files 

for i in $samples; do 
echo $i
STAR --runThreadN 50 --genomeDir /home/STAR-GENOME/ --readFilesIn $dir"$r"_R1_001.fastq.gz $dir"$r"_R2_001.fastq.gz --readFilesCommand zcat --outFileNamePrefix $starout$i --outSAMtype BAM SortedByCoordinate
done


