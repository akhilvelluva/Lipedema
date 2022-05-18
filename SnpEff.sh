#!/bin/bash
###################################### variant annotations
dir="/raw/fastq/bam/vcf/"                  # The name of directory where you have all vcf files

###################################### snpEff Building databases

java -jar snpEff.jar download -v GRCh37.75 # Example to install the SnpEff database for the human genome

###################################### 
ls $dir|grep vcf|sed 's/.vcf//g'| while read r;do java -Xmx5024m -jar snpEff.jar database-name $dir"$r".vcf >$dir"$r"-ANN.vcf ;done

