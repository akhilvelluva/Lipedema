#!/bin/bash
###################################### freebayes Variant detecting
dir="/raw/fastq/bam"                  # The name of directory where you have all aligned bam files

###################################### samtools “index”

ls $dir | grep .bam|while read r;do samtools index $r;done

###################################### 
ls $dir|grep bam|grep -v bai|sed 's/.bam//g'| while read r;do freebayes -f Reference.fasta $dir"$r".bam >$dir"$r".vcf ;done

