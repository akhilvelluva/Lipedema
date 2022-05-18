#!/bin/bash

dir="/home/star-out/"                # The name of directory where you have all the BAM files
htseqout= "/home/htseq/"             # The name of directory where you need to store the htseq count

###################################### Samtools Index

ls $dir| grep .bam|while read r;do samtools index $r;done

###################################### Htseq-count
###################################### -s whether the data is from a strand-specific assay (default: yes)
ls $dir|grep bam|grep -v bai|sed 's/.bam//g'| while read r;do htseq-count -f bam -s no -r pos -a 30 $dir"$r".bam annotation.gff3 >$htseqout$r.txt;done


