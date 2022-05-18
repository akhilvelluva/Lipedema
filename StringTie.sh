#!/bin/bash

dir="/home/star-out/"                        # The name of directory where you have all the BAM files
stringtieout= "/home/stringtie/"             # The name of directory where you need to store the stringtie results
gtf= "/home/download/annotation.gtf"         # Genome annotation gtf format


###################################### stringtie
######################################

ls $dir|grep -v "bam.bai"|sed 's/.bam//g'|stringtie $dir"$r".bam -o $stringtieout$r-out.gtf -G $gtf -A $stringtieout$r-gene_abund.tab -C $stringtieout$r-cov_refs.gtf -j 10 -c 5 -a 20 -p 18;done

