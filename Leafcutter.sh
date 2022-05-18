#!/bin/bash

dir="/home/star-out/"                        # The name of directory where you have all the BAM files
stringtieout= "/home/stringtie/"             # The name of directory where you need to store the stringtie results
leafcutter= "/home/download/annotation.gtf"  # leafcutter output folder




for bamfile in `ls $dir*bam`; do echo "Converting $bamfile to $bamfile.junc"; bash bam2junc.sh $dir$bamfile $leafcutter$bamfile.junc >> $leafcutterinfo_juncfiles.txt; done

