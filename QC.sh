#!/bin/bash
################################################ QC
#
#
################################################ This script required two tools
################################################ fastqc
################################################ mutiqc
dir="/raw/fastq"                                # The name of directory where you have all the raw fastq files 

ls $dir| grep fastq|while read r;do  fastqc $r;done

cd $dir
multiqc .
