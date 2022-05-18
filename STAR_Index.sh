#!/bin/bash
treads="1"    # Number of treads
genomeDir= "/home/STAR-GENOME/" # Dirctory to output index
genome=  "/Downloads/genome.fa" # Genome file in fasta format
gtf=  "/Downloads/genome.gtf"   # Gene transfer format, file format used to hold information about gene structure
sjdb="49"                       # read length minus 1 example  If you have 100b reads, the ideal value of --sjdbOverhang is 99
mkdir -p $genomeDir

STAR --runThreadN $treads --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbOverhang $sjdb

