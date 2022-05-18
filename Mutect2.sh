#!/bin/bash
######################################
######################################
********************************
*********************************************Call somatic mutations using GATK4 Mutect2
ref="/home/downloads/genome.fa"                     # Refrence genome in fasta format
dir="/home/star-out/"                               # The name of directory where you have all the BAM files
SplitNCigarReads= "/home/SplitNCigarReads/"         # The name of directory where you need to store the SplitNCigarReads BAM
PON= "/home/download/PON/"                          # The name of directory where you need to store the PanelOfNormals
Mutect= "/home/download/Mutect/"                    # The name of directory where you need to store the Mutect Results
##############################################
##############################################Step 1. CreateSequenceDictionary
samtools faidx $ref
gatk --java-options "-Xmx200G" CreateSequenceDictionary -R $ref

#############################################Step 2 MarkDuplicates (Picard)

samples=$(ls $dir| sed 's/.bam//g'|grep -v ".bai")
echo $samples
echo $dir
for i in $samples; do
picard MarkDuplicates I=$dir$i.bam O=$dir$i-MKDUP.bam M=$dir$i.matrix
done

#############################################Step 3 SplitNCigarReads. 
#STAR assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK). So we use the GATK’s ReassignOneMappingQuality read filter to reassign all good alignments to the default value of 60.

samples=$(ls $dir|grep "MKDUP.bam"| sed 's/.bam//g') #Grep only the marked duplicated bam files 
echo $samples
echo $dir
for i in $samples; do
gatk --java-options "-Xmx150G" SplitNCigarReads -R $r -I $dir$i.bam -O $SplitNCigarReads$i-SplitNCigarReads.bam
done
############################################# Step 4 CreateSomaticPanelOfNormals
############################################# RUN only the Normal SplitNCigarReads sample files
 
samples=$(ls $SplitNCigarReads|egrep -w "MKDUP-SplitNCigarReads"|grep "LipoControl"| grep -v ".bai"|sed 's/.bam//g') 
echo $samples
echo $SplitNCigarReads
for i in $samples; do
gatk --java-options "-Xmx150G" Mutect2 -R $ref -I $SplitNCigarReads$i.bam -O $PON$i-PON.vcf.gz
done

############################################# Step 5. Call somatic mutations (Mutect2)
############################################# Split sample into control and cancer 
Csamples=$(ls $SplitNCigarReads|egrep -w "MKDUP-SplitNCigarReads"|grep "LipoControl"| grep -v ".bai"|sed 's/.bam//g')
Lsamples=$(ls $SplitNCigarReads|egrep -w "MKDUP-SplitNCigarReads"|grep "Lipoma"| grep -v ".bai"|sed 's/.bam//g')
pon=$(ls $PON|grep "MKDUP-SplitNCigarReads-PON.vcf.gz"|grep -v ".tbi"|grep -v ".stats")

for p in $pon;do for l in $Lsamples; do for c in $Csamples; do
gatk --java-options "-Xmx150G" Mutect2 -R $ref -I $dir$l.bam -I $dir$c.bam -pon $pondir$p -normal `echo $c|sed 's/-MKDUP-SplitNCigarReads//g'` -O $Mutect`echo $c|sed 's/-MKDUP-SplitNCigarReads//g'`.Mutect2.vcf.gz 
done 

#############################################Step 6. FilterMutectCalls

samples=$(ls $Mutect| grep "Mutect2.vcf.gz"|grep -v ".stats"|grep -v ".tbi")
echo $samples
#echo $dir
for i in $samples; do
gatk --java-options "-Xmx200G" FilterMutectCalls  -R $ref -V $dir$i -O $Mutect`echo $i|sed 's/.Mutect2.vcf.gz/_Filtered_Mutect2.vcf.gz/g'`
done


