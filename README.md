
## QC
Quality control for the raw reads was performed to access the sequence quality, GC content, the presence of adaptors, and the read length.

Please refer to the script: qc.sh

## Mapping
We have used BWA, a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. The BWA-MEM algorithm, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate.

Please refer to the script : BWA-MEM.sh

## Variant detecting
We used freebayes, a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.


Please refer to the script : freebayes.sh

## Genomic variant annotations and functional effect prediction
With the help of SnpEff we have performed the genetic variant annotation and functional effects. It also annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).


Please refer to the script : SnpEff.sh

## gnomAD allele frequency
We developed a custom script that helps to grep the gnomAD allele frequency from the variant of interest.


Please refer to the script : gnomAD.sh

# Rna seq analysis
## Spliced Transcripts Alignment with STAR
To determine where on the genome our reads originated from, we aligned our reads to the reference genome using STAR (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments


Please refer to the script : STAR_Index.sh & STAR.sh 
## Counting reads
HTSeq is a Python package that calculates the number of mapped reads to each gene.


Please refer to the script : HTSeq.sh
## Transcript assembly and quantification
We have assembled the RNA-Seq alignments into potential transcripts and quantified them using StringTie.

Please refer to the scriptt : StringTie.sh
## Differential expression Analysis
We have used the R package DESeq2 to analyze count data from high-throughput sequencing assays and test for differential expressions.


Please refer to the script : StringTie.sh
## Gene ontology enrichment
Gene ontology enrichment analysis was performed with the R package GOfuncR (version 1.14.0) for up-and down-regulated DEGs. GO nodes with a family-wise error rate (FWER) <0.05 were considered significantly enriched. To check for unspecific GO enrichment analysis results, the four control subjects were split into two groups, and differential gene expression and GO enrichment analyses were performed for those two control groups.

Please refer to the Pipeline: https://rpubs.com/Akhil_Velluva/GOfuncR
## Quantification of RNA splicing
The tool Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples.

Please refer to the script : Leafcutter.sh

## Somatic variants discovery
We have used Mutect2 to call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.

Please refer to the script : Mutect2.sh

## Differentially expressed genes in lipedema study
Please refer to the file : DEG.xls

