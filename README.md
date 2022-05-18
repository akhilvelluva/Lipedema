# Lipedema

## QC
Quality control for the raw reads were performed to access the sequence quality, GC content,
the presence of adaptors and the read length
Please refer the script : qc.sh

## Mapping
We have used BWA, a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. The
BWA-MEM algorithm, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate.
Please refer the script : BWA-MEM.sh

## Variant detecting
We used freebayes, a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.


Please refer the script : freebayes.sh

## Genomic variant annotations and functional effect prediction
 With the help of SnpEff we have performd the genetic variant annotation and functional effects. It also annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).


Please refer the script : SnpEff.sh

## gnomAD allele frequency
We have build a custom script which help to grep the gnomAD allele frequency from the instered varients.


Please refer the script : SnpEff.sh

# Rna seq analysis
## Spliced Transcripts Alignment with STAR
To determine where on the genome our reads originated from, we will align our reads to the reference genome using STAR (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments


Please refer the script : STAR.sh
## Counting reads
HTSeq is a Python package that calculates the number of mapped reads to each gene.


Please refer the script : HTSeq.sh
## Transcript assembly and quantification
We have assembled the RNA-Seq alignments into potential transcripts and quatified using StringTie.


Please refer the script : StringTie.sh
## Differential expression Analysis
We have used the R package DESeq2 to analyse count data from high-throughput sequencing assays and test for differential expressions.


Please refer the script : StringTie.sh
## Gene ontology enrichment
Gene ontology enrichment analysis was performed with the R package GOfuncR (version 1.14.0) for up- and down-regulated DEGs. GO nodes with a family wise error rate (FWER) <0.05 were considered significantly enriched. To check for unspecific GO enrichment analysis results, the four control subjects were split in two groups and differential gene expression and GO enrichment analyses were performed for those two control groups.

Please refer to the Pipeline: https://rpubs.com/Akhil_Velluva/GOfuncR
## Leafcutter for the splicing
## mutect for somatic variant calling.

the 837 genes we identified to be DE in lipedema
