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
## Genomic variant annotations and functional effect prediction
 With the help of SnpEff, Genetic variant annotation and functional effect were predicted. It also annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).
## gnomAD allele frequency
We have build a custom script which help to grep the gnomAD allele frequency

# Rna seq analysis
mapping with star
htseq for gene counts
potentially strelka for the tpms
DESeq/or what Linnaeus used, go enrichment
Leafcutter for the splicing
mutect for somatic variant calling.

the 837 genes we identified to be DE in lipedema
