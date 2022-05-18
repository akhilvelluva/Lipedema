Variant Calling Workflow

01_variant_call.sh
- description: Peform variant calling for a given list of BAM files
- input: bam
- tool: gatk HaplotypeCaller
- reference: GRCh38

02_extract.sh
- description: Extract the variants per variant types i.e. [SNP] [INDEL]
- input: vcf.gz
- tool: gatk SelectVariants
- reference: GRCh38

03_annotate_snp.sh
- description: Annotate the variants 
- input: vcf.gz
- tool: SnpEff
- reference: GRCh38

04_extract_filter_snp.sh
- description: Convert the vcf file into a tabular format and filter the results based from BIOTYPE,IMPACT,DP and QUAL
- input: vcf
- tool: SnpSift (extractFields), awk, vcfEffOnePerLine.pl

05_extract_field_gnomad.sh
- description: Prepare and create a Gnomad reference for the filetered variants
- input: vcf from gnomad
- tool: gatk IndexFeatureFile, VariantsToTable

06_gnomad_lookup.sh
- description: Check for shared and unshared variants with the given CHROM-POS-REF-ALT of the variants in Gnomad and the variants to look up
- input: gnomad.table (from step 05), vcf (from step 04)
- tools: intersect.py (custom script)