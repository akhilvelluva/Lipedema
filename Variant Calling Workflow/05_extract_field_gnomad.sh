#!/bin/bash
input="gnomad.genomes.r2.1.1.sites.vcf.bgz"
ref="GRCh38.p13"
output="gnomad"
mkdir -p $output

echo "__START_COMMAND__"
gatk IndexFeatureFile -I $input
gatk VariantsToTable -V $input -F CHROM -F POS -F ID -F REF -F ALT -F AC -F AN -F AF -F nhomalt -O Gnomad.table --show-filtered
echo "__END_COMMAND__"