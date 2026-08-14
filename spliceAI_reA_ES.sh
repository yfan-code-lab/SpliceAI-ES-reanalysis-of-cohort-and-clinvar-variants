#!/usr/bin/env bash
###############################################################################
# SpliceAI re-analysis of exome sequencing (ES) variants
#
# Per-sample pipeline (suitable for cohort VCFs from P170 onwards):
#   1. Extract OMIM-annotated variants from the gene-specific-target (gst) VCF
#   2. Filter to PASS, MSigDb != ".", and rare alleles (G.AF < 0.01 or missing)
#   3. Annotate with SpliceAI via vcfanno
#   4. Extract SpliceAI scores into a per-sample TSV
#   5. Keep variants with any SpliceAI delta score > 0.1 (*_spliceai_sig.tsv)
#   6. Merge all per-sample files into combined.tsv, then filter to delta > 0.5
###############################################################################

set -euo pipefail

# ── configuration ────────────────────────────────────────────────────────────
# Directory containing the input VCFs, the vcfanno config and the header file.
WORK_DIR=$(pwd)
# Where per-sample result files are collected.
OUTPUT_DIR="$WORK_DIR/output"
# Reference sample used to build the header of combined.tsv
# (set this to any one of your sample IDs).
REF_SAMPLE_ID="SAMPLE_ID"

mkdir -p "$OUTPUT_DIR"

# Requires the `vcfanno` conda environment (conda, vcfanno, bcftools, bgzip, tabix)
conda activate vcfanno

# Sample list file (first argument): one line per sample, e.g. "<ID> <path/to.gst.vcf>"
while read id gstvcf
do
        tab=$(echo $gstvcf | sed 's/gst.vcf/tab/')
        cat $tab | awk '$8!="."' > ${id}_omim.tab
        cp ${id}_omim.tab "$OUTPUT_DIR"
        bgzip -c $gstvcf > ${id}.gst.vcf.gz
        tabix ${id}.gst.vcf.gz
        bcftools view -i 'INFO/MSigDb!="." & FILTER="PASS" & (INFO/G.AF<0.01 | INFO/G.AF=".")' ${id}.gst.vcf.gz -Oz -o ${id}_omim.vcf.gz
        tabix ${id}_omim.vcf.gz
        vcfanno -p 4 vcfanno_spliceai_conf.toml ${id}_omim.vcf.gz > ${id}_omim_spliceai.vcf
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO/DP\t%INFO/MSigDb\t%G.controls_AF\t%G.AF\t%G.AF_male\t%IDT_new.AC\t%IDT.AC\t%INH.AC\t%SpliceAI\n' ${id}_omim_spliceai.vcf > ${id}_spliceai.tsv
        cat header > ${id}_spliceai_sig.tsv
        cat ${id}_spliceai.tsv | sed 's/|/\t/g' | awk -F'\t' '$17 > 0.1 || $18 > 0.1 || $19 > 0.1 || $20 > 0.1' >> ${id}_spliceai_sig.tsv
        cp ${id}_spliceai_sig.tsv "$OUTPUT_DIR"
done < "$1"

# ── merge all per-sample files ───────────────────────────────────────────────
# Create output file with header
echo -n "ID" | cat - <(head -n1 ${REF_SAMPLE_ID}_spliceai_sig.tsv) > combined.tsv

# Process each file
for f in *_sig.tsv; do
    # Extract ID from filename (remove '_spliceai_sig.tsv')
    id=${f/_spliceai_sig.tsv/}

    # Add ID column and append data lines (skip header)
    awk -v ID="$id" 'NR>1 {print ID "\t" $0}' "$f" >> combined.tsv
done
## formatting issue resolved in tsv file.

cat combined.tsv | awk '$21 > 0.5 || $18 > 0.5 || $19 > 0.5|| $20 > 0.5' >> combined_spliceai05.tsv
