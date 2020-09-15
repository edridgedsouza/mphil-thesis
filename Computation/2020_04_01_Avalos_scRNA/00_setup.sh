#!/bin/bash -x

# GEO dataset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134722
# For paper Single cell transcriptome atlas of the Drosophila larval brain
# (Avalos et al, 2019)

# First manually download both the supplementary files
mkdir -p data
rm -rf ./data/GSM*

tar -C ./data -xvf data/GSE134722_RAW.tar
rename 's/FirstInstarLarvalBrain//' data/GSM*
rename 's/Condition//' data/GSM*
rename 's/10X_//' data/GSM*

for run in $(ls data/*_genes.tsv.gz)
do
    sample=$(basename $run _genes.tsv.gz)
    mkdir -p data/${sample}
    mv data/${sample}_genes.tsv.gz data/${sample}/features.tsv.gz
    mv data/${sample}_barcodes.tsv.gz data/${sample}/barcodes.tsv.gz
    mv data/${sample}_matrix.mtx.gz data/${sample}/matrix.mtx.gz
done



# To reset, do:
# rm -rf ./data/GSM*

# Sanity checks:
# diff <(zcat *Normal_finalaggr*barcodes* | tr '-' ' ' | awk '{print $1}'| sort | uniq) <(zcat *Normal_sample*barcodes* | tr '-' ' ' | awk '{print $1}' | sort | uniq)
# diff <(zcat *Starvation_finalaggr*barcodes* | tr '-' ' ' | awk '{print $1}'| sort | uniq) <(zcat *Starvation_sample*barcodes* | tr '-' ' ' | awk '{print $1}' | sort | uniq)
# Should come up blank
# Basically, superset of all barcodes in the individual samples will be equivalent to the barcodes in the aggregate file, IF you ignore the suffix after the dash.
# This is because the suffix tells you when two identical barcodes come from different libraries--collisions occur when you have multiple samples.

# So, you should be safe to use the finalaggr datasets.
