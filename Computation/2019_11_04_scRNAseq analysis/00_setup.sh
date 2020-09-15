#!/bin/bash -x

# GEO dataset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95025
# For paper The Drosophila embryo at single-cell transcriptome resolution.
# (Karaiskos et al, 2017)

# First manually download both the supplementary files from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134722&format=file
mkdir -p data
tar -C ./data -xvf data/GSE95025_RAW.tar
