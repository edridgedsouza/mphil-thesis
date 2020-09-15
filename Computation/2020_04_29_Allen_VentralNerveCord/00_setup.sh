#!/bin/bash -x

# GEO dataset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141807
# For paper 	A single-cell transcriptomic atlas of the adult Drosophila ventral nerve cord
# (Allen et al, 2020)

# First manually download both the supplementary files from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141807&format=file
mkdir -p data
tar -C ./data -xvf data/GSE141807_RAW.tar

# rm -f data/GSM*
