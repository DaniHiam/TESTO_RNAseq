# TESTO_RNAseq
12 wk resistance program in females ages 18-35

Source: https://github.com/DaniHiam/TESTO_RNAseq

# Introduction

Project Summary

Title / Aim: .

Methodology:
* Females conducted 12 weeks resistance training (1) PRE (n=22), and (2) POST (n=19) 
* At the end of the experiment, RNA was isolated from skeletal muscle both PRE and POST 12 weeks. 
* Total RNA’s were extracted using the Qiagen Universal AllPrep Kit (250), (Cat# 74106, Hilden, Germany).
* Nanodrop, Qubit and Tapestation performed before sending to Macrogen for analysis.
* We will be comparing genes between PRE to POST resistance training program.

## PRE PROCESSING
Fastqc (v0.11.9) was used to inspect sequence quality.
Multiqc (v1.13) was used to tabulate sequence quality, and mapping statistics
The human transcriptome was downloaded from ENSEMBL version 38. "Homo_sapiens.GRCh38.cdna.all.fa.gz" 
Kallisto (0.46.1) was used to map RNA-seq reads to the transcriptome.

## DATA ANALYSIS
Data were read into R v4.2.1 and duplicate lane data were aggregated, and transcript level counts were aggregated to gene level counts.
Genes with an average of less than 10 reads across samples were excluded from downstream analysis.
DESeq (1.32.0) was used with default settings to assay differential expression between control and treatment groups for all tissues [6].
Pathway analysis was performed with reactome gene sets obtained from MSigDB and converted to mouse gene identifiers with the msigdbr package 
(performed on 16-02-2022) [7,8,9].

Differential pathway analysis was performed with the "mitch" bioconductor package [10].

Genes and pathways with false discovery rate (FDR)<0.05 were considered significant.
