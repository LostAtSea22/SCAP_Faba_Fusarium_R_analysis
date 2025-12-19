#' ---
#' title: "Evaluating the efficacy of FioreDono Primers for capturing Oomycetes using the Land Institute Root microbiome samples"
#' subtitle: "00.loading_packages_v1.r"
#' author: "Laura Los"
#' date_written: 2025-08-01
#' output: html_notebook
#' ---

#' Evaluating the efficacy of FioreDono Primers for capturing Oomycetes using the Land Institute Root microbiome samples (oomycete metabarcoding) 20250725 LL
#' 
#' Sequencing data generated using Illumina iSeq machine. (20250725 LL)
#' 
#' this script was run on the computer BU213-9YC1Z24 (20250725 LL)
#' 
#' saved as rmarkdown using `rmarkdown::render("00.loading_packages_v1.r", output_dir = here("scripts"))`
#' 
#' 

#' Manually input pathway to reference database required for phyloseq analysis
classification.ref.ITS <- "D:/reference_data/All_Eukaryotes/UNITE_10.0_sh_general_release_dynamic_all_04.04.2024.fasta"

#' Load packages
cran_packages <- c("ggplot2", "stats", "scales", "gridExtra", "hexbin", "knitr", "lme4", 
                   "lmerTest", "pheatmap", "picante", "plyr", 
                   "RColorBrewer", "reshape2", "dplyr", "tidyr", "vegan", "here", "writexl", "readxl")
inst <- cran_packages %in% installed.packages()
if(any(!inst)) {
  install.packages(cran_packages[!inst])
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bioc_packages <- c("apeglm", "BiocStyle", "Biostrings", "dada2", "DECIPHER", "decontam", 
                   "DESeq2", "phangorn", "phyloseq", "ShortRead", "vsn")
inst <- bioc_packages %in% installed.packages()
if(any(!inst)) {
  BiocManager::install(bioc_packages[!inst], ask = F)
}

sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

#Save session information
session_info <- capture.output(sessionInfo())
writeLines(session_info, here("session_info", paste(Sys.Date(), "session_info.txt", sep = "_")))