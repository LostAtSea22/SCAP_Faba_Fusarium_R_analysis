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

# Assign Variables
path.0.project <- file.path(here()) # home directory, makes use of the package, "here"
path.1.data <- file.path(path.0.project, "data", "1.data")

path.2.cut <- here("data", "2.cut")
if(!dir.exists(path.2.cut)) dir.create(path.2.cut)
names.2.cut.R1 <- file.path(path.2.cut, basename(names.1.data.R1))
names.2.cut.R2 <- file.path(path.2.cut, basename(names.1.data.R2))

path.3.filter <- here("data", "3.filter")
if(!dir.exists(path.3.filter)) dir.create(path.3.filter)
names.3.filter.R1 <- file.path(path.3.filter, basename(names.2.cut.R1))
names.3.filter.R2 <- file.path(path.3.filter, basename(names.2.cut.R2))

path.4.dada <- here("data", "4.dada")
if(!dir.exists(path.4.dada)) dir.create(path.4.dada)
