#' Manually input pathway to reference database required for phyloseq analysis
classification.ref.ITS <- "D:/reference_data/All_Eukaryotes/UNITE_10.0_sh_general_release_dynamic_all_04.04.2024.fasta"
fusarioid.db <- "D:/reference_data/Fusaria/FUSARIOID-ID_v20251231/fusarioidIDversion31122025.fasta"

#' Load packages
cran_packages <- c("ggplot2", "stats", "scales", "gridExtra", "hexbin", "knitr", "lme4", 
                   "lmerTest", "pheatmap", "picante", "plyr", "beepr",
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

# Assign Variables for directories
path.0.project <- file.path(here()) # home directory, makes use of the package, "here"

# Set u beepr to notify you when Rstudio completes a step
# notify_when_idle <- function() {while(TRUE) {if (!any(grepl("R", ps -A$output()))) {beep() break} Sys.sleep(5)}}
beep("mario")
