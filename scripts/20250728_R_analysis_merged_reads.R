#' ---
#' title: "Evaluating the efficacy of FioreDono Primers for capturing Oomycetes using the Land Institute Root microbiome samples"
#' author: "Laura Los"
#' date: 2025-08-01
#' output: html_notebook
#' ---

#' Evaluating the efficacy of FioreDono Primers for capturing Oomycetes using the Land Institute Root microbiome samples (oomycete metabarcoding) 20250725 LL
#' Sequencing data generated using Illumina iSeq machine. (20250725 LL)
#' this script was run on the computer BU213-9YC1Z24 (20250725 LL)

# Load packages
cran_packages <- c("ggplot2", "gridExtra", "hexbin", "knitr", "lme4", 
                   "lmerTest", "pheatmap", "picante", "plyr", 
                   "RColorBrewer", "reshape2", "tidyr", "vegan", "here")
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

# # Use cutadapt to keep only reads with primer match, and to perform other filtering (N's, length criteria)
# # see ~/R_scripts/20240717_TLI_iSeq_cutadapt.R
# # in Ubuntu terminal: "R CMD BATCH ~/R_scripts/20240717_TLI_iSeq_cutadapt.R"
# # Cutadapt version 3.5 with Python 3.10.12 (updated 20240717 LL)
# 
# ## Running Cutadapt
# 
# ## in Ubuntu terminal (in the folder that contains the script): 
# # R CMD BATCH Microbiomes_SequenceProcessing_cutadapt.R
# 
# set.seed(1024)
# 
# ### Load data
# path.0.project <- file.path("/mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/")
# path.2.data <- file.path("/mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/2.data")
# setwd(path.0.project)
# list.files(path.0.project)
# 
# # Sort ensures forward/reverse reads are in same order
# names.2.data.R1 <- sort(list.files(path.2.data, pattern = "_R1_001.fastq.gz"))
# names.2.data.R2 <- sort(list.files(path.2.data, pattern = "_R2_001.fastq.gz"))
# 
# # expand names to include full path
# names.2.data.R1 <- file.path(path.2.data, names.2.data.R1)
# names.2.data.R2 <- file.path(path.2.data, names.2.data.R2)
# 
# # Identify primer sequences #done
# f <- "GCGGAAGGATCATTACCAC" #20250725_LL
# f.rc <- "GTGGTAATGATCCTTCCGC" #20250725_LL
# r <- "TCTTCATCGDTGTGCGAGC" #20250725_LL
# r.rc <- "GCTCGCACAHCGATGAAGA" #20250725_LL
# flags.R1 <- paste("-g", f, "-a", r.rc) 
# flags.R2 <- paste("-G", r, "-A", f.rc) 
# 
# # Run Cutadapt
# cutadapt <- "/usr/bin/cutadapt"
# system2(cutadapt, args = "--version")
# # 3.5
# # Cutadapt flags: 
# # "-g" (Read 1) / "-G" (Read 2) = regular 5' adapter (includes trimming before adapter)
# # "-a" (Read 1) / "-A" (Read 2) = regular 3' adapter (includes trimming after adapter)
# 
# path.3.cut <- file.path(path.0.project, "3.cut")
# if(!dir.exists(path.3.cut)) dir.create(path.3.cut)
# names.3.cut.R1 <- file.path(path.3.cut, basename(names.2.data.R1))
# names.3.cut.R2 <- file.path(path.3.cut, basename(names.2.data.R2))
# 
# names.3.cut.R1[1] #[1] means first element in this variable
# names.3.cut.R2[1]
# names.2.data.R1[1]
# names.2.data.R2[1]
# 
# for(i in seq_along(names.2.data.R1)) {
#    system2(cutadapt, args = c(flags.R1, flags.R2, 
#                              "-n", 2, #Remove up to "-n" adapters from each read
#                              "--discard-untrimmed", #Discard reads that do not contain an adapter
#                              "--pair-filter=any", #Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
#                              "--minimum-length", 50, #Discard reads shorter than 50, quality filtering
#                              "--max-n", 0, #I dont understand what 0 means here
#                              "-o", shQuote(names.3.cut.R1[i]), 
#                              "-p", shQuote(names.3.cut.R2[i]), # output file names
#                              shQuote(names.2.data.R1[i]), 
#                              shQuote(names.2.data.R2[i]))) # input file names
# }


# Count number of reads in which the primer is found BEFORE cutadapt applied (this step takes some time)
primerHits <- function(primer, names) {
  nhits <- vcountPattern(primer, sread(readFastq(names)), fixed = FALSE, max.mismatch = 1)
  return(sum(nhits > 0))
}

r1s <- sapply(c(f, r.rc), primerHits, names = names.2.data.R1)
r1s
# GCGGAAGGATCATTACCAC GCTCGCACAHCGATGAAGA (20250725 LL)
#             2125006               18704 

r2s <- sapply(c(r, f.rc), primerHits, names = names.2.data.R2)
r2s
# TCTTCATCGDTGTGCGAGC GTGGTAATGATCCTTCCGC (20250725 LL)
#             3661001               16290 

# Count number of reads in which the primer is found AFTER cutadapt applied (this step takes some time)
path.3.cut <- file.path(path.0.project, "3.cut")
if(!dir.exists(path.3.cut)) dir.create(path.3.cut)
names.3.cut.R1 <- file.path(path.3.cut, basename(names.2.data.R1))
names.3.cut.R2 <- file.path(path.3.cut, basename(names.2.data.R2))

sapply(c(f, r.rc), primerHits, names = names.3.cut.R1)
# GCGGAAGGATCATTACCAC GCTCGCACAHCGATGAAGA 
#                   0                   0 
sapply(c(r, f.rc), primerHits, names = names.3.cut.R2)
# TCTTCATCGDTGTGCGAGC GTGGTAATGATCCTTCCGC 
#                   0                   0 

## Filter & trim on quality ####
path.4.filter <- file.path(path.0.project, "4.filter")
if(!dir.exists(path.4.filter)) dir.create(path.4.filter)

names.4.filter.R1 <- file.path(path.4.filter, basename(names.3.cut.R1))
names.4.filter.R2 <- file.path(path.4.filter, basename(names.3.cut.R2))
length(names.4.filter.R1)

# number of samples = 15 (20250725 LL)
library(stats)
library(ggplot2)
qual.plot_R1 <- plotQualityProfile(names.3.cut.R1[sample(seq(1:10), 6, replace = FALSE)])
# plotQualityProfile(names.3.cut.R1[1]) #this works 20250728 LL
# X-axis = base position (cycle)
# Y-axis = quality score (Q)
# Darker blue = median quality
# Orange line = 25th percentile (you care about this for filtering)

print(qual.plot_R1)

names.3.cut.R1
ggsave("plotQualityProfile_6random_R1_20250728.pdf", qual.plot_R1)

# 20240812 LL - Looks OK, odd drop around the 20th cycle consistently...
#20250728 LL - the number of reads per sample were inconcsistent (chalk this up to not repeating the qPCR to get it within the standard curve) otherwise the quality plots were very similar to the previous quality plots where there is a dip in quality around 20 bases into the read and then again a dip for the last 20 bases. Another difference is there are two samples that in the random subset that have multiple dips throughout the cycles... 

qual.plot_R2 <- plotQualityProfile(names.3.cut.R2[sample(seq(1:10), 6, replace = FALSE)]) #did work #picking a random subset of 6 samples, we only have 11. Matt simplified this sample(seq(1:11), 6) ## don't want to look at negative control or mock reference samples here
names.3.cut.R2
print(qual.plot_R2)
ggsave("plotQualityProfile_6random_R2_20250728.pdf", qual.plot_R2)

filter <- filterAndTrim(names.3.cut.R1, names.4.filter.R1, 
                        names.3.cut.R2, names.4.filter.R2, 
                        maxN = 0,
                        maxEE = c(2, 2),
                        truncLen = c(125, 120), #truncLen: Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded. #this may need to be changed to 121 since this is the max length for the reads?
                        #if truncLen = c(125, 120) = Keep only the first 125 bases from the forward read
                        # Keep only the first 120 bases from the reverse read
                        compress = TRUE) #EE = sum(10^(-Q/10)), therefore maxEE = 2 means that if you sum all the bases Q values the max value that will pass filtering will be a sum of 2 (Q10= +0.1, Q20= +0.01, Q30= +0.001)  
filter

# reads.in reads.out
# 3-xo_S2_L001_R1_001.fastq.gz           100313     97147
# 3_S1_L001_R1_001.fastq.gz              107912    104843
# 6-xo_S4_L001_R1_001.fastq.gz            97287     93578
# 6_S3_L001_R1_001.fastq.gz               95228     92738
# 8-xo_S6_L001_R1_001.fastq.gz           104905     95719
# 8_S5_L001_R1_001.fastq.gz              119098    116158
# Ae-DMSO_S12_L001_R1_001.fastq.gz       173349    168662
# Ae-xo_S11_L001_R1_001.fastq.gz         301900    294100
# MC-DMSO_S8_L001_R1_001.fastq.gz        174246    170318
# MC-xo_S7_L001_R1_001.fastq.gz          184702    180370
# NTC-DMSO_S14_L001_R1_001.fastq.gz        6478      6376
# NTC-xo_S13_L001_R1_001.fastq.gz        113596    111422
# Pi-DMSO_S10_L001_R1_001.fastq.gz       135316    131884
# Pi-xo_S9_L001_R1_001.fastq.gz          220547    212901
# Undetermined_S0_L001_R1_001.fastq.gz    93434     70669

#20250728 hmmm my negative control still has a lot of reads... thats no good :(
#20250728 therefore I will skipthe following step

# negative control sample reduced to 0 reads, no filtered file created because all reads removed (updated 20240718 LL)
# remove NTC from file name list
# names.4.filter.R1 <- names.4.filter.R1[c(1:9,11)]
# names.4.filter.R2 <- names.4.filter.R2[c(1:9,11)]
# names.4.filter.R1
# names.4.filter.R2

sum(filter[, 1])
# 1,695,349 reads in (20240718 LL)
# 2,028,311 reads in (20250728)
sum(filter[, 2])
# 1,637,523 reads out (20240718 LL)
# 1,946,885 reads out (20250728)
sum(filter[, 1]) - sum(filter[, 2])
# 57,826 reads culled (20240718 LL)
# 81,426 reads culled (20250728)

print(names.4.filter.R1)
length(names.4.filter.R2)

filtered.qual.plot.R1 <- plotQualityProfile(names.4.filter.R1[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R1)
ggsave("plotQualityProfile_post_filter_6_random_plot_R1.pdf", filtered.qual.plot.R1, width= 7, height = 7)
# looks a lot better 20250728

filtered.qual.plot.R2 <- plotQualityProfile(names.4.filter.R2[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R2)
ggsave("plotQualityProfile_post_filter_6_random_plot_R2.pdf", filtered.qual.plot.R2, width= 7, height = 7)
# 20250728 R2 reads seem to be lower quality, still seeing a dip around 20 cycles/bases in

# This is where I start getting different numbers than Matt's analysis ####
# DADA2 Analysis ####
## learnErrors - Learn and Plot Errors (uses HMM) ####
err.R1 <- learnErrors(names.4.filter.R1, multithread = TRUE, randomize = TRUE) #based on dada2 workflow, we are "learning" the relationship between bases and error rate (looks for unique sequences and if it is identical to a sequence that it coming up 500 times except for 1 base, that 1 base is determined to be an "error")
#103,339,375 total bases in 826,715 reads from 5 samples will be used for learning the error rates. (20240718 LL)
#100,895,875 total bases in 807,167 reads from 6 samples will be used for learning the error rates.20250728
plotErrors(err.R1)
# Warning message: In scale_y_log10() : log-10 transformation introduced infinite values. (20240718 LL)
# 20250728:
# Warning messages:
#   1: In scale_y_log10() :
#   log-10 transformation introduced infinite values.
# 2: Removed 82 rows containing missing values or values outside the scale range
# (`geom_line()`).

err.R2 <- learnErrors(names.4.filter.R2, multithread = TRUE, randomize = TRUE) #should we be using set.seed for this? #askmatt
#114,073,920 total bases in 950,616 reads from 6 samples will be used for learning the error rates. (20240718 LL)
#114,314,880 total bases in 952,624 reads from 6 samples will be used for learning the error rates.
plotErrors(err.R2)
# plot not produced: Warning message: In scale_y_log10() : log-10 transformation introduced infinite values. (20240718 LL) ####

# 20250728:
# Warning messages:
#   1: In scale_y_log10() :
#   log-10 transformation introduced infinite values.
# 2: Removed 82 rows containing missing values or values outside the scale range
# (`geom_line()`).

## derep - dereplicating amplicon sequences ####
derep.R1 <- derepFastq(names.4.filter.R1, verbose = TRUE)
# 20250728 Encountered 4401 unique sequences from 70669 total sequences read.
names(derep.R1) <- substr(names(derep.R1), 1, 11)
derep.R2 <- derepFastq(names.4.filter.R2, verbose = TRUE)
# 20250728 Encountered 5257 unique sequences from 70669 total sequences read
names(derep.R2) <- substr(names(derep.R2), 1, 11)

### dada - Infer sequence variants ####
dada.R1 <- dada(derep.R1, err = err.R1, multithread = TRUE)  #changes data based on likely errors, 
dada.R1[[5]] #look at sample 5 to see the result
# 20240718 496 sequence variants were inferred from 20926 input unique sequences. Same number of input seq as Matt (20240718 LL)
# 20250728 72 sequence variants were inferred from 2741 input unique sequences.

dada.R2 <- dada(derep.R2, err = err.R2, multithread = TRUE)
dada.R2[[5]]
# 20240718 450 sequence variants were inferred from 21488 input unique sequences. Same number of input seq as Matt (20240718 LL)
# 20250728 61 sequence variants were inferred from 3358 input unique sequences.

# ideally, would track whether / how many reads are culled here. culled = discarded (animals not chosen for breeding)

# mergePairs - merge R1 & R2 --SKIPPED 20240813 LL ####
merge <- mergePairs(dada.R1, derep.R1, #this may fail because iseq length. nmatch is present therefore looks like it works
                   dada.R2, derep.R2,
                   maxMismatch = 2, verbose = TRUE) # in this protocol we wait to merge R1 and R2, keeping quality info for much longer, and only merging Read1 and read2 after you have corrected for sequence errors.
head(merge[[1]], 15) #print the first 15 lines of merge
head(merge[[1]], 1)

# might be losing a lot of reads here.
# remove large objects that are not needed anymore - not working with a lot of data for this so I added # infront of the following four lines
# rm(dada.R1)
# rm(dada.R2)
# rm(derep.R1)
# rm(derep.R2)

## Make sequence tables from dada.R1 (ie just forward reads - not merged)
seqtab <- makeSequenceTable(dada.R1) #makeSequenceTable is a dada2 function that generates a sequence table from dada2 object
#20240718 dada.R1 used instead of merge because sequences did not merge well
dim(seqtab) #dimensions of an object, if a table it will return the vector c(nrow, ncol), ie number of rows and number of columns.
# 20240718 10 samples X 2,489 sequence variants (20240718 LL)
# 20250728 15 samples x 386 sequence variants
sum(seqtab)
# 20240718   797,006 reads (20240718 LL)
# 20250728 1,946,112 reads 
unname(seqtab[, 1:10]) # removes the column names from the first 10 columns of the sequence table
seq.len.distr.tab <- table(nchar(getSequences(seqtab))) #provides a frequency table of the lengths of the sequences in the sequence table
# getSequences() is a dada2 function, that extracts the sequences (the unique sequence variants, or ASVs),
# nchar() function calculates the number of characters in a string, 
# table() function creates a frequency table (ie. how many sequences have each possible length).
barplot(seq.len.distr.tab, main = "Frequency of Sequence Lengths", xlab = "Sequence Length", ylab = "Frequency", col = "blue")
pdf(file = "Freq_Seq_Length_from_seqtab_R1.pdf")
# 20250728 all the sequences are 125 bases long - I think this is because I didnt merge the reads because there was no overlap

## Make sequence table from dada.R2 data
seqtab.R2 <- makeSequenceTable(dada.R2)
dim(seqtab.R2)
sum(seqtab.R2)

## Make sequence tables from merged reads
seqtab <- makeSequenceTable(merge) #makeSequenceTable is a dada2 function that generates a sequence table from dada2 object
dim(seqtab) #dimensions of an object, if a table it will return the vector c(nrow, ncol), ie number of rows and number of columns.
# 20240718 10 samples X 2,489 sequence variants (20240718 LL)
# 20250728 15 samples x 386 sequence variants from dada.R1
# 20250728 15 samples x 701 sequence variants from merge
sum(seqtab)
# 20240718   797,006 reads (20240718 LL)
# 20250728 1,946,112 reads from dada.R1
# 20250728 1,945,689 reads from dada.R2
# 20250728 1,436,298 reads from merge
library(dada2)


unname(seqtab[, 1:10]) # removes the column names from the first 10 columns of the sequence table
seq.len.distr.tab <- table(nchar(getSequences(seqtab))) #provides a frequency table of the lengths of the sequences in the sequence table
                                                        # getSequences() is a dada2 function, that extracts the sequences (the unique sequence variants, or ASVs),
                                                        # nchar() function calculates the number of characters in a string, 
                                                        # table() function creates a frequency table (ie. how many sequences have each possible length).
barplot(seq.len.distr.tab, main = "Frequency of Sequence Lengths", xlab = "Sequence Length", ylab = "Frequency", col = "blue")
pdf(file = "Freq_Seq_Length_from_seqtab_merge.pdf")
# 20250728 from merge we get a variety of lengths and frequencies 

# cull unexpectedly short or long reads
# skipped because I do not know whether to expect length heterogeneity among oomycetes
## seqtab.length <- seqtab[, nchar(colnames(seqtab)) %in% seq(251, 257)]
## dim(seqtab.length)
# 72 samples X 18,892 sequence variants
## sum(seqtab.length)
# 7,812,529 reads

#seqtab.IgnoreEnds <- collapseNoMismatch(seqtab.length, minOverlap = 200)
#dim(seqtab.IgnoreEnds)
# this got hung up and I skipped it (MB)

## Remove chimeras
seqtab.chimera <- removeBimeraDenovo(seqtab,
                                     multithread = TRUE, verbose = TRUE, method = "consensus")
# 20240808 Identified 1,621 bimeras out of 2,489 input sequences. (20240808 LL) comapred to 1,666 bimeras out of 2,537 input sequences identified by MB
# 20250728 Identified 515 bimeras out of 701 input sequences.

dim(seqtab.chimera)
# 20240808 10 samples X 868 sequence variants (20240808 LL, rows are samples, columns are unique sequences)
# 20250728 15 samples x 186 sequence variants

sum(seqtab.chimera)
# 20240808   716,317 reads (20240808 LL)
# 20250728 1,384,835 reads
1 - (sum(seqtab.chimera) / sum(seqtab))
# 20240808 0.1012401 proportion of chimeric reads #AskMAtt is this high? this seems kinda high (10%)
# 20250728 0.03583031 proportion of chimeric reads - this doesnt seem too bad to me LL


## 2025-08-01 I separated out the scripts up until this point #####################

## Trim off conserved flanking regions
path.5.ITSx <- file.path(path.0.project, "5.ITSx")
if(!dir.exists(path.5.ITSx)) dir.create(path.5.ITSx)
ITSx.input <- DNAStringSet(colnames(seqtab.chimera))
ASV.names <- seq(dim(seqtab.chimera)[2])
names(ITSx.input) <- ASV.names #assigns ASV,names (which is a sequence of numbers) to the columns of ITSx.input
writeXStringSet(ITSx.input, file.path(path.5.ITSx, "ITSx.input.fasta"), format = "fasta")

# in Ubuntu:
# ITSx --license
# v. 1.1.3 (20240808 LL)
# v. 1.1.3 (20250728 LL)

# permit ITSx to flag reads as belonging to other than fungi:
# ITSx -i /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/5.ITSx/ITSx.input.fasta -o /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/5.ITSx/ITSx_out --allow_reorder T --allow_single_domain 1e-3,0 --partial 50
# 20240809 output suggests mostly fungi (589/864 = 68.17%) and plants, oomycetes (14/864 = 1.62%) (20240809 LL)
# 20250728 looking at the file: ITSx_out.summary.txt:
#   20250728 output suggests mostly oomycetes (63/186 = 33.87%) and some limited fungi (4/186 = 2.15%) 
#   20250728 therefore this suggests that FioreDonno primers may do a better job extracting just the oomycetes, but why od I have so few sequences?

# # pull out just the oomycete reads:
# cd /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/5.ITSx/
# grep -A 1 "|O" ITSx_out.ITS1.full_and_partial.fasta > ITSx.oomycete.fasta
# sed -i 's/|/ /' ITSx.oomycete.fasta
# sed -i 's/\s.*$//' ITSx.oomycete.fasta
# sed -i 's/--//' ITSx.oomycete.fasta

# 20250728 note: ITS2 is empty but thsi makes sense since the Fiore Donno Primers capture only ITS1

oomycetes <- readDNAStringSet(file.path(path.5.ITSx, "ITSx.oomycete.fasta"))
width(oomycetes)
length(oomycetes)
# 20250728 63
length(ITSx.input)
length(ITSx.input) - length(oomycetes)
# 854 of the original 868 sequence variants were not recognizable as oomycetes ITS reads, 14 categorized as oomycetes (20240812 LL)
# 20250728 123 of the original 186 sequence variants were not recognized as oomycete ITS reads, 63 categorized as oomycetes

## Assign taxonomy ####
classification.ref.ITS <- "D:/reference_data/All_Eukaryotes/UNITE_10.0_sh_general_release_dynamic_all_04.04.2024.fasta"
set.seed(1024)
taxtab.ITS.oomycetes <- assignTaxonomy(oomycetes, refFasta = classification.ref.ITS, multithread = TRUE, minBoot = 80)
dim(taxtab.ITS.oomycetes)
# 20240812 14 rows (14 sequence variants, with 3 having unassigned kingdoms etc.), 7 columns (taxonomic ranks) for the taxtab.ITS database (20240812 LL)
# 20250728 63 rows (63 sequence variants), 7 columns (taxonomic ranks) for the taxtab.ITS database
taxtab.ITS.oomycetes <- data.frame(taxtab.ITS.oomycetes)
head(taxtab.ITS.oomycetes)

oomycetes.df <- data.frame(oomycetes)
# row names tell which columns of seqtab.chimera (made earlier on in analysis, when removing chimeras) should be retained, and which should be dropped as non-oomycetes
dim(seqtab.chimera)
# 20240812 10 X 868, so 868 sequence variants (20240812 LL)
# 20250728 15 x 186, so 186 sequence variants
seqtab.oomycetes <- seqtab.chimera[, as.numeric(rownames(oomycetes.df))]
dim(seqtab.oomycetes)
# 20240812 10 X 14, extract 14 variant sequences from chimera-removed sequence table as oomycetes (20240812 LL)
# 20250728 15 x 63, extract 63 variant sequences from chimera-removed sequence table as oomycetes

# problem using actual sequences as ASV identifiers, because post-ITSx there are duplicate sequences...
sum(duplicated(oomycetes.df$oomycetes)) #gives the total count of duplicated sequences. (20240812 LL)
# 20240812 1 #AskMatt what does this mean?
# 20250728 1 #this is good right, no duplicated sequences? 

# 20240812 Keep using the pre-ITSx sequences instead #ASkMatt why?
# 20250728 I think its becasue ITSx is trimming the reads so they are only within the ITSx region,
#     we want to retain as much info as possible but ITSx was useful for identifying our taxa

rownames(taxtab.ITS.oomycetes) <- colnames(seqtab.oomycetes)

# Make phyloseq objects
# metadata <- metadata[-10, ] 
# 20240812 removing the negative control I had in my data (which was row 10) (20240812 LL)
# 20250728 I commented this out because I am worried about contamination in my negative control so I'm leaving it in 

# # in order to use phyloseq we need all the sample names to match for the otu table, the tax table and the sample data
# colnames(seqtab.oomycetes)         # from otu_table
# rownames(taxtab.ITS.oomycetes)     # from tax_table
# rownames(metadata)                 # from sample_data
# # if the colnames or row names are sequences intead of sample names you will need to transpose the dataframes accordingly:
# seqtab.oomycetes <- t(seqtab.oomycetes)
# colnames(seqtab.oomycetes)
# taxtab.ITS.oomycetes <- t(taxtab.ITS.oomycetes)
# rownames(taxtab.ITS.oomycetes)

# assigned to metadata rownames based off of: colnames(seqtab.oomycetes)
colnames(seqtab.oomycetes)

rownames(metadata) <- c("3-xo_S2_L00", "3_S1_L001_R", "6-xo_S4_L00",
                        "6_S3_L001_R", "8-xo_S6_L00", "8_S5_L001_R",
                        "Ae-DMSO_S12", "Ae-xo_S11_L", "MC-DMSO_S8_",
                        "MC-xo_S7_L0", "NTC-DMSO_S1", "NTC-xo_S13_",
                        "Pi-DMSO_S10", "Pi-xo_S9_L0")
oomycetes.ps <- phyloseq(otu_table(seqtab.oomycetes, taxa_are_rows = FALSE),
                         tax_table(as.matrix(taxtab.ITS.oomycetes)), sample_data(metadata))

oomycetes.ps
# read the OTU table
# Extract abundance matrix from the phyloseq object
OTU.oomycetes = as(otu_table(oomycetes.ps), "matrix")
## open the OTU table
dim(OTU.oomycetes)
# 20240812 10 samples X 14 taxa (20240812 LL)
# 20250728 14 samples x 63 taxa 
rownames(OTU.oomycetes)

path.6.phyloseq <- file.path(path.0.project, "6.phyloseq")
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)
save(oomycetes.ps, file = file.path(path.6.phyloseq, "oomycetes.ps.RData"))

# save R environment to lab workstation.

#install and load writexl package to save taxtabs
library(writexl)
write_xlsx(taxtab.ITS.oomycetes, file.path(path.0.project, "taxtab.ITS.oomycetes.xlsx"))

write.table(tax_table(oomycetes.ps), file = file.path(path.6.phyloseq, "oomycetes.ps.taxonomy.txt"))
write.table(otu_table(oomycetes.ps), file = file.path(path.6.phyloseq, "oomycetes.ps.abundances.txt"))

# 20250728 - Laura from 20240812 performed more analysis because I was having a problem... 
#   where Several of the ASVs called oomycetes by ITSx are classified as fungi (Acrocalymma sp.) (look at taxtab.ITS) ####
#   I did not have this problem so I concluded the analysis at this stage.
#   just kidding I wanted to get a better breakdown of the species in each sample 

# create a table containing taxonomy and abundance: ####
oomycete.taxtab <- as.data.frame(tax_table(oomycetes.ps)) # this has your sequences as rownames and the taxonomic assignments as colnames
colnames(oomycete.taxtab)
rownames(oomycete.taxtab)

oomycetes.otutab <- as.data.frame(otu_table(oomycetes.ps))
colnames(oomycetes.otutab) #sequences
rownames(oomycetes.otutab) #sample_IDs

#20250728 therefore I will transpose otutab so that sequences are rownames
oomycetes.otutab.t <- t(oomycetes.otutab)
colnames(oomycetes.otutab.t) #sample_IDs
rownames(oomycetes.otutab.t) #sequences

#merge otutab and taxtab by rownames
oomycete.taxtab$ASV <- rownames(oomycete.taxtab)

oomycetes.otutab.t <- as.data.frame(oomycetes.otutab.t)
oomycetes.otutab.t$ASV <- rownames(oomycetes.otutab.t)

species.abundance <- merge(oomycetes.otutab.t, oomycete.taxtab, by = "ASV")

head(species.abundance)

write_xlsx(species.abundance, file.path(path.0.project, "species.abundance.xlsx"))

head(species.abundance)

# create figure of oomycete genera counts ####

library(tidyr)   # for pivot_longer()
library(dplyr)   # for %>%, mutate()
library(ggplot2) # for plotting

taxa.long <- species.abundance %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Genus, sample, count)) %>% 
  replace_na(replace = list(Genus =  "other"))
  

print(taxa.long, n = 50)

genus.counts <- taxa.long %>% 
  group_by(sample, Genus) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(genus.counts, n = 50)
str(genus.counts)
length(unique(genus.counts$sample))
#14
length(unique(genus.counts$Genus))
#7
7*14
#98

summarise

# ➋ Plot
genus.plot <- ggplot(genus.counts, aes(x = sample, y = count, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
genus.plot
ggsave(plot = genus.plot, "genus_counts.png")

# create figure of oomycete species counts ####

library(tidyr)   # for pivot_longer()
library(dplyr)   # for %>%, mutate()
library(ggplot2) # for plotting

species.long <- species.abundance %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Genus, Species, sample, count)) %>%   
  mutate(taxa = coalesce(Species, Genus)) %>% # if Species is NA, then it takes data from Genus column 
  replace_na(replace = list(taxa = "other"))

print(species.long, n = 100)

species.counts <- species.long %>% 
  group_by(sample, taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(species.counts, n = 50)
str(species.counts)
length(unique(species.counts$sample))
#14
length(unique(species.counts$taxa))
#15
15*14
#210
str(species.counts)
#210 rows

# ➋ Plot
taxa.plot <- ggplot(species.counts, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
taxa.plot
ggsave(plot = taxa.plot, "genus_counts.png")

# investigate Phytophthora infestans taxonomic assignments #### 

# # pull out just the oomycete reads:
# cd /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/5.ITSx/
# cp ITSx_out.ITS1.full_and_partial.fasta ITSx.ITS1.all.fasta
# sed -i 's/|/ /' ITSx.ITS1.all.fasta
# sed -i 's/\s.*$//' ITSx.ITS1.all.fasta
# sed -i 's/--//' ITSx.ITS1.all.fasta

all.ITSx <- readDNAStringSet(file.path(path.5.ITSx, "ITSx.ITS1.all.fasta"))
width(all.ITSx)
length(all.ITSx)
# 20250730 109
length(ITSx.input)
# 188
length(ITSx.input) - length(all.ITSx)
# 79
109-63
# 20250730 79 of the original 188 sequence variants were not recognized as ITS reads, 63 categorized as oomycetes, 46 non-oomyctes remaining

## Assign taxonomy ####
classification.ref.ITS <- "D:/reference_data/All_Eukaryotes/UNITE_10.0_sh_general_release_dynamic_all_04.04.2024.fasta"
set.seed(1024)
taxtab.ITS.all.ITSx <- assignTaxonomy(all.ITSx, refFasta = classification.ref.ITS, multithread = TRUE, minBoot = 80)
dim(taxtab.ITS.all.ITSx)
# 20240812 14 rows (14 sequence variants, with 3 having unassigned kingdoms etc.), 7 columns (taxonomic ranks) for the taxtab.ITS database (20240812 LL)
# 20250728 63 rows (63 sequence variants), 7 columns (taxonomic ranks) for the taxtab.ITS database
taxtab.ITS.all.ITSx <- data.frame(taxtab.ITS.all.ITSx)
head(taxtab.ITS.all.ITSx)

all.ITSx.df <- data.frame(all.ITSx)
# row names tell which columns of seqtab.chimera (made earlier on in analysis, when removing chimeras) should be retained, and which should be dropped as non-all.ITSx
dim(seqtab.chimera)
# 20240812 10 X 868, so 868 sequence variants (20240812 LL)
# 20250728 15 x 186, so 186 sequence variants
seqtab.all.ITSx <- seqtab.chimera[, as.numeric(rownames(all.ITSx.df))]
dim(seqtab.all.ITSx)
# 20240812 10 X 14, extract 14 variant sequences from chimera-removed sequence table as all.ITSx (20240812 LL)
# 20250728 15 x 63, extract 63 variant sequences from chimera-removed sequence table as all.ITSx

# problem using actual sequences as ASV identifiers, because post-ITSx there are duplicate sequences...
sum(duplicated(all.ITSx.df$all.ITSx)) #gives the total count of duplicated sequences. (20240812 LL)
# 20240812 1 #AskMatt what does this mean?
# 20250728 1 #this is good right, no duplicated sequences? 

# 20240812 Keep using the pre-ITSx sequences instead #ASkMatt why?
# 20250728 I think its becasue ITSx is trimming the reads so they are only within the ITSx region,
#     we want to retain as much info as possible but ITSx was useful for identifying our taxa

rownames(taxtab.ITS.all.ITSx) <- colnames(seqtab.all.ITSx)

# Make phyloseq objects
# metadata <- metadata[-10, ] 
# 20240812 removing the negative control I had in my data (which was row 10) (20240812 LL)
# 20250728 I commented this out because I am worried about contamination in my negative control so I'm leaving it in 

# # in order to use phyloseq we need all the sample names to match for the otu table, the tax table and the sample data
# colnames(seqtab.all.ITSx)         # from otu_table
# rownames(taxtab.ITS.all.ITSx)     # from tax_table
# rownames(metadata)                 # from sample_data
# # if the colnames or row names are sequences intead of sample names you will need to transpose the dataframes accordingly:
# seqtab.all.ITSx <- t(seqtab.all.ITSx)
# colnames(seqtab.all.ITSx)
# taxtab.ITS.all.ITSx <- t(taxtab.ITS.all.ITSx)
# rownames(taxtab.ITS.all.ITSx)

# assigned to metadata rownames based off of: colnames(seqtab.all.ITSx)
colnames(seqtab.all.ITSx)

rownames(metadata) <- c("3-xo_S2_L00", "3_S1_L001_R", "6-xo_S4_L00",
                        "6_S3_L001_R", "8-xo_S6_L00", "8_S5_L001_R",
                        "Ae-DMSO_S12", "Ae-xo_S11_L", "MC-DMSO_S8_",
                        "MC-xo_S7_L0", "NTC-DMSO_S1", "NTC-xo_S13_",
                        "Pi-DMSO_S10", "Pi-xo_S9_L0")
all.ITSx.ps <- phyloseq(otu_table(seqtab.all.ITSx, taxa_are_rows = FALSE),
                         tax_table(as.matrix(taxtab.ITS.all.ITSx)), sample_data(metadata))

all.ITSx.ps
# read the OTU table
# Extract abundance matrix from the phyloseq object
OTU.all.ITSx = as(otu_table(all.ITSx.ps), "matrix")
## open the OTU table
dim(OTU.all.ITSx)
# 20240812 10 samples X 14 taxa (20240812 LL)
# 20250728 14 samples x 63 taxa 
rownames(OTU.all.ITSx)

path.6.phyloseq <- file.path(path.0.project, "6.phyloseq")
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)
save(all.ITSx.ps, file = file.path(path.6.phyloseq, "all.ITSx.ps.RData"))

# save R environment to lab workstation.

#install and load writexl package to save taxtabs
library(writexl)
write_xlsx(taxtab.ITS.all.ITSx, file.path(path.0.project, "taxtab.ITS.all.ITSx.xlsx"))

write.table(tax_table(all.ITSx.ps), file = file.path(path.6.phyloseq, "all.ITSx.ps.taxonomy.txt"))
write.table(otu_table(all.ITSx.ps), file = file.path(path.6.phyloseq, "all.ITSx.ps.abundances.txt"))

# 20250728 - Laura from 20240812 performed more analysis because I was having a problem... 
#   where Several of the ASVs called all.ITSx by ITSx are classified as fungi (Acrocalymma sp.) (look at taxtab.ITS) ####
#   I did not have this problem so I concluded the analysis at this stage.
#   just kidding I wanted to get a better breakdown of the species in each sample 

# create a table containing taxonomy and abundance: ####
all.ITSx.taxtab <- as.data.frame(tax_table(all.ITSx.ps)) # this has your sequences as rownames and the taxonomic assignments as colnames
colnames(all.ITSx.taxtab)
rownames(all.ITSx.taxtab)

all.ITSx.otutab <- as.data.frame(otu_table(all.ITSx.ps))
colnames(all.ITSx.otutab) #sequences
rownames(all.ITSx.otutab) #sample_IDs

#20250728 therefore I will transpose otutab so that sequences are rownames
all.ITSx.otutab.t <- t(all.ITSx.otutab)
colnames(all.ITSx.otutab.t) #sample_IDs
rownames(all.ITSx.otutab.t) #sequences

#merge otutab and taxtab by rownames
all.ITSx.taxtab$ASV <- rownames(all.ITSx.taxtab)

all.ITSx.otutab.t <- as.data.frame(all.ITSx.otutab.t)
all.ITSx.otutab.t$ASV <- rownames(all.ITSx.otutab.t)

species.abundance.all.ITS <- merge(all.ITSx.otutab.t, all.ITSx.taxtab, by = "ASV")

head(species.abundance.all.ITS)

write_xlsx(species.abundance.all.ITS, file.path(path.0.project, paste("species.abundance.allITS.xlsx", Sys.Date())))

head(species.abundance.all.ITS)

species.ITS.long <- species.abundance.all.ITS %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Genus, Species, sample, count)) %>%   
  mutate(taxa = coalesce(Species, Genus)) %>% # if Species is NA, then it takes data from Genus column 
  replace_na(replace = list(taxa = "other"))

print(species.ITS.long, n = 100)

species.ITS.counts <- species.ITS.long %>% 
  group_by(sample, taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(species.ITS.counts, n = 50)
str(species.ITS.counts)
length(unique(species.ITS.counts$sample))
#14
length(unique(species.ITS.counts$taxa))
#22
22*14
#308
str(species.ITS.counts)
#308 rows

# ➋ Plot
taxa.ITS.plot <- ggplot(species.ITS.counts, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
taxa.ITS.plot
filename = paste("taxa_ITS_counts.png", Sys.time())
ggsave(filename = filename, plot = taxa.ITS.plot)

