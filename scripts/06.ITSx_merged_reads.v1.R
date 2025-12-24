#' date: 2025-08-05
#' author: "Laura Los"
#' title: "06.ITSx triming conserved flanking regions on merged reads"
#' version_notes: "

#------------

## Trim off conserved flanking regions ####

# create folder
path.5.ITSx <- file.path(here("data", "5.ITSx_merged_reads"))
if(!dir.exists(path.5.ITSx)) dir.create(path.5.ITSx)

#import sequence tables files
seqtab.merged.chimera.rm <- read.csv(here("data", "4.dada", "2025-12-23SeqTable_chimeras_rm_merged_pairs.csv"))
rownames(seqtab.merged.chimera.rm) <- seqtab.merged.chimera.rm$X
seqtab.merged.chimera.rm <- seqtab.merged.chimera.rm[,-1]
rownames(seqtab.merged.chimera.rm)
colnames(seqtab.merged.chimera.rm[,1:10])

#prepare data for ITSx
ITSx.input <- DNAStringSet(colnames(seqtab.merged.chimera.rm))

ASV.names <- seq(dim(seqtab.merged.chimera.rm)[2])
names(ITSx.input) <- ASV.names #assigns ASV,names (which is a sequence of numbers) to the columns of ITSx.input
writeXStringSet(ITSx.input, file.path(path.5.ITSx, "ITSx.input.fasta"), format = "fasta")

# use system 2 to access Windows subsystem for linux (wsl) and check "ITSx" version (save output to "data/5.ITSx") 
ITSx <- "/home/losl/miniconda3/bin/ITSx"
cat(paste(system2("wsl", args = c(ITSx, "--license"), stdout = TRUE, stderr = TRUE), collapse = "\n"), "\n")

ITSx_version <- capture.output(cat(paste(system2("wsl", args = c(ITSx, "--license"), stdout = TRUE, stderr = TRUE), collapse = "\n"), "\n"))
writeLines(ITSx_version, file.path(path.5.ITSx, paste(Sys.Date(), "ITSx_version.txt", sep = "_")))

# Build the full argument vector (each token a separate element!)

# Helper function: convert Windows paths like "D:/dir/file.fastq.gz" to WSL "/mnt/d/dir/file.fastq.gz"
win_to_wsl <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  sub("^([A-Za-z]):/", "/mnt/\\L\\1/", p, perl = TRUE)
}

ITSx.input.file <- file.path(path.5.ITSx, "ITSx.input.fasta")
ITSx.input.file  <- win_to_wsl(ITSx.input.file)
ITSx.input.file

ITS.output.prefix <- file.path(path.5.ITSx, "ITSx_out")
ITS.output.prefix <- win_to_wsl(ITS.output.prefix)


#' # Run ITSx in Ubuntu
#' permit ITSx to flag reads as belonging to other than fungi
#' copy the output (in R console) of `cat(paste(ITSx.args))` and run in Ubuntu/WSL
#' **example:** /home/losl/miniconda3/bin/ITSx -i /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/data/5.ITSx/ITSx.input.fasta -o /mnt/d/Sequencing_Data_iSeq/working_copies/20250618_FS10002142_21_BTR67815-2931/data/5.ITSx/ITSx_out --allow_reorder T --allow_single_domain 1e-3,0 --partial 50
ITSx.args <- c(
  ITSx,
  "-i", ITSx.input.file,
  "-o", ITS.output.prefix,
  "--allow_reorder", "T",
  "--allow_single_domain", "1e-3,0",
  "--partial", 50
)

cat(paste(ITSx.args))

#' # pull out just the fungal reads:
#' copy the output of the following code and paste it into Ubuntu/WSL to change into the ITSx directory

cd.to.ITsx.dir <- c(
  "cd",
  win_to_wsl(path.5.ITSx)
)
cat(paste(cd.to.ITsx.dir))

#' paste the following commented-out code into Ubuntu/WSL to have it extract the fungal reads
# grep -A 1 "|f" ITSx_out.ITS1.full_and_partial.fasta > ITSx.fungal.fasta
# sed -i 's/|/ /' ITSx.fungal.fasta
# sed -i 's/\s.*$//' ITSx.fungal.fasta
# sed -i 's/--//' ITSx.fungal.fasta


fungals <- readDNAStringSet(file.path(path.5.ITSx, "ITSx.fungal.fasta"))
width(fungals)
length(fungals)
# 20250728 63 (merged reads)
# 20250805 55 (non-merged reads)
# 20250805 66 (merged reads)
length(ITSx.input)
# 20250805 190 (merged reads)
length(ITSx.input) - length(fungals)
# 20250805 124 (merged reads)
# 20250805 124 of the original 190 sequence variants were not recognized as fungal ITS reads, 66 categorized as fungals


# 2025-12-23
# > cat(paste(cd.to.ITsx.dir))
# cd /mnt/d/Sequencing_Data_MiSeq/working_copies/251215_M06430_0030_000000000-M9JJFC/SCAP_Faba_Fusarium_R_analysis/data/5.ITSx_merged_reads
# > fungals <- readDNAStringSet(file.path(path.5.ITSx, "ITSx.fungal.fasta"))
# > width(fungals)
# integer(0)
# > length(fungals)
# [1] 0