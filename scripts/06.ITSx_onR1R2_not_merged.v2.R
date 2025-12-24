#' date: 2025-08-05
#' author: "Laura Los"
#' title: "06.ITSx trimming conserved flanking regions for Forwards and reverse reads (not merged)"
#' version_notes: "

#------------

## Trim off conserved flanking regions ####
  
# create folder
path.5.ITSx <- file.path(here("data", "5.ITSx_R1R2_not_merged"))
if(!dir.exists(path.5.ITSx)) dir.create(path.5.ITSx)

#import sequence tables files
seqtab.R1.R2.chimera.rm <- read.csv(here("data", "4.dada", "2025-08-01SeqTable_chimeras_rm_R1R2.csv"))
rownames(seqtab.R1.R2.chimera.rm) <- seqtab.R1.R2.chimera.rm$X
seqtab.R1.R2.chimera.rm <- seqtab.R1.R2.chimera.rm[,-1]
rownames(seqtab.R1.R2.chimera.rm)
colnames(seqtab.R1.R2.chimera.rm[,1:10])

#prepare data for ITSx
ITSx.input <- DNAStringSet(colnames(seqtab.R1.R2.chimera.rm))
ASV.names <- seq(dim(seqtab.R1.R2.chimera.rm)[2])
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

#' # pull out just the oomycete reads and clean fasta headers:
#' copy the output of the following code and paste it into Ubuntu/WSL to change into the ITSx directory

cd.to.ITsx.dir <- c(
  "cd",
  win_to_wsl(path.5.ITSx)
)
cat(paste(cd.to.ITsx.dir))

#' paste the following commented-out code into Ubuntu/WSL to have it extract the oomycete reads
# grep -A 1 "|O" ITSx_out.ITS1.full_and_partial.fasta > ITSx.oomycete.fasta
# sed -i 's/|/ /' ITSx.oomycete.fasta
# sed -i 's/\s.*$//' ITSx.oomycete.fasta
# sed -i 's/--//' ITSx.oomycete.fasta

#' paste the following commented-out code into Ubuntu/WSL to have it remove excess text in the header
cp ITSx_out.ITS1.full_and_partial.fasta ITSx1.fap.cleanheader.fasta
sed -i 's/|/ /' ITSx1.fap.cleanheader.fasta
sed -i 's/\s.*$//' ITSx1.fap.cleanheader.fasta
sed -i 's/--//' ITSx1.fap.cleanheader.fasta

oomycetes <- readDNAStringSet(file.path(path.5.ITSx, "ITSx.oomycete.fasta"))
width(oomycetes)
length(oomycetes)
# 20250728 63 (merged reads)
# 20250805 55 (non-merged reads)
length(ITSx.input)
#20250805 428
length(ITSx.input) - length(oomycetes)
#20250805 373
# 854 of the original 868 sequence variants were not recognizable as oomycetes ITS reads, 14 categorized as oomycetes (20240812 LL)
# 20250728 123 of the original 186 sequence variants were not recognized as oomycete ITS reads, 63 categorized as oomycetes
# 20250805 373 of the original 428 sequence variants were not recognized as oomycete ITS reads, 55 categorized as oomycetes