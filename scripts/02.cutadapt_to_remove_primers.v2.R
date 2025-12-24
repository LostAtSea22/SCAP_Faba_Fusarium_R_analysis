# create output folder for cutadapt:
path.2.cut <- here("data", "2.cut")
if(!dir.exists(path.2.cut)) dir.create(path.2.cut)
names.2.cut.R1 <- file.path(path.2.cut, basename(names.1.data.R1))
names.2.cut.R2 <- file.path(path.2.cut, basename(names.1.data.R2))

# Helper: convert Windows paths like "D:/dir/file.fastq.gz" to WSL "/mnt/d/dir/file.fastq.gz"
win_to_wsl <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  sub("^([A-Za-z]):/", "/mnt/\\L\\1/", p, perl = TRUE)
}

# Check cutadapt is reachable via WSL (should print a version line)
cutadapt <- "/usr/bin/cutadapt"
cat(paste(system2("wsl", args = c(cutadapt, "--version"), stdout = TRUE, stderr = TRUE), collapse = "\n"), "\n")

# # Primer sequences

# Identify primer sequences (20251219 LL)
f <- "GTCATCGGCCACGTCGACTCTGG" #fefF
f.rc <- dada2:::rc(f)
r <- "CCTTDCCGAGCTCRGCGGCTTCC" #fefR
r.rc <- dada2:::rc(r)

f
f.rc
r
r.rc

#' use the forward and reverse primer sequences in cutadapt

flags.R1 <- paste("-g", f, "-a", r.rc) 
flags.R2 <- paste("-G", r, "-A", f.rc) 

# Make sure your file paths are Linux-style for WSL
inR1  <- win_to_wsl(names.1.data.R1)
inR2  <- win_to_wsl(names.1.data.R2)
outR1 <- win_to_wsl(names.2.cut.R1)
outR2 <- win_to_wsl(names.2.cut.R2)

# Optional: where to store per-sample logs (on Windows side, no conversion needed)
log_dir <- here::here("logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

for (i in seq_along(inR1)) {
  msg <- sprintf("[%d/%d] Running cutadapt on %s", i, length(inR1), basename(inR1[i]))
  message(msg)
  
  # Build the full argument vector (each token a separate element!)
  args <- c(
    cutadapt,
    flags.R1,               # R1: 5' and 3' adapters
    flags.R2,               # R2: 5' and 3' adapters
    "-n", 2,                #Remove up to "-n" adapters from each read
    "--discard-untrimmed",  #Discard reads that do not contain an adapter
    "--pair-filter=any",    #Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
    "--minimum-length", 50, #Discard reads shorter than 50, quality filtering
    "--max-n", 0,           # 0 means discard reads containing any N's
    "-o", outR1[i],         # output file names
    "-p", outR2[i],
    inR1[i], inR2[i]        # input file names
  )
  
  # Capture both stdout (contains the JSON report) and stderr (progress + messages)
  # stdout_txt <- tempfile(fileext = ".json")
  stderr_txt <- file.path(log_dir, sprintf("cutadapt_%s.stderr.txt", tools::file_path_sans_ext(basename(names.2.cut.R1[i]))))
  
  res <- tryCatch(
    {
      out <- system2("wsl", args = args, stdout = TRUE, stderr = TRUE)
      attr(out, "status") <- if (is.null(attr(out, "status"))) 0L else attr(out, "status")
      list(out = out, status = attr(out, "status"))
    },
    error = function(e) list(out = paste("ERROR:", conditionMessage(e)), status = 1L)
  )
  
  # Write logs
  writeLines(res$out, con = stderr_txt)
  
  # Check return code and output files
  if (res$status != 0L || !file.exists(names.2.cut.R1[i]) || !file.exists(names.2.cut.R2[i])) {
    stop(sprintf("Cutadapt failed for %s; see %s", basename(inR1[i]), stderr_txt))
  }
  
  # Quick success indicator (file sizes)
  s1 <- file.info(names.2.cut.R1[i])$size
  s2 <- file.info(names.2.cut.R2[i])$size
  message(sprintf("  -> OK. Output sizes: R1=%s bytes, R2=%s bytes", format(s1, big.mark=","), format(s2, big.mark=",")))
  flush.console()
}

# Count number of reads in which the primer is found BEFORE cutadapt applied (this step takes some time)
primerHits <- function(primer, names) {
  nhits <- vcountPattern(primer, sread(readFastq(names)), fixed = FALSE, max.mismatch = 1)
  return(sum(nhits > 0))
}

r1s <- sapply(c(f, r.rc), primerHits, names = names.1.data.R1)
r1s
# GCGGAAGGATCATTACCAC GCTCGCACAHCGATGAAGA (20250725 LL)
#             2125006                   18704 

# GTCATCGGCCACGTCGACTCTGG GGAAGCCGCYGAGCTCGGHAAGG (20251219_LL)
#             1127082                  451021 

r2s <- sapply(c(r, f.rc), primerHits, names = names.1.data.R2)
r2s
# TCTTCATCGDTGTGCGAGC GTGGTAATGATCCTTCCGC (20250725 LL)
#             3661001               16290 

# CCTTDCCGAGCTCRGCGGCTTCC CCAGAGTCGACGTGGCCGATGAC (20251219_LL)
#             1584023               219019 

# Count number of reads in which the primer is found AFTER cutadapt applied (this step takes some time)
sapply(c(f, r.rc), primerHits, names = names.2.cut.R1)
# GCGGAAGGATCATTACCAC GCTCGCACAHCGATGAAGA (20250725 LL)
#                   0                   0 

# GTCATCGGCCACGTCGACTCTGG GGAAGCCGCYGAGCTCGGHAAGG (20251219_LL)
# 1                       3 
sapply(c(r, f.rc), primerHits, names = names.2.cut.R2)
# TCTTCATCGDTGTGCGAGC GTGGTAATGATCCTTCCGC (20250725 LL)
#                   0                   0 

# CCTTDCCGAGCTCRGCGGCTTCC CCAGAGTCGACGTGGCCGATGAC (20251219_LL)
# 0                       1 
