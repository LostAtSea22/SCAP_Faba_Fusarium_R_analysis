#' # Assign taxonomy using phyloseq

#' # set file path
path.6.phyloseq <- file.path(here("data", "6.phyloseq_on_all_merged_reads"))
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)

#' ## Import ASVs (amplicon sequence variants)
#' I wanted to assign taxonomy to all reads
#' therefore I will extract the ASVs from the column headings of the dada2 table (seqtab)
seqtab <- read.csv(file.path(here("data", "4.dada", "2026-01-07SeqTable_chimeras_rm_merged_pairs.csv")))
seqtab <- t(seqtab)
seqtab <- as.data.frame(seqtab)
seqtab$ASV <- rownames(seqtab)
rownames(seqtab) <- NULL
ASVs.df <- seqtab[-1,"ASV"]
#' this ends up being a 1 column dataframe where each row is an ASV

#' # Assign taxonomy against UNITE Fungal Database to each ASV using dada2::assignTaxonomy
#' this step takes a little time
set.seed(1024)
taxtab <- assignTaxonomy(ASVs.df, refFasta = classification.ref.ITS, multithread = TRUE, minBoot = 80)
#' **NOTE**: the variarble, `classification.ref.ITS` is loaded in the script `00.loading_packages.v2.R`

taxtab <- data.frame(taxtab) 
dim(taxtab)
#' rows are sequence variants (ASVs), 7 columns expected (7 taxonomic ranks) 
print(paste(Sys.Date(), "-", dim(taxtab)[1], "sequence variants (ASVs) and", dim(taxtab)[2], "taxonomic ranks"))
# 2026-01-07 - 5043 sequence variants (ASVs) and 7 taxonomic ranks

#' save taxtab as R object (since this step takes some time)
save(taxtab, file = file.path(path.6.phyloseq, "taxtab.RData"))
