#` # Assign taxonomy using phyloseq
#' # set file path
path.6.phyloseq <- file.path(here("data", "6.phyloseq_on_all_reads_ITS_trimmed"))
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)

#' ## ASVs (amplicon sequence variants) output from ITS
#' note these are shorter than the ASVs in seqtab (this affects how we process things downstream, because phyloseq only accepts identical matches in the seqtab and taxtab)
#' this ends up being a 1 column dataframe where each row is an ASV
ITS.ASVs <- readDNAStringSet(file.path(here("data","5.ITSx_R1R2_not_merged", "ITSx1.fap.cleanheader.fasta"))) #Amplicon Sequence Variants (ITS.ASVs)
ITS.ASVs.df <- data.frame(ITS.ASVs)
ITS.ASVs.df
rownames(ITS.ASVs.df)
#' ITSx has some duplicated ASVs in its output because ITSx trims to ITS1 (or ITS2) therefore after trimming the ends, some ASVs end up being identical 
#' initially I thought we needed to extract only the unique values, but instead I decided to keep them all since we ultimately use the seqtab ASVs (the full ASV identified by dada) as the rownames for taxtab
#' problem using actual sequences as ASV identifiers, because post-ITSx there are duplicate sequences...

#' # Assign taxonomy against UNITE Fungal Database to each ASV using dada2::assignTaxonomy
#' this step takes a little time
set.seed(1024)
taxtab <- assignTaxonomy(ITS.ASVs, refFasta = classification.ref.ITS, multithread = TRUE, minBoot = 80)
taxtab <- data.frame(taxtab) 
dim(taxtab)
#' rows are sequence variants (ASVs), 7 columns expected (7 taxonomic ranks) 
# 20240812 14 rows (14 sequence variants, with 3 having unassigned kingdoms etc.), 7 columns (taxonomic ranks) for the taxtab.ITS database (20240812 LL)
# 20250728 63 rows (63 sequence variants), 7 columns (taxonomic ranks) for the taxtab.ITS database
# 20250805 on merged oomycetes 66 rows (63 sequence variants), 7 columns (taxonomic ranks) for the taxtab.ITS database

#' save taxtab as R object (since this step takes some time)
save(taxtab, file = file.path(path.6.phyloseq, "taxtab.RData"))

