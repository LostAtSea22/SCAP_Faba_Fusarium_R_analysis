# Import dada2 objects ####
derep.R1 <- readRDS(here("data", "4.dada", "derep.R1.rds"))
derep.R2 <- readRDS(here("data", "4.dada", "derep.R2.rds"))
dada.R1 <- readRDS(here("data", "4.dada", "dada.R1.rds"))
dada.R2 <- readRDS(here("data", "4.dada", "dada.R2.rds"))

# mergePairs - merge R1 & R2 #################
merge <- mergePairs(dada.R1, derep.R1, #this may fail because iseq length. nmatch is present therefore looks like it works
                    dada.R2, derep.R2,
                    maxMismatch = 2, verbose = TRUE) # in this protocol we wait to merge R1 and R2, keeping quality info for much longer, and only merging Read1 and read2 after you have corrected for sequence errors.
head(merge[[1]], 15) #print the first 15 lines of merge
head(merge[[1]], 1)

# 20260107 LL
# 31266 paired-reads (in 677 unique pairings) successfully merged out of 59401 (in 1510 pairings) input.
# Duplicate sequences in merged output.
# 29060 paired-reads (in 696 unique pairings) successfully merged out of 55144 (in 1564 pairings) input.
# Duplicate sequences in merged output.
# 33274 paired-reads (in 735 unique pairings) successfully merged out of 58118 (in 1682 pairings) input.
# Duplicate sequences in merged output.
# 42366 paired-reads (in 942 unique pairings) successfully merged out of 75259 (in 2152 pairings) input.
# Duplicate sequences in merged output.
# 35748 paired-reads (in 774 unique pairings) successfully merged out of 81871 (in 1904 pairings) input.
# Duplicate sequences in merged output.
# 32429 paired-reads (in 814 unique pairings) successfully merged out of 61478 (in 1772 pairings) input.
# Duplicate sequences in merged output.
# 41853 paired-reads (in 820 unique pairings) successfully merged out of 77122 (in 1899 pairings) input.
# Duplicate sequences in merged output.
# 30693 paired-reads (in 711 unique pairings) successfully merged out of 61330 (in 1866 pairings) input.
# Duplicate sequences in merged output.
# 5507 paired-reads (in 153 unique pairings) successfully merged out of 12114 (in 352 pairings) input.
# Duplicate sequences in merged output.
# 21650 paired-reads (in 519 unique pairings) successfully merged out of 62809 (in 1376 pairings) input.
# Duplicate sequences in merged output.
# 15240 paired-reads (in 624 unique pairings) successfully merged out of 47798 (in 1603 pairings) input.
# Duplicate sequences in merged output.
# 33 paired-reads (in 10 unique pairings) successfully merged out of 183480 (in 2035 pairings) input.
# 106 paired-reads (in 2 unique pairings) successfully merged out of 106 (in 2 pairings) input.
# 24 paired-reads (in 1 unique pairings) successfully merged out of 24 (in 1 pairings) input.



# 63254 paired-reads (in 217 unique pairings) successfully merged out of 96976 (in 269 pairings) input.
# Duplicate sequences in merged output.
# 64627 paired-reads (in 190 unique pairings) successfully merged out of 104679 (in 236 pairings) input.
# Duplicate sequences in merged output.
# 84855 paired-reads (in 162 unique pairings) successfully merged out of 93318 (in 226 pairings) input.
# Duplicate sequences in merged output.
# 79924 paired-reads (in 131 unique pairings) successfully merged out of 92543 (in 208 pairings) input.
# Duplicate sequences in merged output.
# 92260 paired-reads (in 86 unique pairings) successfully merged out of 95584 (in 128 pairings) input.
# Duplicate sequences in merged output.
# 102727 paired-reads (in 58 unique pairings) successfully merged out of 116043 (in 104 pairings) input.
# Duplicate sequences in merged output.
# 168626 paired-reads (in 19 unique pairings) successfully merged out of 168657 (in 20 pairings) input.
# 294025 paired-reads (in 20 unique pairings) successfully merged out of 294080 (in 23 pairings) input.
# Duplicate sequences in merged output.
# 152766 paired-reads (in 58 unique pairings) successfully merged out of 170274 (in 113 pairings) input.
# 162860 paired-reads (in 59 unique pairings) successfully merged out of 180271 (in 131 pairings) input.
# 6286 paired-reads (in 5 unique pairings) successfully merged out of 6376 (in 6 pairings) input.
# 111402 paired-reads (in 2 unique pairings) successfully merged out of 111405 (in 3 pairings) input.
# 129 paired-reads (in 7 unique pairings) successfully merged out of 131831 (in 43 pairings) input.
# 258 paired-reads (in 17 unique pairings) successfully merged out of 212647 (in 76 pairings) input.
# 52305 paired-reads (in 107 unique pairings) successfully merged out of 70460 (in 196 pairings) input.
# Duplicate sequences in merged output.

# for “63,254 paired-reads … merged out of 96,976 input”
  # ~65% of your passing, denoised F/R pairs overlapped cleanly and produced a merged amplicon. The rest didn’t merge (usually insufficient overlap or mismatches in the overlap).
# for “in 217 unique pairings … out of 269 pairings”
  # You had 269 distinct forward-ASV × reverse-ASV combinations present; 217 of those combinations yielded a valid merged sequence. The others failed for the reasons above.

## Make sequence tables from merged read pairs ####
seqtab.merge <- makeSequenceTable(merge) #makeSequenceTable is a dada2 function that generates a sequence table from dada2 object
dim(seqtab.merge) #dimensions of an object, if a table it will return the vector c(nrow, ncol), ie number of rows and number of columns.
# 20240718 10 samples X 2,489 sequence variants
# 20250728 15 samples x 386 sequence variants
# 20250801 15 samples x 707 sequence variants
# 20260107 14 samples x 5372 sequence variants

sum(seqtab.merge)
# 20240718   797,006 reads (20240718 LL)
# 20250728 1,946,112 reads 
# 20250801 1,436,304 reads
# 20260107   319,249 reads (not very many... based on this and the information above, I think many of my reads failed to overlap)

unname(seqtab.merge[, 1:10]) # removes the column names from the first 10 columns of the sequence table
seq.merge.len.distr.tab <- table(nchar(getSequences(seqtab.merge))) #provides a frequency table of the lengths of the sequences in the sequence table
# getSequences() is a dada2 function, that extracts the sequences (the unique sequence variants, or ASVs),
# nchar() function calculates the number of characters in a string, 
# table() function creates a frequency table (ie. how many sequences have each possible length).
barplot(seq.merge.len.distr.tab, main = "Frequency of Sequence Lengths", xlab = "Sequence Length", ylab = "Frequency", col = "blue")
pdf(file = here("data", paste(Sys.Date(), "Freq_Seq_Length_from_seqtab_merge.pdf", sep = "")))
# 20250728 all the sequences are 125 bases long - I think this is because I didnt merge the reads because there was no overlap
# 20250728 sequences vary from 125-233 bases long because of merge step
# 20260107 sequence length varies from 225-335, nothing close to 500bp which is the expected amplicon length... therefore I am not confident about the overlapping
# 20260107 there is also a large bar at 225 and ~330, this makes sense because the forward reads were trimmed to 225nt, and the reverse reads were trimmed to 125, therefore no reads could actually be longer than 350
write.csv(seqtab.merge, file = here("data", "4.dada", paste(Sys.Date(), "SeqTable_mereged_pairs.csv", sep = "")))

## Make sequence tables from dada.R1 (ie just forward reads - not merged) ####
seqtab.R1 <- makeSequenceTable(dada.R1) #makeSequenceTable is a dada2 function that generates a sequence table from dada2 object
#20240718 dada.R1 used instead of merge because sequences did not merge well
dim(seqtab.R1) #dimensions of an object, if a table it will return the vector c(nrow, ncol), ie number of rows and number of columns.
# 20240718 10 samples X 2,489 sequence variants (20240718 LL)
# 20250728 15 samples x 386 sequence variants
# 20250801 15 samples x 386 sequence variants
# 20260107 14 samples x 11,513 sequence variants (now this is a number closer to what I would expect!)

sum(seqtab.R1)
# 20240718   797,006 reads (20240718 LL)
# 20250728 1,946,112 reads 
# 20250801 1,946,112 reads 
# 20260107   844,932 reads

unname(seqtab.R1[, 1:10]) # removes the column names from the first 10 columns of the sequence table
seq.R1.len.distr.tab <- table(nchar(getSequences(seqtab.R1))) #provides a frequency table of the lengths of the sequences in the sequence table
plot.R1 <- barplot(seq.R1.len.distr.tab, main = "Frequency of Sequence Length", xlab = "Sequence Length", ylab = "Frequency", col = "blue")
plot.R1 #returns just a printout of 0.7 20260107
str(seq.R1.len.distr.tab) #shows that all sequences are a length of 225 20260107
pdf(file = here("data", paste(Sys.Date(), "Freq_Seq_Length_from_seqtab_R1.pdf", sep = "")))
# 20250728 all the sequences are 125 bases long - I think this is because I didnt merge the reads because there was no overlap


## Make sequence table from dada.R2 (ie just reverse reads - not merged) ####
seqtab.R2 <- makeSequenceTable(dada.R2)
dim(seqtab.R2)
# 20250801 15 347 ASVs
# 20260107 14 9555 ASVs
sum(seqtab.R2)
# 20250801 1,945,715
# 20260107   846,175

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

## Remove chimeras from merged reads ##################################
seqtab.merge.chimera <- removeBimeraDenovo(seqtab.merge, multithread = TRUE, verbose = TRUE, method = "consensus")
# 20240808 Identified 1,621 bimeras out of 2,489 input sequences. (20240808 LL) comapred to 1,666 bimeras out of 2,537 input sequences identified by MB
# 20250728 Identified 515 bimeras out of 701 input sequences.
# 20250801 Identified 517 bimeras out of 707 input sequences.
# 20260107 Identified 329 bimeras out of 5372 input sequences 6.12%

dim(seqtab.merge)

dim(seqtab.merge.chimera)
# 20240808 10 samples X 868 sequence variants (ASVs)
# 20250728 15 samples x 186 sequence variants
# 20250801 15 samples x 190 sequence variants
# 20260107 14 samples x 5043 sequence variants

sum(seqtab.merge.chimera)
# 20240808   716,317 reads
# 20250728 1,384,835 reads
# 20250801 1,384,827 reads
# 20260107   305,371 reads

1 - (sum(seqtab.merge.chimera) / sum(seqtab.merge))
# 20240808 0.1012401 proportion of chimeric reads #AskMAtt is this high? this seems kinda high (10%)
# 20250728 0.03583031 proportion of chimeric reads - this doesnt seem too bad to me LL
# 20250801 0.03583991 proportion of chimeric reads - this doesnt seem too bad to me LL
# 20260107 0.04347077 proportion of chimeric reads - not bad


write.csv(seqtab.merge.chimera, file = here("data", "4.dada", paste(Sys.Date(), "SeqTable_chimeras_rm_merged_pairs.csv", sep = "")))

## Remove chimeras from R1 reads #############################
seqtab.R1.chimera <- removeBimeraDenovo(seqtab.R1, multithread = TRUE, verbose = TRUE, method = "consensus")
# 20250801 Identified 155 bimeras out of 384 input sequences.
# 20260107 Identified 253 bimeras out of 11513 input sequences

dim(seqtab.R1.chimera)
# 20240808 10 samples X 868 sequence variants (20240808 LL, rows are samples, columns are unique sequences)
# 20250728 15 samples x 186 sequence variants
# 20260107 14 11260

sum(seqtab.R1.chimera)
# 20250801 1,902,730 reads
# 20260107 795971

1 - (sum(seqtab.R1.chimera) / sum(seqtab.R1))
# 20250801 0.02223887 proportion of chimeric reads
# 20260107 0.05794667

write.csv(seqtab.R1.chimera, file = here("data", "4.dada", paste(Sys.Date(), "SeqTable_chimeras_rm_R1.csv", sep = "")))

## Remove chimeras from R2 reads ###############################
seqtab.R2.chimera <- removeBimeraDenovo(seqtab.R2,
                                     multithread = TRUE, verbose = TRUE, method = "consensus")
# 20250801 Identified 148 bimeras out of 347 input sequences.
# 20260107 Identified 129 bimeras out of 9555 input sequences

dim(seqtab.R2.chimera)
# 20250801 15 samples x 199 sequence variants
# 20260107 14           9,426

sum(seqtab.R2.chimera)
# 20250801 1,919,694 reads
# 20260107 822,790

1 - (sum(seqtab.R2.chimera) / sum(seqtab.R2))
# 20250801 0.01337349 proportion of chimeric reads
# 20260107 0.02763613

write.csv(seqtab.R2.chimera, file = here("data", "4.dada", paste(Sys.Date(), "SeqTable_chimeras_rm_R2.csv", sep = "")))

## join the forwards and reverse reads into a large dataframe so that I can assign phylogenies to all of them ##################
ncol(seqtab.R1.chimera)
# 229
# 20260107 11260
ncol(seqtab.R2.chimera)
# 199
# 20260107 9426

seqtab.R1.R2.chimera <- cbind(seqtab.R1.chimera, seqtab.R2.chimera)
ncol(seqtab.R1.R2.chimera)
# 428
# 20260107 20686

sum(ncol(seqtab.R1.chimera), ncol(seqtab.R2.chimera))
# 20260107 20686

write.csv(seqtab.R1.R2.chimera, file = here("data", "4.dada", paste(Sys.Date(), "SeqTable_chimeras_rm_R1R2.csv", sep = "")))

