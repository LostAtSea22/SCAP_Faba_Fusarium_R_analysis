# This script takes some time.
# This is where I start getting different numbers than Matt's analysis ####
# DADA2 Analysis ####
## learnErrors - Learn and Plot Errors (uses HMM) ####
err.R1 <- learnErrors(names.3.filter.R1, multithread = TRUE, randomize = TRUE) #based on dada2 workflow, we are "learning" the relationship between bases and error rate (looks for unique sequences and if it is identical to a sequence that it coming up 500 times except for 1 base, that 1 base is determined to be an "error")
#103,339,375 total bases in 826,715 reads from 5 samples will be used for learning the error rates. (20240718 LL)
#100,895,875 total bases in 807,167 reads from 6 samples will be used for learning the error rates.20250728
# 20251219 105,832,875 total bases in 846,663 reads from 10 samples will be used for learning the error rates.

plotErrors(err.R1)

# Warning message: In scale_y_log10() : log-10 transformation introduced infinite values. (20240718 LL)
# 20250728:
# Warning messages:
#   1: In scale_y_log10() :
#   log-10 transformation introduced infinite values.
# 2: Removed 82 rows containing missing values or values outside the scale range
# (`geom_line()`).

err.R2 <- learnErrors(names.3.filter.R2, multithread = TRUE, randomize = TRUE) #should we be using set.seed for this? #askmatt
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
derep.R1 <- derepFastq(names.3.filter.R1, verbose = TRUE)
# 20250728 Encountered 4401 unique sequences from 70669 total sequences read. This is only true for Undetermined_S0_L001_R1_001.fastq.gz
derep.R1
names(derep.R1) <- substr(names(derep.R1), 1, 11)
derep.R2 <- derepFastq(names.3.filter.R2, verbose = TRUE)
# 20250728 Encountered 5257 unique sequences from 70669 total sequences read
names(derep.R2) <- substr(names(derep.R2), 1, 11)

### dada - Infer sequence variants ####
dada.R1 <- dada(derep.R1, err = err.R1, multithread = TRUE)  #changes data based on likely errors, 
dada.R1[[5]] #look at sample 5 to see the result
# 20240718 496 sequence variants were inferred from 20926 input unique sequences. Same number of input seq as Matt (20240718 LL)
# 20250728 72 sequence variants were inferred from 2741 input unique sequences.
# 20250801 70 sequence variants were inferred from 2741 input unique sequences.
# 20250801 Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

dada.R2 <- dada(derep.R2, err = err.R2, multithread = TRUE)
dada.R2[[5]]
# 20240718 450 sequence variants were inferred from 21488 input unique sequences. Same number of input seq as Matt (20240718 LL)
# 20250728 61 sequence variants were inferred from 3358 input unique sequences.
# 20250801 63 sequence variants were inferred from 3358 input unique sequences.
# 20250801 Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# ideally, would track whether / how many reads are culled here. culled = discarded (from animals not chosen for breeding)

path.4.dada <- here("data", "4.dada")
if(!dir.exists(path.4.dada)) dir.create(path.4.dada)

saveRDS(derep.R1, here("data", "4.dada", "derep.R1.rds"))
saveRDS(derep.R2, here("data", "4.dada", "derep.R2.rds"))

saveRDS(dada.R1, here("data", "4.dada", "dada.R1.rds"))
saveRDS(dada.R2, here("data", "4.dada", "dada.R2.rds"))
