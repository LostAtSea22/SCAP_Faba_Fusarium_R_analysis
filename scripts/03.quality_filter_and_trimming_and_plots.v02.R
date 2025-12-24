## Filter & trim on quality ####
## Create a filter directory ####
path.3.filter <- here("data", "3.filter")
if(!dir.exists(path.3.filter)) dir.create(path.3.filter)

names.3.filter.R1 <- file.path(path.3.filter, basename(names.2.cut.R1))
names.3.filter.R2 <- file.path(path.3.filter, basename(names.2.cut.R2))
length(names.3.filter.R1)

# number of samples = 14 (20251219 LL)

## Plot quality profiles of 6 random samples to determine filtering and trimming thresholds ####
qual.plot_R1 <- plotQualityProfile(names.2.cut.R1[sample(seq(1:10), 6, replace = FALSE)]) #this takes some time
# plotQualityProfile(names.2.cut.R1[1]) #this works 20250728 LL
# X-axis = base position (cycle)
# Y-axis = quality score (Q)
# Darker blue = median quality
# Orange line = 25th percentile (you care about this for filtering)

print(qual.plot_R1)
path.3.filter
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_6random_R1.pdf", sep = ""), qual.plot_R1)

# 20240812 LL - Looks OK, odd drop around the 20th cycle consistently...
# 20250728 LL - the number of reads per sample were inconcsistent (chalk this up to not repeating the qPCR to get it within the standard curve) otherwise the quality plots were very similar to the previous quality plots where there is a dip in quality around 20 bases into the read and then again a dip for the last 20 bases. Another difference is there are two samples that in the random subset that have multiple dips throughout the cycles... 

qual.plot_R2 <- plotQualityProfile(names.2.cut.R2[sample(seq(1:10), 6, replace = FALSE)]) #did work #picking a random subset of 6 samples, we only have 11. Matt simplified this sample(seq(1:11), 6) ## don't want to look at negative control or mock reference samples here
names.2.cut.R2
print(qual.plot_R2)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_6random_R2.pdf", sep = ""), qual.plot_R2)


# Filter reads based on quality plots ####
filter <- filterAndTrim(names.2.cut.R1, names.3.filter.R1, 
                        names.2.cut.R2, names.3.filter.R2, 
                        maxN = 0,
                        maxEE = c(2, 2),
                        truncLen = c(125, 120), #truncLen: Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded. #this may need to be changed to 121 since this is the max length for the reads?
                        #if truncLen = c(125, 120) = Keep only the first 125 bases from the forward read
                        # Keep only the first 120 bases from the reverse read
                        compress = TRUE) #EE = sum(10^(-Q/10)), therefore maxEE = 2 means that if you sum all the bases Q values the max value that will pass filtering will be a sum of 2 (Q10= +0.1, Q20= +0.01, Q30= +0.001)  
filter

#reads.in reads,out
# reads.in reads.out
# 052-LB1_S14_L001_R1_001.fastq.gz    99495     82954
# 052-LB2_S15_L001_R1_001.fastq.gz   209012     76782
# 052-LB3_S16_L001_R1_001.fastq.gz   108747     90065
# 052-LB4_S17_L001_R1_001.fastq.gz   139446    111895
# 052-LB5_S18_L001_R1_001.fastq.gz   148304    117850
# 056-LB1_S19_L001_R1_001.fastq.gz   114951     92846
# 056-LB2_S20_L001_R1_001.fastq.gz   133483    110224
# 056-LB3_S21_L001_R1_001.fastq.gz   106829     86767
# 056-LB4_S22_L001_R1_001.fastq.gz    27497     19578
# 056-LB5_S23_L001_R1_001.fastq.gz    89609     76976
# 131-LB1_S24_L001_R1_001.fastq.gz    75680     62669
# MC_S26_L001_R1_001.fastq.gz        207612    191717
# PCRnc_S27_L001_R1_001.fastq.gz        194       113
# Xnc-LB1_S25_L001_R1_001.fastq.gz       84        34


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
# names.3.filter.R1 <- names.3.filter.R1[c(1:9,11)]
# names.3.filter.R2 <- names.3.filter.R2[c(1:9,11)]
# names.3.filter.R1
# names.3.filter.R2

# How many reads were culled during the filter step? ####
sum(filter[, 1])
# 1,695,349 reads in (20240718 LL)
# 2,028,311 reads in (20250728)
sum(filter[, 2])
# 1,637,523 reads out (20240718 LL)
# 1,946,885 reads out (20250728)
sum(filter[, 1]) - sum(filter[, 2])
# 57,826 reads culled (20240718 LL)
# 81,426 reads culled (20250728)

print(names.3.filter.R1)
length(names.3.filter.R2)

# replot the Quality Profile to see if filtering improved plots ####
filtered.qual.plot.R1 <- plotQualityProfile(names.3.filter.R1[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R1)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_post_filter_6_random_plot_R1.pdf", sep = ""), filtered.qual.plot.R1, width= 7, height = 7)
# looks a lot better 20250728

filtered.qual.plot.R2 <- plotQualityProfile(names.3.filter.R2[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R2)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_post_filter_6_random_plot_R2.pdf", sep = ""), filtered.qual.plot.R2, width= 7, height = 7)
# 20250728 R2 reads seem to be lower quality, still seeing a dip around 20 cycles/bases in
