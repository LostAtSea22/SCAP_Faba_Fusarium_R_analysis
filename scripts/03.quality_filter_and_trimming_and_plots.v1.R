## Filter & trim on quality ####
## Create a filter directory ####
path.3.filter <- here("data", "3.filter")
if(!dir.exists(path.3.filter)) dir.create(path.3.filter)

names.3.filter.R1 <- file.path(path.3.filter, basename(names.2.cut.R1))
names.3.filter.R2 <- file.path(path.3.filter, basename(names.2.cut.R2))
length(names.3.filter.R1)

# number of samples = 15 (20250725 LL)
# number of samples = 14 (20260106 LL)

## Plot quality profiles of 6 random samples to determine filtering and trimming thresholds ####
qual.plot_R1 <- plotQualityProfile(names.2.cut.R1[sample(seq(1:10), 6, replace = FALSE)]) #this takes some time
# plotQualityProfile(names.2.cut.R1[1]) #this works 20250728 LL
# X-axis = base position (cycle)
# Y-axis = quality score (Q)
# Darker blue = median quality
# Orange line = 25th percentile (you care about this for filtering)

print(qual.plot_R1)
path.3.filter
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_6random_R1_cut1.pdf", sep = ""), qual.plot_R1)

# 20240812 LL - Looks OK, odd drop around the 20th cycle consistently...
# 20250728 LL - the number of reads per sample were inconcsistent (chalk this up to not repeating the qPCR to get it within the standard curve) otherwise the quality plots were very similar to the previous quality plots where there is a dip in quality around 20 bases into the read and then again a dip for the last 20 bases. Another difference is there are two samples that in the random subset that have multiple dips throughout the cycles... 
# 20260106 LL - looks pretty rough after 200bp

## CUTADAPT RERUN OUTPUT Plot quality profiles of 6 random samples to determine filtering and trimming thresholds ####
qual.plot_R1_cut2 <- plotQualityProfile(names.2.cut2.R1[sample(seq(1:10), 6, replace = FALSE)]) #this takes some time

print(qual.plot_R1_cut2)
path.3.filter
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_6random_R1_cut2.pdf", sep = ""), qual.plot_R1_cut2)

# 20240812 LL - Looks OK, odd drop around the 20th cycle consistently...
# 20250728 LL - the number of reads per sample were inconcsistent (chalk this up to not repeating the qPCR to get it within the standard curve) otherwise the quality plots were very similar to the previous quality plots where there is a dip in quality around 20 bases into the read and then again a dip for the last 20 bases. Another difference is there are two samples that in the random subset that have multiple dips throughout the cycles... 
# 20260106 LL - looks pretty rough after 200bp

qual.plot_R2 <- plotQualityProfile(names.2.cut.R2[sample(seq(1:10), 6, replace = FALSE)]) #did work #picking a random subset of 6 samples, we only have 11. Matt simplified this sample(seq(1:11), 6) ## don't want to look at negative control or mock reference samples here
names.2.cut.R2
print(qual.plot_R2)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_6random_R2.pdf", sep = ""), qual.plot_R2)

# 20260106 THOSE look even worse...

# Filter reads based on quality plots ####
filter <- filterAndTrim(names.2.cut.R1, names.3.filter.R1, 
                        names.2.cut.R2, names.3.filter.R2, 
                        maxN = 0,
                        maxEE = c(2, 2),
                        truncLen = c(225, 125), #truncLen: Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded. #this may need to be changed to 121 since this is the max length for the reads?
                        #if truncLen = c(125, 120) = Keep only the first 125 bases from the forward read
                        # Keep only the first 120 bases from the reverse read
                        # determined by quality plots
                        compress = TRUE) #EE = sum(10^(-Q/10)), therefore maxEE = 2 means that if you sum all the bases Q values the max value that will pass filtering will be a sum of 2 (Q10= +0.1, Q20= +0.01, Q30= +0.001)  
filter

# 20260106 LL
# reads.in reads.out
# 052-LB1_S14_L001_R1_001.fastq.gz    99495     62520
# 052-LB2_S15_L001_R1_001.fastq.gz   209012     58090
# 052-LB3_S16_L001_R1_001.fastq.gz   108747     61423
# 052-LB4_S17_L001_R1_001.fastq.gz   139446     79271
# 052-LB5_S18_L001_R1_001.fastq.gz   148304     85513
# 056-LB1_S19_L001_R1_001.fastq.gz   114951     64778
# 056-LB2_S20_L001_R1_001.fastq.gz   133483     80269
# 056-LB3_S21_L001_R1_001.fastq.gz   106829     64743
# 056-LB4_S22_L001_R1_001.fastq.gz    27497     13658
# 056-LB5_S23_L001_R1_001.fastq.gz    89609     65484
# 131-LB1_S24_L001_R1_001.fastq.gz    75680     50709
# MC_S26_L001_R1_001.fastq.gz        207612    183750
# PCRnc_S27_L001_R1_001.fastq.gz        194       106
# Xnc-LB1_S25_L001_R1_001.fastq.gz       84        28

#20250728 hmmm my negative control still has a lot of reads... thats no good :(
#20250728 therefore I will skip the following step

# 20260106 handful of reads in the negative control - curious if there was some contamintaion in the extraction kit
# 20260106 therefore I will skip the following step

# negative control sample reduced to 0 reads, no filtered file created because all reads removed (updated 20240718 LL)
# remove NTC from file name list
# names.3.filter.R1 <- names.3.filter.R1[c(1:9,11)]
# names.3.filter.R2 <- names.3.filter.R2[c(1:9,11)]
# names.3.filter.R1
# names.3.filter.R2

# How many reads were culled during the filter step? ####
sum(filter[, 1])
# 1,695,349 reads in (20240718 LL)
# 2,028,311 reads in (20250728 LL)
# 1,460,943 reads in (20260106 LL)

sum(filter[, 2])
# 1,637,523 reads out (20240718 LL)
# 1,946,885 reads out (20250728)
#   800,380 reads out (20260106 LL) #lots of reads filtered.. concerning
#   870,342 reads out (20260106 LL) #lots of reads filtered.. concerning

sum(filter[, 1]) - sum(filter[, 2])
#  57,826 reads culled (20240718 LL)
#  81,426 reads culled (20250728)
# 660,563 reads culled (20260106 LL)
# 590,601 reads culled (20260106 LL)

# replot the Quality Profile to see if filtering improved plots ####
filtered.qual.plot.R1 <- plotQualityProfile(names.3.filter.R1[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R1)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_post_filter_6_random_plot_R1.pdf", sep = ""), filtered.qual.plot.R1, width= 7, height = 7)
# looks a lot better 20250728
# looks quite a big better 20260106 LL except Xnc

filtered.qual.plot.R2 <- plotQualityProfile(names.3.filter.R2[sample(seq(1:15), 6, replace = FALSE)])
print(filtered.qual.plot.R2)
ggsave(paste(path.3.filter, "/", Sys.Date(), "plotQualityProfile_post_filter_6_random_plot_R2.pdf", sep = ""), filtered.qual.plot.R2, width= 7, height = 7)

# 20250728 R2 reads seem to be lower quality, still seeing a dip around 20 cycles/bases in
# 20260106 R2 do appear lower


# Export objects required downstream --------------------------------------
saveRDS(names.3.filter.R1, file.path(here("data", "3.filter"), "names.3.filter.R1.RDS"))
saveRDS(names.3.filter.R2, file.path(here("data", "3.filter"), "names.3.filter.R2.RDS"))

