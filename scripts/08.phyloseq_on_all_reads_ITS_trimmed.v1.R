#` # Create Abundance (OTU/ASV) Table using phyloseq
#' **NOTE**: the variarble, `classification.ref.ITS` is assigned in the script `00.loading_packages.v2.R`
#' # set file path
path.6.phyloseq <- file.path(here("data", "6.phyloseq_on_all_reads_ITS_trimmed"))
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)

#' # Import data

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

#' import the Sequence Table which is essentially a matrix where the rownames are your sample names, each column is an ASV and the values are the relative counts or hits for each ASV in the sample
#' note the seqtab contains all ASVs assigned by dada2 with chimeras removed, from mereged pairs
seqtab <- read.csv(file.path(here("data", "4.dada", "2025-08-01SeqTable_chimeras_rm_R1R2.csv")))
seqtab # the column X contains sample names, need to make this the row names
rownames(seqtab) <- seqtab[,"X"] #make column X values the rownames
seqtab <- select(seqtab, -X) #remove column X
dim(seqtab)
print(paste(Sys.Date(), "-", dim(seqtab)[1],"samples and", dim(seqtab)[2], "sequence variants"))

#' # load taxonomy table created from ITSx ASVs analysed against UNITE Fungal Database using dada2::assignTaxonomy
load(file.path(path.6.phyloseq, "taxtab.RData"))
dim(taxtab)
print(paste(Sys.Date(), "-", dim(taxtab)[1],"sequence variants (ASVs) and", dim(taxtab)[2], "taxonomic ranks"))
#' rows are sequence variants (ASVs), 7 columns expected (7 taxonomic ranks) 

#' # Extract ASVs from the SeqTable based on ASVs in ITSx (ex. if the ITSx file contains only oomycete ASVs, kind of a filtering step) 
#' this step can be skipped if you are interested in complete taxonomy in the sample
#' this was the paste script:
seqtab.ITS.ASVs <- seqtab[, as.numeric(rownames(ITS.ASVs.df))]
#' row names of ITS.ASVs.df are numbers that reference the corresponding ASVS column from seqtab (made earlier on in dada2 step, when removing chimeras)
#' this step creates a subset of seqtab that corresponds to the ASVs that were filtered during the ITS step (ex. oomycetes)
dim(seqtab.ITS.ASVs)
print(paste(Sys.Date(), "-", dim(seqtab.ITS.ASVs)[1],"samples and", dim(seqtab.ITS.ASVs)[2], "sequence variants from ITS filtering step"))


#' # Remove Negative Control in data 
# metadata <- metadata[-10, ] 
# 20240812 removing the negative control I had in my data (which was row 10) (20240812 LL)
# 20250728 I commented this out because I am worried about contamination in my negative control so I'm leaving it in 

#' # Phyloseq
#' to use phyloseq we need:
#'    - all the sample names to match for the otu table (seqtab.ITS.ASVs), the tax table (taxtab) and the sample data (metadata)
#'    - **OTU/ASV** Table (or Abundance Table):
#'        - This table contains the counts or abundances of operational taxonomic units (OTUs) or amplicon sequence variants (ASVs) across different samples.
#'        - It should be a numeric matrix where rows represent taxa (OTUs/ASVs) and columns represent samples, or vice versa, depending on how it's imported (specified by taxa_are_rows argument in otu_table()).
#'    - **Taxonomy Table**:
#'        - This file provides the taxonomic classification for each OTU/ASV present in the abundance table.
#'        - It should be a character matrix with rows named by the OTU/ASV identifiers and columns named by the taxonomic ranks (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species).
#'    - **Sample Data** (Metadata) File:
#'        - This file contains descriptive information about each sample, such as experimental conditions, environmental factors, or subject characteristics
#'        - It should be a data frame with rows named by sample identifiers and columns representing various metadata categories.
#'    - Phylogenetic Tree (Optional but Recommended):
#'        - A phylogenetic tree representing the evolutionary relationships between the OTUs/ASVs can be included for analyses that rely on phylogenetic distances.
#'        - This is typically in Newick format.

#' Preparing ASV Abundance Table
seqtab.ITS.ASVs <- as.matrix(seqtab.ITS.ASVs)
seqtab.ITS.ASVs <- t(seqtab.ITS.ASVs)
str(seqtab.ITS.ASVs)

#' remove undetermined reads (reads that illuina was unable to assign to a specific index, including PhiX reads) from data (probably should have done this earlier on in the analysis???)
colnames(seqtab.ITS.ASVs)               
seqtab.ITS.ASVs <- subset(seqtab.ITS.ASVs, select = -c(Undetermine))
colnames(seqtab.ITS.ASVs)

#' Preparing ASV Taxonomy Table (for tax_table)
#' the ASVs in taxtab are trimmed versions of the ASVs in seqtab.ITS.ASVs (becasue ITSx trims to ITS1 or ITS2)
#' but phyloseq requires the two matrices to have identical rownames therefore we are assigning the seqtab.ITS.ASVs ASVs (rownames) to the taxtab rownames
rownames(taxtab) <- rownames(seqtab.ITS.ASVs)
rownames(taxtab)                          

# Preparing Metadata (for sample_data)
metadata <- data.frame(metadata)
rownames(metadata) <- metadata$sample 
rownames(metadata)

# Confirm that the rownames of seqtab and taxtab are identical (ASV sequences)
all(rownames(seqtab.ITS.ASVs) == rownames(taxtab))
#should return TRUE

# Order and confirm that the colnames of seqtab and rownames of metadata are identical (sample IDs)
metadata <- metadata[order(rownames(metadata)), ]
seqtab.ITS.ASVs <- seqtab.ITS.ASVs[,order(colnames(seqtab.ITS.ASVs))]
all(colnames(seqtab.ITS.ASVs) == rownames(metadata))
#should return TRUE

#' Create phyloseq object for downstream generation of OTU table
ps.object <- phyloseq(otu_table(seqtab.ITS.ASVs, taxa_are_rows = TRUE),
                         tax_table(as.matrix(taxtab)), sample_data(metadata))
str(ps.object)

# Save tables from ps.object
save(ps.object, file = file.path(path.6.phyloseq, "ps.object.RData")) # save phyloseq object as .RData 
write.csv(tax_table(ps.object), file = file.path(path.6.phyloseq, "ps.object.taxonomy.csv"))
write.csv(otu_table(ps.object), file = file.path(path.6.phyloseq, "ps.object.abundances.csv"))

#' create a table containing taxonomy and abundance: ####
#' Extract OTU table (abundance matrix) as dataframe from phyloseq object
OTU.table <- as.data.frame(otu_table(ps.object))
colnames(OTU.table) #sequences
rownames(OTU.table) #sample_IDs
print(paste(Sys.Date(), "-", dim(OTU.table)[2], "samples and", dim(OTU.table)[1], "ASVs"))

taxonomy <- as.data.frame(tax_table(ps.object)) # this has your sequences as rownames and the taxonomic assignments as colnames
print(paste(Sys.Date(), "-", dim(taxonomy)[1], "ASVs and", dim(taxonomy)[2], "Taxonomic Kingdoms"))

# confirm that the OTU Table and taxonomy table have matching rownames (ASV sequences)
all(rownames(taxonomy) == rownames(OTU.table))
#should return TRUE

colnames(taxonomy) #sample_IDs
rownames(taxonomy) #sequences
colnames(OTU.table) #sample_IDs
rownames(OTU.table) #sequences

#merge otutab and taxtab by rownames
taxonomy$ASV <- rownames(taxonomy)
OTU.table$ASV <- rownames(OTU.table)

species.abundance <- merge(OTU.table, taxonomy, by = "ASV")
head(species.abundance)

# save merged Taxonomy Abundance File
save(species.abundance, file = file.path(path.6.phyloseq, "species.abundance.RData")) # save species.abundance object as .RData 
write_xlsx(species.abundance, file.path(path.6.phyloseq, "species.abundance.xlsx"))

