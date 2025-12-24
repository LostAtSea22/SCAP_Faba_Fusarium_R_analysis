#` # Create Abundance (OTU/ASV) Table using phyloseq
#' **NOTE**: the variarble, `classification.ref.ITS` is assigned in the script `00.loading_packages.v2.R`
#' # set file path
path.6.phyloseq <- file.path(here("data", "6.phyloseq_on_all_merged_reads"))
if(!dir.exists(path.6.phyloseq)) dir.create(path.6.phyloseq)

#' # Import Data

#' ## Import ASVs (amplicon sequence variants)
#' I wanted to assign taxonomy to all reads (not just those that ITSx assigns to oomycetes)
#' therefore I will extract the ASVs from the column headings of the dada2 table (seqtab)
seqtab <- read.csv(file.path(here("data", "4.dada", "2025-08-01SeqTable_chimeras_rm_merged_pairs.csv")))
seqtab <- t(seqtab)
seqtab <- as.data.frame(seqtab)
seqtab$ASV <- rownames(seqtab)
rownames(seqtab) <- NULL
ASVs.df <- seqtab[-1,"ASV"]
ASVs.df <- as.data.frame(ASVs.df)
nrow(ASVs.df)

#' import the Sequence Table which is essentially a matrix where the rownames are your sample names, each column is an ASV and the values are the relative counts or hits for each ASV in the sample
#' note the seqtab contains all ASVs assigned by dada2 with chimeras removed, from mereged pairs
seqtab <- read.csv(file.path(here("data", "4.dada", "2025-08-01SeqTable_chimeras_rm_merged_pairs.csv")))
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

#' # rename seqtab to seqtab.ASVs (artifact of the alternate scripts for processing post ITSx ASVs)
seqtab.ASVs <- seqtab
dim(seqtab.ASVs)
print(paste(Sys.Date(), "-", dim(seqtab.ASVs)[1],"samples and", dim(seqtab.ASVs)[2], "sequence variants from ITS filtering step"))


#' # Remove Negative Control in data 
# metadata <- metadata[-10, ] 
# 20240812 removing the negative control I had in my data (which was row 10) (20240812 LL)
# 20250728 I commented this out because I am worried about contamination in my negative control so I'm leaving it in 

#' # Phyloseq
#' to use phyloseq we need:
#'    - all the sample names to match for the otu table (seqtab.ASVs), the tax table (taxtab) and the sample data (metadata)
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
seqtab.ASVs <- as.matrix(seqtab.ASVs)
seqtab.ASVs <- t(seqtab.ASVs)
str(seqtab.ASVs)

#' remove undetermined reads (reads that illuina was unable to assign to a specific index, including PhiX reads) from data (probably should have done this earlier on in the analysis???)
colnames(seqtab.ASVs)               
seqtab.ASVs <- subset(seqtab.ASVs, select = -c(Undetermine))
colnames(seqtab.ASVs)

#' Preparing ASV Taxonomy Table (for tax_table)
#' the ASVs in taxtab are trimmed versions of the ASVs in seqtab.ASVs (becasue ITSx trims to ITS1 or ITS2)
#' but phyloseq requires the two matrices to have identical rownames therefore we are assigning the seqtab.ASVs ASVs (rownames) to the taxtab rownames
rownames(taxtab) <- rownames(seqtab.ASVs)
rownames(taxtab)                          

# Preparing Metadata (for sample_data)
metadata <- data.frame(metadata)
rownames(metadata) <- metadata$sample 
rownames(metadata)

# Confirm that the rownames of seqtab and taxtab are identical (ASV sequences)
all(rownames(seqtab.ASVs) == rownames(taxtab))
# should return TRUE

# Order and confirm that the colnames of seqtab and rownames of metadata are identical (sample IDs)
metadata <- metadata[order(rownames(metadata)), ]
seqtab.ASVs <- seqtab.ASVs[,order(colnames(seqtab.ASVs))]
all(colnames(seqtab.ASVs) == rownames(metadata))
# should return TRUE

#' Create phyloseq object for downstream generation of OTU table
ps.object <- phyloseq(otu_table(seqtab.ASVs, taxa_are_rows = TRUE),
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
save(species.abundance, file = file.path(path.6.phyloseq, "species.abundance.RData")) # save phyloseq object as .RData 
write_xlsx(species.abundance, file.path(path.6.phyloseq, "species.abundance.xlsx"))
