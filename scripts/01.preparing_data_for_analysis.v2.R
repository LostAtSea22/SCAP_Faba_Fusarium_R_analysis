script_name <- "01.preparing_data_for_analysis.r"
print(paste("script name:", script_name))

#' this script imports the data into R
#' 
#' duration: none of these steps take very long

# Finding your data
path.0.project <- file.path(here()) # home directory, makes use of the package, "here"
path.1.data <- file.path(path.0.project, "data", "1.data")
list.files(path.1.data)

# 28 files (20251219 LL)

# Sort ensures forward/reverse reads are in same order
names.1.data.R1 <- sort(list.files(path.1.data, pattern = "_R1_001.fastq.gz")) #sort(...): This function sorts the list of filenames alphabetically
names.1.data.R2 <- sort(list.files(path.1.data, pattern = "_R2_001.fastq.gz")) #pattern = : this function ensures only files that match this pattern are included in the list.

get.sample.name <- function(names) {strsplit(basename(names), "_")[[1]][1]}
sample.names <- unname(sapply(names.1.data.R1, get.sample.name))
sample.names
length(sample.names)
# 15 samples (20250725 LL)

metadata <- read_xlsx(here("data", "metadata.xlsx"))
str(metadata) #str() function stands for "structure" and provides a compact, human-readable description of any R object.
# 14 observations X 5 variables (20250725 LL)
head(metadata)

# expand names to include full path
names.1.data.R1 <- file.path(path.1.data, names.1.data.R1)
names.1.data.R2 <- file.path(path.1.data, names.1.data.R2)
