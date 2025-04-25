## Open RStudio ================================================================
### Istall and Load R Libraries ------------------------------------------------
# These are the libraries that will be used for this pipeline. Not all will be
# used for each persons particular project, but it does not hurt to have them
# all installed, loaded, and ready to be used when needed.
# Install  BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Install all of the libraries needed through BiocManager.
BiocManager::install("dada2", ask = FALSE)
BiocManager::install("phyloseq", ask = FALSE)
BiocManager::install("msa", ask = FALSE)
BiocManager::install("DECIPHER", ask = FALSE)
BiocManager::install("rBLAST", ask = FALSE)
BiocManager::install("ShortRead", ask = FALSE)

# Install and other libraries you may need (or install through
# "Install Packages" window). Libraries will only need to be installed once.
# If you get a message saying some packages have more recent versions available,
# and asking if you want to update them, chose "1: ALL".
install.packages("digest")
install.packages("tidyverse")
install.packages("seqinr")
install.packages("ape")
install.packages("vegan")
install.packages("patchwork")
install.packages("remotes")
install.packages("R.utils")
install.packages("phylotools")
install.packages("data.table")
remotes::install_github("ropensci/bold", upgrade = TRUE)
remotes::install_github("ropensci/taxize", upgrade = TRUE)
remotes::install_github("fkeck/refdb", upgrade = TRUE)
remotes::install_github("tobiasgf/lulu", upgrade = TRUE)
remotes::install_github("boldsystems-central/BOLDconnectR", upgrade = TRUE)
install.packages("rMSA", repos = "https://mhahsler.r-universe.dev")

## File Housekeeping ===========================================================
# Set up your working directory. If you created your new project in the
# directory you want as your working directory, you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
# Note: If you have spaces or special characters in the path to your working
# directory, you don't need a character escape (i.e. no \ preceding spaces or
# special characters).
setwd(
  "/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME"
)
# Save project name as an object
project_name <- basename(getwd())

# Create all the subdirectories we will use
# Define the directory names
dir_names <- c(
  "data/raw",
  "data/working/trimmed_sequences",
  "data/results",
  "ref"
)
# Create the directories using sapply
sapply(dir_names, dir.create, recursive = TRUE)

# Find all the read files in the project directory, save their paths, and
# confirm. Basespace saves the reads in sample-specific folders, using
# "recursive = TRUE" allows us to find all read files in the working directory
raw_reads <- list.files(pattern = ".fastq.gz", recursive = TRUE)
head(raw_reads)

# Copy the read files to the "data/raw" directory, and confirm that they are
# there.
file.copy(raw_reads, "data/raw", recursive = TRUE)
head(list.files("data/raw"))
