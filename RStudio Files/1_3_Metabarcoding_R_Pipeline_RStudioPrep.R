## Open RStudio ================================================================
### Istall and Load R Libraries ------------------------------------------------
# These are the libraries that will be used for this pipeline. Not all will be
# used for each persons particular project, but it does not hurt to have them
# all installed, loaded, and ready to be used when needed.
# Install  BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
# Install Dada2
BiocManager::install("dada2", version = "3.15")
# Install Phyloseq
BiocManager::install("phyloseq")
# Install DECIPHER
BiocManager::install("DECIPHER")

# Install and other libraries you may need (or install through
# "Install Packages" window). Libraries will only need to be installed once.
install.packages("digest")
install.packages("tidyverse")
install.packages("seqinr")
install.packages("ape")
install.packages("ade4")
install.packages("filesstrings")

# Load all the libraries that will be needed for this pipeline. These will have
# to be reloaded every time you restart RStudio
library(dada2)
library(digest)
library(phyloseq)
library(tidyverse)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)
library(filestrings)

## File Housekeeping ===========================================================
# Set up your working directory. Substitute your own path for the one below.
# Note: If you have spaces or special characters in the path to your working
# directory, you don't need a character escape (i.e. no \ preceding spaces or
# special characters).
setwd ("/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")

# Find all the read files in the project directory, save their paths, and
# confirm. Basespace saves the reads in sample-specific folders, using
# "recursive = TRUE" allows us to find all read files in the working directory
raw_reads <- list.files(pattern = ".fastq.gz", recursive = TRUE)
head(raw_reads)

# Create all the subdirectories we will use
dir.create("data/raw", recursive = TRUE)
dir.create("data/working/trimmed_sequences", recursive=TRUE)
dir.create("data/results")

# Copy the read files to the "data/raw" directory, and confirm that they are
# there.
file.copy(raw_reads, "data/raw", recursive=TRUE)
head(list.files("data/raw"))