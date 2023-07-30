# 2 RUN CUTADAPT ###############################################################

## Load Libraries ==============================================================
# Load all R packages you may need, if not coming directly from
# "1_Metabarcoding_R_Pipeline_ComputerPrep".

library(dada2)
library(digest)
library(phyloseq)
library(tidyverse)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)
library(filesstrings)

## File Housekeeping ===========================================================

# Set your working directory to you project directory, again, if not coming
# directly from "1_Metabarcoding_R_Pipeline_ComputerPrep".
setwd ("/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")


# Make a list of all the files in your "data/raw" folder.
reads.to.trim <- list.files("data/raw")
head(reads.to.trim)
# Separate files by read direction (R1,R2), and save each
reads.to.trim.F <- reads.to.trim[str_detect(reads.to.trim, "R1_001.fastq.gz")]
reads.to.trim.R <- reads.to.trim[str_detect(reads.to.trim, "R2_001.fastq.gz")]

# Separate the elements of "reads.to.trim.F" by underscore, and save the first
# element as "sample.names".
sample.names <- sapply(strsplit(basename(reads.to.trim.F), "_"), `[`, 1)
head (sample.names)

# Define the path to your primer definition fasta file, if you have more than
# one potential primer to trim. This path will be different for each user.

# At LAB, we use both Nextera and iTru sequencing primers. Currently, our Truseq
# primers include spacers between the sequencing primer and the amplicon primer
# (our Nextera primers do not contain these spacers). Make sure the primer
# definition file includes primers with these spacers attached, otherwise all
# reads will be discared as untrimmed. We have primer definition files for the
# standard COI, 12S MiFish, and 18S_V4 iTru primer (with spacers) and for COI and 12S
# MiFish nextera primers (without spacers).

# If your reads are short, and there is potential for readthrough, you need to
# tell cutadapt to look for primers on the 3' end of each read, as well. These
# primers will be ther reverse complement of the normal primers. They also will
# not be anchored, so the files don't need to include any spacers, and if they
# are not found, the read will still be kept. If you know that there will not
# be any readthrough, you don't have to include the two paths to the RC primers.
# For the path to the primer files, replace "PRIMERF" or "PRIMERR" with the name
# of the forward and reverse primer file, respectively.

path.to.Fprimers <- "'/Users/macdonaldk/Dropbox (Smithsonian)/Metabarcoding/Metabarcoding_Pipeline_in_RStudio/primer_files/PRIMERF.fas'"
path.to.Rprimers <- "'/Users/macdonaldk/Dropbox (Smithsonian)/Metabarcoding/Metabarcoding_Pipeline_in_RStudio/primer_files/PRIMERR.fas'"
path.to.FprimersRC <- "'/Users/macdonaldk/Dropbox (Smithsonian)/Metabarcoding/Metabarcoding_Pipeline_in_RStudio/primer_files/PRIMERF_RC.fas'"
path.to.RprimersRC <- "'/Users/macdonaldk/Dropbox (Smithsonian)/Metabarcoding/Metabarcoding_Pipeline_in_RStudio/primer_files/PRIMERR_RC.fas'"

## Run Cutadapt ================================================================

# Save the path to the cutadapt executable file. Your path will be different.
cutadapt_binary <- "/Users/macdonaldk/miniconda3/envs/cutadapt/bin/cutadapt"

# The following for loop runs cutadapt on paired samples, one pair at a time.

# If you are not using a primer definition fasta file, and are only attempting
# to trim a single primer from R1 and a single primer from R2, replace
# "paste0("file:",path.to.Fprimers)" with the primer sequence. Include a "^"
#appended to the 5' end. Do the same for both "-g" and "-G". For example:
# "-g ^FORWARDPRIMERSEQUENCE" in quotations, followed by a comma, and
# "-G ^REVERSEPRIMERSEQUENCE" in quotations, followed by a comma.


for (i in seq_along(sample.names)) {
  system2(
    cutadapt_binary, args = c(
      "-e 0.2 --discard-untrimmed --minimum-length 30 -n 3 --cores=0",
      "-g", paste0("file:",path.to.Fprimers),
      "-a", paste0("file:",path.to.RprimersRC),
      "-G", paste0("file:",path.to.Rprimers),
      "-A", paste0("file:",path.to.FprimersRC),
      "-o", paste0("data/working/trimmed_sequences/",sample.names[i],"_trimmed_R1.fastq.gz"),
      "-p", paste0("data/working/trimmed_sequences/",sample.names[i],"_trimmed_R2.fastq.gz"),
      paste0("data/raw/",reads.to.trim.F[i]), paste0("data/raw/",reads.to.trim.R[i])
      )
    )
}

# We are including our default parameters for cutadapt. You can change these
# parameters if you have prefer others. 

# -e 0.2 allows an error rate of 0.2 (20% of primer basepairs can me wrong)

# --minimum-length 30 removes all reads that are not at least 30 bp. However,
# as currently implemented in cutadapt, this does not always work correctly,
# and sometimes it removes the sequence of the reads, but not the name, leaving
# empty reads. We deal with this later in the pipeline.

# -n 3 This is the minimum number of basepairs that the primer must overlap the
# read to be counted and removed. Three is the default. This only affects the 
# 3' primer. The 5' primer is anchored, which means it must be found in its
# entirety, or the read is removed

# -N 0 tells cutadapt how many cores to use while trimming. 0 sets cutadapt to
# automatically detect the number of cores.

# -g and -G are the paths to the 5' primers with spacers
# -a and -A are the paths to the 3' RC primers. If you are certain there is no
# possibility of read-through, you can omit -a and -A

# -o is the output for trimmed R1 reads
# -p is the output for the trimmed paired R2 reads

# The final line contains the two paired input files
