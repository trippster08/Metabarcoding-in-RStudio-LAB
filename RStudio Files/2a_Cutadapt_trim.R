# 2 RUN CUTADAPT ###############################################################

## Load Libraries ==============================================================
# Load all R packages you may need, if not coming directly from
# "1_Metabarcoding_R_Pipeline_ComputerPrep".

library(tidyverse)
library(ShortRead)

## File Housekeeping ===========================================================

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd(
  "/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME"
)

# Make a list of all the files in your "data/raw" folder.
reads_to_trim <- list.files("data/raw")
head(reads_to_trim)
# Separate files by read direction (R1,R2), and save each
reads_to_trim_F <- reads_to_trim[str_detect(reads_to_trim, "R1_001.fastq.gz")]
reads_to_trim_R <- reads_to_trim[str_detect(reads_to_trim, "R2_001.fastq.gz")]
# Look to ensure that there are the same number of F and R reads
length(reads_to_trim_F)
length(reads_to_trim_R)
# Separate the elements of "reads_to_trim_F" by underscore, and save the first
# element as "sample_names".
sample_names_raw <- sapply(strsplit(basename(reads_to_trim_F), "_"), `[`, 1)
head(sample_names_raw)

# Count the number of reads in each sample.
sequence_counts_raw <- sapply(
  paste0("data/raw/", reads_to_trim_F),
  function(file) {
    fastq_data <- readFastq(file)
    length(fastq_data)
  }
)
# Name these counts with your sample names
names(sequence_counts_raw) <- sample_names_raw
head(sequence_counts_raw)

# Define the path to your primer definition fasta file, if you have more than
# one potential primer to trim. This path will be different for each user.

# At LAB, we use both Nextera and iTru sequencing primers. Currently, our Truseq
# primers include spacers between the sequencing primer and the amplicon primer
# (our Nextera primers do not contain these spacers). Make sure the primer
# definition file includes primers with these spacers attached, otherwise all
# reads will be discared as untrimmed. We have primer definition files for the
# standard COI, 12S MiFish, and 18S_V4 iTru primer (with spacers) and for COI
# and 12S MiFish nextera primers (without spacers). Our primer definition files
# will be downloaded when you download the pipeline. To use a pair of these
# files, just replace "PRIMERF" or "PRIMERR" with the name of the actual primer
# file in the directory "primers" in the paths defined below. If you have
# different primers to remove, open an existing primer file and follow the
# formatting in creating your new primer files.

# If your reads are short, and there is potential for readthrough, you need to
# tell cutadapt to look for primers on the 3' end of each read, as well. These
# primers will be ther reverse complement of the normal primers. They also will
# not be anchored, so the files don't need to include any spacers, and if they
# are not found, the read will still be kept.

# If you know that there will not be any readthrough, you can remove the two
# paths to the RC primers, and two entire lines from the cutadapt arguments:
# following and including "-a" and following and including "-A". Also remove
# "-n 2", because you don't need to run cutadapt twice, since each read will
# only have one primer.

# Again, for the path to the primer files, replace "PRIMERF" or "PRIMERR" with
# the name of the forward and reverse primer file, respectively.

# THE PATHS SHOWN BELOW ARE EXAMPLES ONLY. PLEASE CHANGE PATH TO YOUR PRIMER FILES.
path_to_Fprimers <- "Metabarcoding-in-RStudio-LAB-main/primers/PRIMERF.fas"
path_to_Rprimers <- "Metabarcoding-in-RStudio-LAB-main/primers/PRIMERR.fas"
path_to_FprimersRC <- "Metabarcoding-in-RStudio-LAB-main/primers/PRIMERF_RC.fas"
path_to_RprimersRC <- "Metabarcoding-in-RStudio-LAB-main/primers/PRIMERR_RC.fas"
# Make sure it worked
path_to_Fprimers

## Run Cutadapt ================================================================

# Save the path to the cutadapt executable file. Your path will be different.
cutadapt_binary <- "/Users/macdonaldk/mambaforge/envs/cutadapt/bin/cutadapt"

# The following for loop runs cutadapt on paired samples, one pair at a time.

# If you are not using a primer definition fasta file, and are only attempting
# to trim a single primer from R1 and a single primer from R2, replace
# "paste0("file:",path_to_Fprimers)" with the primer sequence. Include a "^"
#appended to the 5' end. Do the same for both "-g" and "-G". For example:
# "-g ^FORWARDPRIMERSEQUENCE" in quotations, followed by a comma, and
# "-G ^REVERSEPRIMERSEQUENCE" in quotations, followed by a comma.
# fmt: skip

### Run cutadapt no 3' trimming ------------------------------------------------
# Run this if you have no read-through in sequences. In other words, you should
# not find any primers on the 3' end of the sequence
for (i in seq_along(sample_names_raw)) {
  system2(
    cutadapt_binary,
    args = c(
      "-e 0.2 --discard-untrimmed --minimum-length 30 --cores=0",
      "-g", paste0("file:", path_to_Fprimers),
      "-G", paste0("file:", path_to_Rprimers),
      "-o", paste0(
        "data/working/trimmed_sequences/",
        sample_names_raw[i],
        "_trimmed_R1.fastq.gz"
      ), "-p", paste0(
        "data/working/trimmed_sequences/",
        sample_names_raw[i],
        "_trimmed_R2.fastq.gz"
      ),
      paste0("data/raw/", reads_to_trim_F[i]),
      paste0("data/raw/", reads_to_trim_R[i])
    )
  )
}

### Run cutadapt WITH 3' trimming ----------------------------------------------
# Run this if you have read-through in sequences. In other words, you may have
# primers on the 3' end of reads
for (i in seq_along(sample_names_raw)) {
  system2(
    cutadapt_binary,
    args = c(
      "-e 0.2 --discard-untrimmed --minimum-length 30 -n 2 -O 3 --cores=0",
      "-g",
      paste0("file:", path_to_Fprimers),
      "-a",
      paste0("file:", path_to_RprimersRC),
      "-G",
      paste0("file:", path_to_Rprimers),
      "-A",
      paste0("file:", path_to_FprimersRC),
      "-o",
      paste0(
        "data/working/trimmed_sequences/",
        sample_names_raw[i],
        "_trimmed_R1.fastq.gz"
      ),
      "-p",
      paste0(
        "data/working/trimmed_sequences/",
        sample_names_raw[i],
        "_trimmed_R2.fastq.gz"
      ),
      paste0("data/raw/", reads_to_trim_F[i]),
      paste0("data/raw/", reads_to_trim_R[i])
    )
  )
}
# Save all the objects we've created so far so we don't have to create these
# again down the road if we leave this project
save(
  reads_to_trim,
  reads_to_trim_F,
  reads_to_trim_R,
  sample_names_raw,
  sequence_counts_raw,
  path_to_Fprimers,
  path_to_Rprimers,
  path_to_FprimersRC,
  path_to_RprimersRC,
  file = "data/working/1_trim.RData"
)

# We are including our default parameters for cutadapt. You can change these
# parameters if you have prefer others.

# -e 0.2 allows an error rate of 0.2 (20% of primer basepairs can me wrong)

# --minimum-length 30 removes all reads that are not at least 30 bp. However,
# as currently implemented in cutadapt, this does not always work correctly,
# and sometimes it removes the sequence of the reads, but not the name, leaving
# empty reads. We deal with this later in the pipeline.

# -O 3 This is the minimum number of basepairs that the primer must overlap the
# read to be counted and removed. Three is the default. This only affects the
# 3' primer. The 5' primer is anchored, which means it must be found in its
# entirety, or the read is removed

# -n 2 This is the number of times to run cutadapt on each sample. You may need
# to run cutadapt twice if you have readthrough and need to remove primers from
# the 3' end as well as the 5' end. Cutadapt will only remove one primer from
# a sample each time it's run.

# --cores=0 tells cutadapt how many cores to use while trimming. 0 sets cutadapt
# to automatically detect the number of cores.

# -g and -G are the paths to the 5' primers with spacers
# -a and -A are the paths to the 3' RC primers. If you are certain there is no
# possibility of read-through, you can omit -a and -A

# -o is the output for trimmed R1 reads
# -p is the output for the trimmed paired R2 reads

# The final line contains the two paired input files
