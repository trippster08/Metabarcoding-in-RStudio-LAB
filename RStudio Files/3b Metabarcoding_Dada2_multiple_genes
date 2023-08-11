# 4 DADA2 ######################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous step. 
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

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd("/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")

# Set a path to the directory with the cutadapt-trimmed reads for each gene. 
# I've used "gene1" and "gene2" here, find and replace with the name of the
# each gene-specific directory used in by cutadapt in the previous step.
path.gene1 <- "data/working/trimmed_sequences/gene1"
path.gene2 <- "data/working/trimmed_sequences/gene2"
# This lists the files inside the selected folder.
list.files(path.gene1)
list.files(path.gene2)

# We will now run the rest of this section in multiple times, once for each
# gene present in the Illumina run. Again, replace each "gene1", "gene2",
# "gene3", etc with your specific gene name.

## Gene1 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs.gene1 <- sort(list.files(path.gene1, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs.gene1 <- sort(list.files(path.gene1, pattern="_R2.fastq.gz", full.names = TRUE))

# Make sure you have the same number of forward reads as reverse reads
length(fnFs.gene1)
length(fnRs.gene1)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(fnFs.gene1)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.gene1.exists <- fnFs.gene1[file.size(fnFs.gene1) > 50 & file.size(fnRs.gene1) > 50]
length(fnFs.gene1.exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.gene1.exists <- fnRs.gene1[file.size(fnFs.gene1) > 50 & file.size(fnRs.gene1) > 50]
length(fnRs.gene1.exists)
file.size(fnFs.gene1.exists)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs.gene1 <- fnFs.gene1.exists
fnRs.gene1 <- fnRs.gene1.exists
length(fnFs.gene1)
length(fnRs.gene1)
file.size(fnFs.gene1)

# Update your samples names
sample.names.gene1 <- sapply(strsplit(basename(fnFs.gene1), "_"), `[`, 1)
length(sample.names.gene1)
head(sample.names.gene1)

## Filter and Trim =============================================================

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets (to aggregate all samples, replace N
# with the number of samples, or with length(fnFs)). For example, "fnFs[1:2]"
# will result in two plots, one for the first sample and one for the second.
# "fnFs[1:17]" will result in 17 plots, one for each sample.  Using
# "aggregate = TRUE" will combine any samples called (for example, "fnFS[1:17]"
# aggregates sample 1 through 17) into a single plot. This results in the same
# the quality plots as Qiime2.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.
qualplotF.gene1 <- plotQualityProfile(
  fnFs.gene1[1:N],
  aggregate = TRUE
)
qualplotF.gene1

# Here we modify our quality plot to better visualize where the quality cut-off
# should be. "Scale_x_continuous" is the minimum and maximum x-axis values to be
# shown.  We use "breaks=seq(a,b,c)", to indicate the first axis tick "a", last
# tick "b", and frequency of ticks "c". The example shown results in a plot that
# starts at 190 bp, ends at 220 bp, and has axis ticks every 2 bp.
qualplotF.gene1 + scale_x_continuous(limits=c(190,220), breaks=seq(190,220,2))

# I also sometimes look at all (or most of) the individual quality plots, to see
# if there are any troublesome samples. Negative controls often have poor
# quality, but you may want to rerun any samples that are substantially worse
# than the rest. 
plotQualityProfile(
  fnFs.gene1[1:N],
  aggregate = FALSE
)

# Examine the reverse reads as you did the forward.
qualplotR.gene1 <- plotQualityProfile(
  fnRs.gene1[1:N],
  aggregate = TRUE
)
qualplotR.gene1
qualplotR.gene1 + scale_x_continuous(limits=c(150,200), breaks=seq(150,200,2))

plotQualityProfile(
  fnRs.gene1[1:N],
  aggregate = FALSE
)

# This creates files for the reads that will be quality filtered with dada2
  # in the next step.  
filtFs.gene1 <- file.path(path.gene1, "filtered", paste0(sample.names.gene1, "_F_filt.fastq.gz"))
filtRs.gene1 <- file.path(path.gene1, "filtered", paste0(sample.names.gene1, "_R_filt.fastq.gz"))

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtFs and filtRs goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtFs.gene1) <- sample.names.gene1
names(filtRs.gene1) <- sample.names.gene1

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(i,j)" is how you tell Dada2 where to truncate all forward (i) and
# reverse (j) reads. Using "0" means reads will not be truncated. 
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 6,6) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

out.gene1 <- filterAndTrim(
  fnFs.gene1,
  filtFs.gene1,
  fnRs.gene1,
  filtRs.gene1,
  truncLen=c(0,0),
  maxN=0,
  maxEE=c(4,4),
  rm.phix=TRUE,
  truncQ=2,
  compress=TRUE,
  multithread=TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out.gene1
head(out.gene1)

# After filtering, if there are any samples that have no remaining reads
# (i.e. reads.out = 0), you will get the following error running learnErrors:
# "Error in derepFastq(fls[[i]], qualityType = qualityType) : Not all provided
# files exist." That is because while these empty sample names still exist in
# filtFs and filtRs, there is no data connected to these names, and dada2
# doesn't like that.

# This step changes filtFs and filtRs to only contain the names of samples with
# reads. Do this only if there are samples in "out" with zero reads.

# You will notice that the number of items in filtFs is now the number of
# samples with reads (i.e. the description for filtFs and filtRs goes from
# "Named chr [1:N]" to "Named chr [1:N-(# of empty samples)]).
exists.gene1 <- file.exists(filtFs.gene1) & file.exists(filtRs.gene1)
exists.gene1
filtFs.gene1 <- filtFs.gene1[exists.gene1]
filtRs.gene1 <- filtRs.gene1[exists.gene1]

# I sometimes look at the quality of the filtered reads. Sometimes my filtering
# parameters don't get rid of all the poor quality at the 3' end, but if quality
# is too poor, they won't merge, so it doesn't seem to hurt read passage.
plotQualityProfile(filtFs.gene1[1:N], aggregate = TRUE)
plotQualityProfile(filtRs.gene1[1:N], aggregate = TRUE)

## Estimating Error Rates and Denoising ========================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errF.gene1 <- learnErrors(
  filtFs.gene1,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errR.gene1 <- learnErrors(
  filtRs.gene1, 
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal defintion
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality.
plotErrors(errF.gene1, nominalQ=TRUE)
plotErrors(errR.gene1, nominalQ=TRUE)

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtFs), and "err =" which is the error file from
# learnErrors (effF).
dadaFs.gene1 <- dada(
  filtFs.gene1, 
  err = errF.gene1, 
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

dadaRs.gene1 <- dada(
  filtRs.gene1, 
  err=errR.gene1, 
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
dadaFs.gene1[[1]]
dadaRs.gene1[[1]]

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (dadaFs), then the forward filtered and truncated reads (filtFs), then the
# same for the reverse reads (dadaRs and filtRs).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged.gene1 <- mergePairs(
  dadaFs.gene1, 
  filtFs.gene1, 
  dadaRs.gene1, 
  filtRs.gene1, 
  minOverlap = 12,
  maxMismatch = 0,
  verbose=TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged.gene1[[1]])
head(merged.gene1[[6]])

## Create and Trim Sequence-Table ==============================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab.gene1 <- makeSequenceTable(merged.gene1)
# This describes the dimensions of the table just made
dim(seqtab.gene1)

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
table(nchar(getSequences(seqtab.gene1)))

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding genes (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab313 <- seqtab[,nchar(colnames(seqtab)) %in% 298:322]
dim(seqtab313)
table(nchar(getSequences(seqtab313)))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab.nochim.gene1 <- removeBimeraDenovo(
  seqtab.gene1,
  method="consensus",
  multithread=TRUE,
  verbose=TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab.nochim.gene1)

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2. I've added a column to the typical
# tutorial version of this that gives us the percentage of reads that made it
# through the process. This is a good quick way to see if something is wrong
# (i.e. only a small proportion make it through).
getN <- function(x) sum(getUniques(x))
track.gene1 <- cbind(
  out.gene1, 
  sapply(dadaFs.gene1, getN), 
  sapply(dadaRs.gene1, getN), 
  sapply(merged.gene1, getN), 
  rowSums(seqtab.nochim.gene1), 
  100*(rowSums(seqtab.nochim.gene1) / out.gene1[,1]))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track.gene1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "%kept")
rownames(track.gene1) <- sample.names.gene1

# Look at the results for the first 6 samples, or all, depending upon the number
# of samples.
head(track.gene1)
track.gene1
# !!!!!!!Note, if the initial Filter and Trim step left any samples with 0
# reads, and you had to use the file.exists command, it will cause a problem
# tracking reads. I'm working on it.


## Export Sequence-Table =======================================================
# This exports a sequence-table: columns of ASV's, rows of samples, and
# values = number of reads. This is the only export you need for downstream
# analyses. You can do anything you want in this pipeline with this table. This
# table can also be easily merged with other tables from the same project but
# from different runs before downstream analyses.

# If you have mulitple Miseqruns for the same project that will need to be
# combined for further analyses, you may want to name this file
# "PROJECTNAME_MISEQRUN_sequence-table.tsv" to differentiate different runs.
# In "5 Metabarcoding_R_Pipeline_RStudio_ImportCombine" we'll show how to
# combine data from separate runs for analyses.
write.table(
  seqtab.nochim.gene1,
  file="data/results/PROJECTNAME_gene1_sequence-table.tsv",
  quote = FALSE,
  sep="\t",
  row.names = TRUE,
  col.names = NA
)


## Gene2 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs.gene2 <- sort(list.files(path.gene2, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs.gene2 <- sort(list.files(path.gene2, pattern="_R2.fastq.gz", full.names = TRUE))

# Make sure you have the same number of forward reads as reverse reads
length(fnFs.gene2)
length(fnRs.gene2)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(fnFs.gene2)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.gene2.exists <- fnFs.gene2[file.size(fnFs.gene2) > 50 & file.size(fnRs.gene2) > 50]
length(fnFs.gene2.exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.gene2.exists <- fnRs.gene2[file.size(fnFs.gene2) > 50 & file.size(fnRs.gene2) > 50]
length(fnRs.gene2.exists)
file.size(fnFs.gene2.exists)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs.gene2 <- fnFs.gene2.exists
fnRs.gene2 <- fnRs.gene2.exists
length(fnFs.gene2)
length(fnRs.gene2)
file.size(fnFs.gene2)

# Update your samples names
sample.names.gene2 <- sapply(strsplit(basename(fnFs.gene2), "_"), `[`, 1)
length(sample.names.gene2)
head(sample.names.gene2)

## Filter and Trim =============================================================

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets (to aggregate all samples, replace N
# with the number of samples, or with length(fnFs)). For example, "fnFs[1:2]"
# will result in two plots, one for the first sample and one for the second.
# "fnFs[1:17]" will result in 17 plots, one for each sample.  Using
# "aggregate = TRUE" will combine any samples called (for example, "fnFS[1:17]"
# aggregates sample 1 through 17) into a single plot. This results in the same
# the quality plots as Qiime2.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.
qualplotF.gene2 <- plotQualityProfile(
  fnFs.gene2[1:N],
  aggregate = TRUE
)
qualplotF.gene2

# Here we modify our quality plot to better visualize where the quality cut-off
# should be. "Scale_x_continuous" is the minimum and maximum x-axis values to be
# shown.  We use "breaks=seq(a,b,c)", to indicate the first axis tick "a", last
# tick "b", and frequency of ticks "c". The example shown results in a plot that
# starts at 190 bp, ends at 220 bp, and has axis ticks every 2 bp.
qualplotF.gene2 + scale_x_continuous(limits=c(190,220), breaks=seq(190,220,2))

# I also sometimes look at all (or most of) the individual quality plots, to see
# if there are any troublesome samples. Negative controls often have poor
# quality, but you may want to rerun any samples that are substantially worse
# than the rest. 
plotQualityProfile(
  fnFs.gene2[1:N],
  aggregate = FALSE
)

# Examine the reverse reads as you did the forward.
qualplotR.gene2 <- plotQualityProfile(
  fnRs.gene2[1:N],
  aggregate = TRUE
)
qualplotR.gene2
qualplotR.gene2 + scale_x_continuous(limits=c(150,200), breaks=seq(150,200,2))

plotQualityProfile(
  fnRs.gene2[1:N],
  aggregate = FALSE
)

# This creates files for the reads that will be quality filtered with dada2
  # in the next step.  
filtFs.gene2 <- file.path(path.gene2, "filtered", paste0(sample.names.gene2, "_F_filt.fastq.gz"))
filtRs.gene2 <- file.path(path.gene2, "filtered", paste0(sample.names.gene2, "_R_filt.fastq.gz"))

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtFs and filtRs goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtFs.gene2) <- sample.names.gene2
names(filtRs.gene2) <- sample.names.gene2

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(i,j)" is how you tell Dada2 where to truncate all forward (i) and
# reverse (j) reads. Using "0" means reads will not be truncated. 
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 6,6) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

out.gene2 <- filterAndTrim(
  fnFs.gene2,
  filtFs.gene2,
  fnRs.gene2,
  filtRs.gene2,
  truncLen=c(0,0),
  maxN=0,
  maxEE=c(4,4),
  rm.phix=TRUE,
  truncQ=2,
  compress=TRUE,
  multithread=TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out.gene2
head(out.gene2)

# After filtering, if there are any samples that have no remaining reads
# (i.e. reads.out = 0), you will get the following error running learnErrors:
# "Error in derepFastq(fls[[i]], qualityType = qualityType) : Not all provided
# files exist." That is because while these empty sample names still exist in
# filtFs and filtRs, there is no data connected to these names, and dada2
# doesn't like that.

# This step changes filtFs and filtRs to only contain the names of samples with
# reads. Do this only if there are samples in "out" with zero reads.

# You will notice that the number of items in filtFs is now the number of
# samples with reads (i.e. the description for filtFs and filtRs goes from
# "Named chr [1:N]" to "Named chr [1:N-(# of empty samples)]).
exists.gene2 <- file.exists(filtFs.gene2) & file.exists(filtRs.gene2)
exists.gene2
filtFs.gene2 <- filtFs.gene2[exists.gene2]
filtRs.gene2 <- filtRs.gene2[exists.gene2]

# I sometimes look at the quality of the filtered reads. Sometimes my filtering
# parameters don't get rid of all the poor quality at the 3' end, but if quality
# is too poor, they won't merge, so it doesn't seem to hurt read passage.
plotQualityProfile(filtFs.gene2[1:N], aggregate = TRUE)
plotQualityProfile(filtRs.gene2[1:N], aggregate = TRUE)

## Estimating Error Rates and Denoising ========================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errF.gene2 <- learnErrors(
  filtFs.gene2,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errR.gene2 <- learnErrors(
  filtRs.gene2, 
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal defintion
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality.
plotErrors(errF.gene2, nominalQ=TRUE)
plotErrors(errR.gene2, nominalQ=TRUE)

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtFs), and "err =" which is the error file from
# learnErrors (effF).
dadaFs.gene2 <- dada(
  filtFs.gene2, 
  err = errF.gene2, 
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

dadaRs.gene2 <- dada(
  filtRs.gene2, 
  err=errR.gene2, 
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
dadaFs.gene2[[1]]
dadaRs.gene2[[1]]

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (dadaFs), then the forward filtered and truncated reads (filtFs), then the
# same for the reverse reads (dadaRs and filtRs).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged.gene2 <- mergePairs(
  dadaFs.gene2, 
  filtFs.gene2, 
  dadaRs.gene2, 
  filtRs.gene2, 
  minOverlap = 12,
  maxMismatch = 0,
  verbose=TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged.gene2[[1]])
head(merged.gene2[[6]])

## Create and Trim Sequence-Table ==============================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab.gene2 <- makeSequenceTable(merged.gene2)
# This describes the dimensions of the table just made
dim(seqtab.gene2)

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
table(nchar(getSequences(seqtab.gene2)))

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding genes (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab313 <- seqtab[,nchar(colnames(seqtab)) %in% 298:322]
dim(seqtab313)
table(nchar(getSequences(seqtab313)))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab.nochim.gene2 <- removeBimeraDenovo(
  seqtab.gene2,
  method="consensus",
  multithread=TRUE,
  verbose=TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab.nochim.gene2)

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2. I've added a column to the typical
# tutorial version of this that gives us the percentage of reads that made it
# through the process. This is a good quick way to see if something is wrong
# (i.e. only a small proportion make it through).
getN <- function(x) sum(getUniques(x))
track.gene2 <- cbind(
  out.gene2, 
  sapply(dadaFs.gene2, getN), 
  sapply(dadaRs.gene2, getN), 
  sapply(merged.gene2, getN), 
  rowSums(seqtab.nochim.gene2), 
  100*(rowSums(seqtab.nochim.gene2) / out.gene2[,1]))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track.gene2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "%kept")
rownames(track.gene2) <- sample.names.gene2

# Look at the results for the first 6 samples, or all, depending upon the number
# of samples.
head(track.gene2)
track.gene2
# !!!!!!!Note, if the initial Filter and Trim step left any samples with 0
# reads, and you had to use the file.exists command, it will cause a problem
# tracking reads. I'm working on it.


## Export Sequence-Table =======================================================
# This exports a sequence-table: columns of ASV's, rows of samples, and
# values = number of reads. This is the only export you need for downstream
# analyses. You can do anything you want in this pipeline with this table. This
# table can also be easily merged with other tables from the same project but
# from different runs before downstream analyses.

# If you have mulitple Miseqruns for the same project that will need to be
# combined for further analyses, you may want to name this file
# "PROJECTNAME_MISEQRUN_sequence-table.tsv" to differentiate different runs.
# In "5 Metabarcoding_R_Pipeline_RStudio_ImportCombine" we'll show how to
# combine data from separate runs for analyses.
write.table(
  seqtab.nochim.gene2,
  file="data/results/PROJECTNAME_gene2_sequence-table.tsv",
  quote = FALSE,
  sep="\t",
  row.names = TRUE,
  col.names = NA
)