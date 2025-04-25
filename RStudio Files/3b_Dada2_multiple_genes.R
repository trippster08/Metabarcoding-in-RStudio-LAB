# 4 DADA2 ######################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need.
library(dada2)
library(digest)
library(tidyverse)
library(seqinr)
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

# Set a path to the directory with the cutadapt-trimmed reads for each gene.
# I've used "GENE1" and "GENE2" here, find and replace with the name of the
# each gene-specific directory used in by cutadapt in the previous step.
path_to_trimmed_GENE1 <- "data/working/trimmed_sequences/GENE1"
path_to_trimmed_GENE2 <- "data/working/trimmed_sequences/GENE2"
# This lists the files inside the selected folder.
list.files(path_to_trimmed_GENE1)
list.files(path_to_trimmed_GENE2)

# We will now run the rest of this section in multiple times, once for each
# gene present in the Illumina run. Again, replace each "GENE1", "GENE2",
# "GENE3", etc with your specific gene name.

## Gene1 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# trimmed_F) and the other for reverse reads (R2, called trimmed_R).
trimmed_F_GENE1 <- sort(list.files(
  path_GENE1,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))
trimmed_R_GENE1 <- sort(list.files(
  path_GENE1,
  pattern = "_R2.fastq.gz",
  full.names = TRUE
))

# Make sure you have the same number of forward reads as reverse reads
length(trimmed_F_GENE1)
length(trimmed_R_GENE1)

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed_GENE1 <- sapply(
  strsplit(basename(trimmed_F_GENE1), "_"),
  `[`,
  1
)
head(sample_names_trimmed_GENE1)
length(sample_names_trimmed_GENE1)

# Now count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed_GENE1 <- sapply(trimmed_F_GENE1, function(file) {
  fastq_data_GENE1 <- readFastq(file)
  length(fastq_data_GENE1)
})
names(sequence_counts_trimmed_GENE1) <- sample_names_trimmed_GENE1
print(sequence_counts_trimmed_GENE1)
# Look at this
sequence_counts_trimmed_GENE1

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_F_GENE1_exists <- trimmed_F_GENE1[
  file.size(trimmed_F_GENE1) > 50 & file.size(trimmed_R_GENE1) > 50
]
length(trimmed_F_GENE1_exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_R_GENE1_exists <- trimmed_R_GENE1[
  file.size(trimmed_F_GENE1) > 50 & file.size(trimmed_R_GENE1) > 50
]
length(trimmed_R_GENE1_exists)
file.size(trimmed_F_GENE1_exists)

# Redefine trimmed_F and trimmed_R as only the existing read files, and check
trimmed_F_GENE1 <- trimmed_F_GENE1_exists
trimmed_R_GENE1 <- trimmed_R_GENE1_exists
length(trimmed_F_GENE1)
length(trimmed_R_GENE1)
file.size(trimmed_F_GENE1)

# Update your samples names
sample_names_GENE1 <- sapply(strsplit(basename(trimmed_F_GENE1), "_"), `[`, 1)
length(sample_names_GENE1)
head(sample_names_GENE1)

## Filter and Trim =============================================================

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets (to aggregate all samples, replace N
# with the number of samples, or with length(trimmed_F)). For example,
# "trimmed_F[1:2]" will result in two plots, one for the first sample and one
# for the second. "trimmed_F[1:17]" will result in 17 plots, one for each
# sample.  Using "aggregate = TRUE" will combine any samples called.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.
qualplotF_GENE1 <- plotQualityProfile(
  trimmed_F_GENE1[1:length(sample_names_trimmed)],
  aggregate = TRUE
)
qualplotF_GENE1

# Export this plot as a pdf.
ggsave(
  "data/results/qualplotR_GENE1.pdf",
  plot = qualplotF_GENE1,
  width = 9,
  height = 9
)

# Here we modify our quality plot to better visualize where the quality cut-off
# should be. "Scale_x_continuous" is the minimum and maximum x-axis values to be
# shown.  We use "breaks=seq(a,b,c)", to indicate the first axis tick "a", last
# tick "b", and frequency of ticks "c". The example shown results in a plot that
# starts at 190 bp, ends at 220 bp, and has axis ticks every 2 bp.
qualplotF_GENE1 +
  scale_x_continuous(limits = c(250, 290), breaks = seq(250, 290, 5))

# Examine the reverse reads as you did the forward.
qualplotR_GENE1 <- plotQualityProfile(
  trimmed_R_GENE1[1:length(trimmed_R)], # nolint: seq_linter.
  aggregate = TRUE
)
qualplotR_GENE1
# Export this plot as a pdf.
ggsave(
  "data/results/qualplotR_GENE1.pdf",
  plot = qualplotR,
  width = 9,
  height = 9
)

qualplotR_GENE1 +
  scale_x_continuous(limits = c(250, 290), breaks = seq(250, 290, 5))

# Save all the objects created to this point in this section
save(
  path_to_trimmed_GENE1,
  trimmed_F_GENE1,
  trimmed_R_GENE1,
  sequence_counts_trimmed_GENE1,
  sample_names_trimmed_GENE1,
  quality_plot_F_GENE1,
  quality_plot_R_GENE1,
  file = "data/working/2_qual_GENE1.Rdata"
)


# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtFs_GENE1 <- file.path(
  path_GENE1,
  "filtered",
  paste0(sample_names_GENE1, "_F_filt.fastq.gz")
)
filtRs_GENE1 <- file.path(
  path_GENE1,
  "filtered",
  paste0(sample_names_GENE1, "_R_filt.fastq.gz")
)

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtFs and filtRs goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtFs_GENE1) <- sample_names_GENE1
names(filtRs_GENE1) <- sample_names_GENE1

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

out_GENE1 <- filterAndTrim(
  trimmed_F_GENE1,
  filtFs_GENE1,
  trimmed_R_GENE1,
  filtRs_GENE1,
  truncLen = c(0, 0),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out_GENE1
head(out_GENE1)

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
exists_GENE1 <- file.exists(filtFs_GENE1) & file.exists(filtRs_GENE1)
exists_GENE1
filtFs_GENE1 <- filtFs_GENE1[exists_GENE1]
filtRs_GENE1 <- filtRs_GENE1[exists_GENE1]

# I sometimes look at the quality of the filtered reads. Sometimes my filtering
# parameters don't get rid of all the poor quality at the 3' end, but if quality
# is too poor, they won't merge, so it doesn't seem to hurt read passage.
plotQualityProfile(filtFs_GENE1[1:N], aggregate = TRUE)
plotQualityProfile(filtRs_GENE1[1:N], aggregate = TRUE)

## Estimating Error Rates and Denoising ========================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errF_GENE1 <- learnErrors(
  filtFs_GENE1,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errR_GENE1 <- learnErrors(
  filtRs_GENE1,
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
plotErrors(errF_GENE1, nominalQ = TRUE)
plotErrors(errR_GENE1, nominalQ = TRUE)

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtFs), and "err =" which is the error file from
# learnErrors (effF).
dadaFs_GENE1 <- dada(
  filtFs_GENE1,
  err = errF_GENE1,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

dadaRs_GENE1 <- dada(
  filtRs_GENE1,
  err = errR_GENE1,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
dadaFs_GENE1[[1]]
dadaRs_GENE1[[1]]

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
merged_GENE1 <- mergePairs(
  dadaFs_GENE1,
  filtFs_GENE1,
  dadaRs_GENE1,
  filtRs_GENE1,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged_GENE1[[1]])
head(merged_GENE1[[6]])

## Create and Trim Sequence-Table ==============================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab_GENE1 <- makeSequenceTable(merged_GENE1)
# This describes the dimensions of the table just made
dim(seqtab_GENE1)

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
table(nchar(getSequences(seqtab_GENE1)))

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding GENEs (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab313 <- seqtab[, nchar(colnames(seqtab)) %in% 298:322]
dim(seqtab313)
table(nchar(getSequences(seqtab313)))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab_nochim_GENE1 <- removeBimeraDenovo(
  seqtab_GENE1,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab_nochim_GENE1)

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2. I've added a column to the typical
# tutorial version of this that gives us the percentage of reads that made it
# through the process. This is a good quick way to see if something is wrong
# (i.e. only a small proportion make it through).
getN <- function(x) sum(getUniques(x))
track_GENE1 <- cbind(
  out_GENE1,
  sapply(dadaFs_GENE1, getN),
  sapply(dadaRs_GENE1, getN),
  sapply(merged_GENE1, getN),
  rowSums(seqtab_nochim_GENE1),
  100 * (rowSums(seqtab_nochim_GENE1) / out_GENE1[, 1])
)

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_GENE1) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim",
  "%kept"
)
rownames(track_GENE1) <- sample_names_GENE1

# Look at the results for the first 6 samples, or all, depending upon the number
# of samples.
head(track_GENE1)
track_GENE1
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
  seqtab_nochim_GENE1,
  file = "data/results/PROJECTNAME_GENE1_sequence-table.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)
# NOTE!!!
# The Sequence-Table in this format is very unwieldy, since each column name is
# an entire ASV. Instead, we can convert each ASV into a short "hash" using
# the md5 encryption model, creating a 32bit representative of each ASV. Each
# hash is essentially unique to the ASV it is representing. We would then
# replace the ASVs in the column headings with their representative md5 hash.
# However, having an ASV hash as a column heading requires the creation of a
# Representative Sequence list, which tells us which hash goes with which ASV.
# gives the user a representative-sequence fasta that contains the ASV, labelled
# with its specfic md5 hash.
# If you want to export a Sequence-Table with a md5 hash instead of ASV sequence
# for each ASV, skip this and go to the next section.

## Create And Use md5 Hash =====================================================
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes. Use whatever
#  table you will later use for your analyses (e.g. seqtab.nochim)
repseq_GENE1 <- getSequences(seqtab_nochim_GENE1)
# We want to look at this list, to make sure you are getting the right thing
head(repseq_GENE1)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq_md5_GENE1 <- c()
for (i in seq_along(repseq_GENE1)) {
  repseq_md5_GENE1[i] <- digest(
    repseq_GENE1[i],
    serialize = FALSE,
    algo = "md5"
  )
}
# Examine the list of feature hashes
head(repseq_md5_GENE1)


# Add md5 hash to the sequence-table from the DADA2 analysis.
seqtab_nochim_md5_GENE1 <- seqtab_nochim_GENE1
colnames(seqtab_nochim_md5_GENE1) <- repseq_md5_GENE1
View(seqtab_nochim_md5_GENE1)

# Create an md5/ASV table, with each row as an ASV and it's representative md5
# hash.
repseq_md5_asv_GENE1 <- tibble(repseq_md5_GENE1, repseq_GENE1)
# Rename column headings
colnames(repseq_md5_asv_GENE1) <- c("md5", "ASV")

head(repseq_md5_asv_GENE1)

## Export Sequence-Table with md5 Hash =========================================
# This exports a sequence-table: columns of ASV's (shown as a md5 hash instead
# of sequence), rows of samples, and values = number of reads. With this table
# you will also need a file that relates each ASV to it's representative md5 hash. We download this in the next section.

write.table(
  seqtab_nochim_md5_GENE1,
  file = "data/results/PROJECTNAME_GENE1_sequence-table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Export Representative Sequences table/fasta =================================
# Here we export our our representative sequences, either as a fasta (with the
# md5 hash as the ASV name), or as a table with ASV and md5 hash as columns.

# This exports all the ASVs in fasta format, with ASV hash as the sequence
# name. This is analogous to the representative sequence output in Qiime2.
write.fasta(
  sequences = as.list(repseq_md5_asv_GENE1$ASV),
  names = repseq_md5_asv_GENE1$md5,
  open = "w",
  as.string = FALSE,
  file.out = "data/results/PROJECTNAME_GENE1_rep-seq.fas"
)

# This exports all the ASVs and their respective md5 hashes as a two-column
# table.
write.table(
  repseq_md5_asv_GENE1,
  file = "data/results/PROJECTNAME_GENE1_representative_sequence_table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Gene2 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# trimmed_F) and the other for reverse reads (R2, called trimmed_R).
trimmed_F_GENE2 <- sort(list.files(
  path_GENE2,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))
trimmed_R_GENE2 <- sort(list.files(
  path_GENE2,
  pattern = "_R2.fastq.gz",
  full.names = TRUE
))

# Make sure you have the same number of forward reads as reverse reads
length(trimmed_F_GENE2)
length(trimmed_R_GENE2)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(trimmed_F_GENE2)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_F_GENE2_exists <- trimmed_F_GENE2[
  file.size(trimmed_F_GENE2) > 50 & file.size(trimmed_R_GENE2) > 50
]
length(trimmed_F_GENE2_exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_R_GENE2_exists <- trimmed_R_GENE2[
  file.size(trimmed_F_GENE2) > 50 & file.size(trimmed_R_GENE2) > 50
]
length(trimmed_R_GENE2_exists)
file.size(trimmed_F_GENE2_exists)

# Redefine trimmed_F and trimmed_R as only the existing read files, and check
trimmed_F_GENE2 <- trimmed_F_GENE2_exists
trimmed_R_GENE2 <- trimmed_R_GENE2_exists
length(trimmed_F_GENE2)
length(trimmed_R_GENE2)
file.size(trimmed_F_GENE2)

# Update your samples names
sample_names_GENE2 <- sapply(strsplit(basename(trimmed_F_GENE2), "_"), `[`, 1)
length(sample_names_GENE2)
head(sample_names_GENE2)

## Filter and Trim =============================================================

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets (to aggregate all samples, replace N
# with the number of samples, or with length(trimmed_F)). For example, "trimmed_F[1:2]"
# will result in two plots, one for the first sample and one for the second.
# "trimmed_F[1:17]" will result in 17 plots, one for each sample.  Using
# "aggregate = TRUE" will combine any samples called (for example, "fnFS[1:17]"
# aggregates sample 1 through 17) into a single plot. This results in the same
# the quality plots as Qiime2.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.
qualplotF_GENE2 <- plotQualityProfile(
  trimmed_F_GENE2[1:N],
  aggregate = TRUE
)
qualplotF_GENE2

# Here we modify our quality plot to better visualize where the quality cut-off
# should be. "Scale_x_continuous" is the minimum and maximum x-axis values to be
# shown.  We use "breaks=seq(a,b,c)", to indicate the first axis tick "a", last
# tick "b", and frequency of ticks "c". The example shown results in a plot that
# starts at 190 bp, ends at 220 bp, and has axis ticks every 2 bp.
qualplotF_GENE2 +
  scale_x_continuous(limits = c(250, 290), breaks = seq(250, 290, 5))

# Examine the reverse reads as you did the forward.
qualplotR_GENE2 <- plotQualityProfile(
  trimmed_R_GENE2[1:N],
  aggregate = TRUE
)
qualplotR_GENE2
qualplotR_GENE2 +
  scale_x_continuous(limits = c(250, 290), breaks = seq(250, 290, 5))


# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtFs_GENE2 <- file.path(
  path_GENE2,
  "filtered",
  paste0(sample_names_GENE2, "_F_filt.fastq.gz")
)
filtRs_GENE2 <- file.path(
  path_GENE2,
  "filtered",
  paste0(sample_names_GENE2, "_R_filt.fastq.gz")
)

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtFs and filtRs goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtFs_GENE2) <- sample_names_GENE2
names(filtRs_GENE2) <- sample_names_GENE2

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

out_GENE2 <- filterAndTrim(
  trimmed_F_GENE2,
  filtFs_GENE2,
  trimmed_R_GENE2,
  filtRs_GENE2,
  truncLen = c(0, 0),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out_GENE2
head(out_GENE2)

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
exists_GENE2 <- file.exists(filtFs_GENE2) & file.exists(filtRs_GENE2)
exists_GENE2
filtFs_GENE2 <- filtFs_GENE2[exists_GENE2]
filtRs_GENE2 <- filtRs_GENE2[exists_GENE2]

# I sometimes look at the quality of the filtered reads. Sometimes my filtering
# parameters don't get rid of all the poor quality at the 3' end, but if quality
# is too poor, they won't merge, so it doesn't seem to hurt read passage.
plotQualityProfile(filtFs_GENE2[1:N], aggregate = TRUE)
plotQualityProfile(filtRs_GENE2[1:N], aggregate = TRUE)

## Estimating Error Rates and Denoising ========================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errF_GENE2 <- learnErrors(
  filtFs_GENE2,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errR_GENE2 <- learnErrors(
  filtRs_GENE2,
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
plotErrors(errF_GENE2, nominalQ = TRUE)
plotErrors(errR_GENE2, nominalQ = TRUE)

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtFs), and "err =" which is the error file from
# learnErrors (effF).
dadaFs_GENE2 <- dada(
  filtFs_GENE2,
  err = errF_GENE2,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

dadaRs_GENE2 <- dada(
  filtRs_GENE2,
  err = errR_GENE2,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
dadaFs_GENE2[[1]]
dadaRs_GENE2[[1]]

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
merged_GENE2 <- mergePairs(
  dadaFs_GENE2,
  filtFs_GENE2,
  dadaRs_GENE2,
  filtRs_GENE2,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged_GENE2[[1]])
head(merged_GENE2[[6]])

## Create and Trim Sequence-Table ==============================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab_GENE2 <- makeSequenceTable(merged_GENE2)
# This describes the dimensions of the table just made
dim(seqtab_GENE2)

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
table(nchar(getSequences(seqtab_GENE2)))

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding genes (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab313 <- seqtab[, nchar(colnames(seqtab)) %in% 298:322]
dim(seqtab313)
table(nchar(getSequences(seqtab313)))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab_nochim_GENE2 <- removeBimeraDenovo(
  seqtab_GENE2,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab_nochim_GENE2)

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2. I've added a column to the typical
# tutorial version of this that gives us the percentage of reads that made it
# through the process. This is a good quick way to see if something is wrong
# (i.e. only a small proportion make it through).
getN <- function(x) sum(getUniques(x))
track_GENE2 <- cbind(
  out_GENE2,
  sapply(dadaFs_GENE2, getN),
  sapply(dadaRs_GENE2, getN),
  sapply(merged_GENE2, getN),
  rowSums(seqtab_nochim_GENE2),
  100 * (rowSums(seqtab_nochim_GENE2) / out_GENE2[, 1])
)

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_GENE2) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim",
  "%kept"
)
rownames(track_GENE2) <- sample_names_GENE2

# Look at the results for the first 6 samples, or all, depending upon the number
# of samples.
head(track_GENE2)
track_GENE2
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
  seqtab_nochim_GENE2,
  file = "data/results/PROJECTNAME_GENE2_sequence-table.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)
# NOTE!!!
# The Sequence-Table in this format is very unwieldy, since each column name is
# an entire ASV. Instead, we can convert each ASV into a short "hash" using
# the md5 encryption model, creating a 32bit representative of each ASV. Each
# hash is essentially unique to the ASV it is representing. We would then
# replace the ASVs in the column headings with their representative md5 hash.
# However, having an ASV hash as a column heading requires the creation of a
# Representative Sequence list, which tells us which hash goes with which ASV.
# gives the user a representative-sequence fasta that contains the ASV, labelled
# with its specfic md5 hash.
# If you want to export a Sequence-Table with a md5 hash instead of ASV sequence
# for each ASV, skip this and go to the next section.

## Create And Use md5 Hash =====================================================
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes. Use whatever
#  table you will later use for your analyses (e.g. seqtab.nochim)
repseq_GENE2 <- getSequences(seqtab_nochim_GENE2)
# We want to look at this list, to make sure you are getting the right thing
head(repseq_GENE2)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq_md5_GENE2 <- c()
for (i in seq_along(repseq_GENE2)) {
  repseq_md5_GENE2[i] <- digest(
    repseq_GENE2[i],
    serialize = FALSE,
    algo = "md5"
  )
}
# Examine the list of feature hashes
head(repseq_md5_GENE2)


# Add md5 hash to the sequence-table from the DADA2 analysis.
seqtab_nochim_md5_GENE2 <- seqtab_nochim_GENE2
colnames(seqtab_nochim_md5_GENE2) <- repseq_md5_GENE2
View(seqtab_nochim_md5_GENE2)

# Create an md5/ASV table, with each row as an ASV and it's representative md5
# hash.
repseq_md5_asv_GENE2 <- tibble(repseq_md5_GENE2, repseq_GENE2)
# Rename column headings
colnames(repseq_md5_asv_GENE2) <- c("md5", "ASV")

head(repseq_md5_asv_GENE2)

## Export Sequence-Table with md5 Hash =========================================
# This exports a sequence-table: columns of ASV's (shown as a md5 hash instead
# of sequence), rows of samples, and values = number of reads. With this table
# you will also need a file that relates each ASV to it's representative md5
# hash. We download this in the next section.

write.table(
  seqtab_nochim_md5_GENE2,
  file = "data/results/PROJECTNAME_GENE2_sequence-table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Export Representative Sequences table/fasta =================================
# Here we export our our representative sequences, either as a fasta (with the
# md5 hash as the ASV name), or as a table with ASV and md5 hash as columns.

# This exports all the ASVs in fasta format, with ASV hash as the sequence
# name. This is analogous to the representative sequence output in Qiime2.
write.fasta(
  sequences = as.list(repseq_md5_asv_GENE2$ASV),
  names = repseq_md5_asv_GENE2$md5,
  open = "w",
  as.string = FALSE,
  file.out = "data/results/PROJECTNAME_GENE2_rep-seq.fas"
)

# This exports all the ASVs and their respective md5 hashes as a two-column
# table.
write.table(
  repseq_md5_asv_GENE2,
  file = "data/results/PROJECTNAME_GENE2_representative_sequence_table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
