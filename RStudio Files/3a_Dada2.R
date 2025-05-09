# 4 DADA2 ######################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous
# step.
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

# Set a path to the directory with the cutadapt-trimmed reads.
path_to_trimmed <- "data/working/trimmed_sequences"

# This lists the files inside the selected folder.
list.files(path_to_trimmed)

# This creates two vectors. One contains the names for forward reads (R1, called
# trimmed_F) and the other for reverse reads (R2, called trimmed_R).
trimmed_F <- sort(list.files(
  path_to_trimmed,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))
trimmed_R <- sort(list.files(
  path_to_trimmed,
  pattern = "_R2.fastq.gz",
  full.names = TRUE
))

# Make sure you have the correct number of samples, and that they match the
# number of sample names you made in the previous section (2a or 2b: Cutadapt).
length(trimmed_F)
length(trimmed_R)

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed <- sapply(strsplit(basename(trimmed_F), "_"), `[`, 1)
head(sample_names_trimmed)
length(sample_names_trimmed)

# Now count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed <- sapply(trimmed_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})
names(sequence_counts_trimmed) <- sample_names_trimmed
# Look at this
sequence_counts_trimmed

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

## Remove empty sample files ===================================================
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_F_exists <- trimmed_F[
  file.size(trimmed_F) > 50 & file.size(trimmed_R) > 50
]
length(trimmed_F_exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_R_exists <- trimmed_R[
  file.size(trimmed_F) > 50 & file.size(trimmed_R) > 50
]
length(trimmed_R_exists)
file.size(trimmed_F_exists)

# Redefine trimmed_F and trimmed_R as only the existing read files, and check
trimmed_F <- trimmed_F_exists
trimmed_R <- trimmed_R_exists
length(trimmed_F)
length(trimmed_R)
file.size(trimmed_F)

# Update your samples names
sample_names_trimmed <- sapply(strsplit(basename(trimmed_F), "_"), `[`, 1)
length(sample_names_trimmed)
head(sample_names_trimmed)

## Filter and Trim =============================================================

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets. For example, "trimmed_F[1:2]"
# will result in two plots, one for the first sample and one for the second.
# "trimmed_F[1:17]" will result in 17 plots, one for the first 17 samples.
# Using "aggregate = TRUE" will combine any samples called.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.
quality_plot_F <- plotQualityProfile(
  trimmed_F[1:length(sample_names_trimmed)],
  aggregate = TRUE
)
# Set up background for building a better quality plot
plot_build <- ggplot_build(quality_plot_F)
x_axis_range <- plot_build$layout$panel_params[[1]]$x.range
max_x <- x_axis_range[2]

# Here we modify our quality plot to better visualize where the quality cut-off
# should be. "Scale_x_continuous" sets the minimum and maximum x-axis values to
# be shown.  We use "breaks=seq(a,b,c)", to indicate the first axis tick "a",
# last tick "b", and frequency of ticks "c". We add a vertical line at every
# tick to make it easier to see where to trim
quality_plot_F_enhanced <- quality_plot_F +
  scale_x_continuous(
    limits = c(0, max_x),
    breaks = seq(0, max_x, 10)
  ) +
  geom_vline(
    xintercept = seq(0, max_x, 10),
    color = "blue",
    linewidth = 0.25
  )
quality_plot_F_enhanced
# Export this plot as a pdf.
ggsave(
  paste0("data/results/", project_name, "_quality_plot_F.pdf"),
  plot = quality_plot_F_enhanced,
  width = 9,
  height = 9
)

# Examine the reverse reads as you did the forward.
quality_plot_R <- plotQualityProfile(
  trimmed_R[1:length(sample_names_trimmed)],
  aggregate = TRUE
)


# Save all the objects created to this point in this section
save(
  path_to_trimmed,
  trimmed_F,
  trimmed_R,
  sequence_counts_trimmed,
  sample_names_trimmed,
  quality_plot_F,
  quality_plot_R,
  file = "data/working/2_qual.Rdata"
)


# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtered_F <- file.path(
  path_to_trimmed,
  "filtered",
  paste0(
    sample_names_trimmed,
    "_F_filt.fastq.gz"
  )
)
filtered_R <- file.path(
  path_to_trimmed,
  "filtered",
  paste0(
    sample_names_trimmed,
    "_R_filt.fastq.gz"
  )
)

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtered_F and filtered_R goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtered_F) <- sample_names_trimmed
names(filtered_R) <- sample_names_trimmed

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(X,Y)" is how you tell Dada2 where to truncate all forward (X) and
# reverse (Y) reads. Using "0" means reads will not be truncated.
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

filtered_summary <- filterAndTrim(
  trimmed_F,
  filtered_F,
  trimmed_R,
  filtered_R,
  truncLen = c(X, Y),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
filtered_summary
head(filtered_summary)

# After filtering, if you have any samples with no reads, you much remove them.
# This step changes filtered_F and filtered_R to only contain the names of
# samples with reads. Do this only if there are samples in
# "filtered_summary" with zero reads.
exists <- file.exists(filtered_F) & file.exists(filtered_R)
exists
filtered_F <- filtered_F[exists]
filtered_R <- filtered_R[exists]

# Set a path to the directory with the dada2-filtered reads.
path_to_filtered <- "data/working/trimmed_sequences/filtered"

# Get sample names for filtered reads
sample_names_filtered <- sapply(strsplit(basename(filtered_F), "_"), `[`, 1)
sample_names_filtered

# Count how many reads remain in each sample after filtering
sequence_counts_filtered <- sapply(filtered_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})

# Name the counts with sample names
names(sequence_counts_filtered) <- sample_names_filtered
sequence_counts_filtered

# Save all the objects created to this point
save(
  path_to_filtered,
  filtered_F,
  filtered_R,
  filtered_summary,
  sample_names_filtered,
  sequence_counts_filtered,
  file = "data/working/3_filtered_summary.Rdata"
)

# Export out as a tsv
write.table(
  filtered_summary,
  file = paste0("data/results/", project_name, "filtered_read_count.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Estimating Error Rates and Denoising ========================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errors_F <- learnErrors(
  filtered_F,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errors_R <- learnErrors(
  filtered_R,
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
error_plots_F <- plotErrors(errors_F, nominalQ = TRUE)
error_plots_F
error_plots_R <- plotErrors(errors_R, nominalQ = TRUE)
error_plots_R
# Export both plots as pdfs
ggsave(
  paste0("data/results/", project_name, "errorplots_F.pdf"),
  plot = error_plots_F,
  width = 9,
  height = 9
)
ggsave(
  paste0("data/results/", project_name, "errorplots_R.pdf"),
  plot = error_plots_R,
  width = 9,
  height = 9
)

# Save the objects created since filtered_summary
save(
  errors_F,
  errors_R,
  error_plots_F,
  error_plots_R,
  file = "data/working/4_errors.Rdata"
)

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).
denoised_F <- dada(
  filtered_F,
  err = errors_F,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

denoised_R <- dada(
  filtered_R,
  err = errors_R,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
denoised_F[[1]]
denoised_R[[1]]

# Save both the objects created in the denoising step
save(
  denoised_F,
  denoised_R,
  file = "data/working/5_denoise.RData"
)

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (denoised_F), then the forward filtered and truncated reads (filtered_F),
# then the same for the reverse reads (denoised_R and filtered_R).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged_reads <- mergePairs(
  denoised_F,
  filtered_F,
  denoised_R,
  filtered_R,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged_reads[[1]])
head(merged_reads[[6]])

## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab <- makeSequenceTable(merged_reads)
# This describes the dimensions of the table just made
dim(seqtab)

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab_nochim)

# Make a list of the ASVs that are considered chimeras, in case you want to look
# at them later
chimeras_list <- isBimeraDenovoTable(
  seqtab,
  multithread = TRUE,
  verbose = TRUE
)
# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2.
repseq_all <- getSequences(seqtab)
# Get a list of just chimera ASVs
repseq_chimera <- repseq_all[chimeras_list]

# Export this as a fasta
write.fasta(
  sequences = as.list(repseq_chimera),
  names = repseq_chimera,
  open = "w",
  as.string = FALSE,
  file.out = paste0("data/results/", project_name, "_rep-seq_chimeras.fas")
)

## Examine Sequence Lengths and Trim ===========================================

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
seq_length_table <- table(nchar(getSequences(seqtab_nochim)))
# Look at the table
seq_length_table
# Export this table as a .tsv
write.table(
  seq_length_table,
  file = paste0("data/results/", project_name, "_ASV_lengths_table.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding genes (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab_nochim_313 <- seqtab_nochim[,
  nchar(colnames(seqtab_nochim)) %in% 298:322
]
dim(seqtab_nochim_313)
table(nchar(getSequences(seqtab_nochim_313)))

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2.This is a good quick way to see if
# something is wrong (i.e. only a small proportion make it through).

# First let remind ourselves what the other sequence count objects look like
sequence_counts_raw
sequence_counts_trimmed
sequence_counts_filtered

# Now lets make a table for the post-filtered samples, including denoised,
# merged, and non-chimeric read counts
getN <- function(x) sum(getUniques(x))
sequence_counts_postfiltered <- as_tibble(
  cbind(
    sapply(denoised_F, getN),
    sapply(denoised_R, getN),
    sapply(merged_reads, getN),
    rowSums(seqtab_nochim)
  ),
  .name_repair = "unique"
) %>%
  mutate(Sample_ID = sample_names_filtered) %>%
  select(
    Sample_ID,
    Denoised_Reads_F = ...1,
    Denoised_Reads_R = ...2,
    Merged_Reads = ...3,
    Non_Chimeras = ...4
  )

# Lets look at these to make sure it makes sense
head(sequence_counts_postfiltered)

# Finally, we are going to put all this read count information together in one
# large table
track_reads <- tibble(
  sequence_counts_raw,
  Sample_ID = names(sequence_counts_raw)
) %>%
  left_join(
    tibble(sequence_counts_trimmed, Sample_ID = names(sequence_counts_trimmed)),
    join_by(Sample_ID)
  ) %>%
  left_join(
    tibble(
      sequence_counts_filtered,
      Sample_ID = names(sequence_counts_filtered)
    ),
    join_by(Sample_ID)
  ) %>%
  left_join(
    sequence_counts_postfiltered,
    join_by(Sample_ID)
  ) %>%
  mutate(Proportion_Kept = Non_Chimeras / sequence_counts_raw) %>%
  select(Sample_ID, everything()) %>%
  select(
    Sample_ID,
    Raw_Reads = sequence_counts_raw,
    Trimmed_Reads = sequence_counts_trimmed,
    Filtered_Reads = sequence_counts_filtered,
    everything()
  )


# Look at the results for the first 6 samples, or all, depending upon the number
# of samples.
head(track_reads)
track_reads

# Export this table as a .tsv
write.table(
  track_reads,
  file = paste0("data/results/", project_name, "_track_reads_table.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

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
  seqtab_nochim,
  file = paste0("data/results/", project_name, "_sequence-table.tsv"),
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

# This makes a new vector containing all the non-chimeric ASV's
# (unique sequences) returned by dada2. We are going to use this list to create
# md5 hashes. Use whatever table you will later use for your analyses
# (e.g. seqtab.nochim)
repseq_nochim <- getSequences(seqtab_nochim)
# We want to look at this list, to make sure you are getting the right thing
head(repseq_nochim)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq_nochim_md5 <- c()
for (i in seq_along(repseq_nochim)) {
  repseq_nochim_md5[i] <- digest(
    repseq_nochim[i],
    serialize = FALSE,
    algo = "md5"
  )
}
# Examine the list of feature hashes
head(repseq_nochim_md5)

# Add md5 hash to the sequence-table from the DADA2 analysis.
seqtab_nochim_md5 <- seqtab_nochim
colnames(seqtab_nochim_md5) <- repseq_nochim_md5
View(seqtab_nochim_md5)

# Export this sequence table with column headings as md5 hashs instead of ASV
# sequences
write.table(
  seqtab_nochim_md5,
  file = paste0("data/results/", project_name, "_sequence-table-md5.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# Create an md5/ASV table, with each row as an ASV and it's representative md5
# hash.
repseq_nochim_md5_asv <- tibble(repseq_nochim_md5, repseq_nochim)
# Rename column headings
colnames(repseq_nochim_md5_asv) <- c("md5", "ASV")

head(repseq_nochim_md5_asv)

## Create Feature-Table from Sequence-Table ====================================
# seqtab.nochim.md5 is a sequence table, with ASV columns and sample rows. Most
# people seem to use Feature-Tables, with ASV rows and sample columns, so here
# we create a Feature-Table

# Transpose the sequence-table, and convert the result into a tibble.
seqtab_nochim_transpose_md5 <- as_tibble(
  t(seqtab_nochim_md5),
  rownames = "ASV"
)

# Check to make sure the table is transposed. The easiest way is just to look at
# the column headings and see if they are now samples (plus "ASV"), as they
# should be.
colnames(seqtab_nochim_transpose_md5)

# Save all the objects created between denoise and here
save(
  merged_reads,
  seqtab,
  seqtab_nochim,
  chimeras_list,
  repseq_all,
  repseq_chimera,
  getN,
  track_reads,
  seq_length_table,
  repseq_nochim,
  repseq_nochim_md5,
  seqtab_nochim_md5,
  repseq_nochim_md5_asv,
  seqtab_nochim_transpose_md5,
  file = "data/working/6_feattab.RData"
)

## Export Feature-Table with md5 Hash =========================================
# This exports a feature-table: row of ASV's (shown as a md5 hash instead
# of sequence), columns of samples, and values = number of reads. With this
# table you will also need a file that relates each ASV to it's representative
# md5 hash. We download this in the next section.

write.table(
  seqtab_nochim_transpose_md5,
  file = paste0("data/results/", project_name, "_feature-table_md5.tsv"),
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
  sequences = as.list(repseq_nochim_md5_asv$ASV),
  names = repseq_nochim_md5_asv$md5,
  open = "w",
  as.string = FALSE,
  file.out = paste0("data/results/", project_name, "_rep-seq.fas")
)

# This exports all the ASVs and their respective md5 hashes as a two-column
# table.
write.table(
  repseq_nochim_md5_asv,
  file = paste0(
    "data/results/",
    project_name,
    "_representative_sequence_md5_table.tsv"
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
