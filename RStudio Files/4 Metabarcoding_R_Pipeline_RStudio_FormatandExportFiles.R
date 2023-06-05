# 4 FORMAT AND EXPORT DADA2 RESULTS ############################################

# This section of the gives lots of options for creating and exporting various
# ways of looking at the DADA2 results. This section is not necessary if you
# just want to use your results directly in taxonomic assignment or phylip or
# another program. Each format below has a description that mentions why you may
# want to use a particular format.

## Load Libraries = ============================================================
# Load all R packages you may need if necessary.
library(dada2)
library(digest)
library(phyloseq)
library(tidyverse)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)
library(filesstrings)

# Set your working directory to you project directory if necessary.
setwd("/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")

## Create and Export Sequence-List Table
# This creates a table containing columns of sample name, ASV, and read
# count. Each row is a separate sample. This is a good table format for
# storage of DADA2 results because it can be easily concatenated
# with other sequence-list tables in Excel or any text-editing software
# (unlike the sequence-table), yet it still contains all the information
# needed from our trimming and denoising steps. It is also a tidy table.

# You will also need to make a Sequence-List table if you want to export a
# feature-to-fasta file later. This fasta file contains every combination of ASV
# and sample, each sequence named with the sample, md5 hash (called a feature in
# Qiime2) and number of reads of that ASV in that sample.

# Convert the sequence-table from your DADA2 results into a tibble,
# converting row names to column 1, labeled "sample". A tibble is a more
# versatile data.frame, but it does not have row headings
# (among other differences, see https://tibble.tidyverse.org/). We'll need
# this to be a tibble for the next step.
seqtab.nochim.tb <- as_tibble(seqtab.nochim, rownames = "sample")

# The sequence-table has a column with sample names, and N columns of ASV's
# containing count values. We want all the count data to be in a single column,
# so we use a tidyr command called "pivot_longer" to make the table "tall",
# which means the table goes from 17x2811 to 47770x3 for example
# (47770 = 2810 x 17. 2810 instead of 2811 because the first column of the
# original table contains sample names, not counts). This makes the table tidier
# (meaning that each column is now a true variable).
seqtab.nochim.tall <- seqtab.nochim.tb %>%
  pivot_longer (
    !sample,
    names_to = "ASV",
    values_to = "count"
  )
# Look at your new table.
head(seqtab.nochim.tall)
# Look at the dimensions of this table.
dim(seqtab.nochim.tall)

# Remove rows with sequence counts = 0. This removes any samples in which a
# particular ASV was not found.
seqtab.nochim.tall.nozero <- subset(seqtab.nochim.tall, count != 0)
# Look at the dimensions of this table to compare to the previous.
dim(seqtab.nochim.tall.nozero)

# Export Sequence-List Table. This includes 3 columns: sample name, ASV
# sequence, and read count for that sample name ASV combo. This is a tidy table
# (each column a variable), and is a good way to include all the possible data
# involved.
write.table(
  seqtab.nochim.tall.nozero,
  file="data/results/PROJECTNAME_SeqList_Tall.tsv",
  quote = FALSE,
  sep="\t",
  row.names = FALSE
)

## Format Qiime2-Formated Results ==============================================
# Dada2 outputs a sequence-table, with columns of ASVs and rows of samples
# (which I call the sequence-table throughout). Qiime2 outputs a table with
# columns of samples and rows of ASVs (which they call a feature-table,
# and henceforth called such here). However, tables with ASVs are not very
# human-readable. To fix this, Qiime2 converts each ASV to a short "hash" using
# the md5 encryption model, creating a 32bit representative of each ASV. Each
# hash is essentially unique to the ASV it is representing. The Qiime2
# feature-table replaces each ASV with it's representative md5 hash. Qiime2 also
# gives the user a representative-sequence fasta that contains the ASV, labelled
# with its specfic md5 hash.

# When we transpose the sequence-table into a feature-table, we will also
# convert this matrix into a tibble (using "as_tibble").
# Using 'rownames = "ASV"' converts the columns names of the original matrix
# (which are ASV's) into a new column of the transposed tibble, called "ASV". We
# are also converting this matrix into a tibble because we need ASV as a column
# later, but this column is of a different data type than the rest of the
# columns, and matrices don't allow mixed data types. Tibbles are also easier
# to manipulate. The sequence-table is more commonly used than the feature-table
# (especially in ecology), I believe, but most downstream programs will allow
# either format, you just have to tell the program which one you use.

# Once we have the hashes, we are going to do multiple things with them. We can
# change the sequence-table so that the column headings are ASV hashes instead
# of ASV's, to make for easier viewing. We can also replace the sequences with
# feature hashes in the first column of the feature-table, to make it mimic the
# "feature-table" output of Qiime2. We can also use both the hash list and
# representative sequence list to output a fasta file of all ASV's,
# mimicking the rep_seq output of Qiime2. Next, we can combine the Sample Name,
# Feature hash, and count of each cell from the feature-table, and output a
# fasta file containing sample-specific sequences, mimicking Matt Kweskins
# program "featuretofasta.py". Finally, we can add md5 hashes as a 4th column
# in our sequence-list table.

# Create Feature-Table =========================================================
# Transpose the sequence-table, and convert the result into a tibble.
seqtab.nochim.transpose <- as_tibble(t(seqtab.nochim), rownames = "ASV")

# Check to make sure the table is transposed. The easiest way is just to look at
# the column headings and see if they are now samples (plus "ASV"), as they
# should be.
colnames(seqtab.nochim.transpose)

## Create And Use md5 Hash =====================================================
# We first create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes. Use whatever
#  table you will later use for your analyses (e.g. seqtab.nochim)
repseq <- getSequences(seqtab.nochim)
# We want to look at this list, to make sure you are getting the right thing
head(repseq)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq.md5 <- c()
for (i in seq_along(repseq)) {
  repseq.md5[i] <- digest(
    repseq[i], 
    serialize=FALSE,
    algo="md5"
  )
}
# Examine the list of feature hashes
head(repseq.md5)
tail(repseq.md5)

## Create Feature-Table With md5 Hash ==========================================
# Add the ASV hash as a column in the transposed tibble and call it feature
# (following Qiime2 convention). This table will have a column of ASV sequences,
# a column of md5 hashes that represent these ASV's, and counts of each of the
#  ASV's for each sample. %>% is how you pipe output in tidyverse. Here I am
#  piping the new tibble to a command that moves the last column (the hash's we
# just added) to the first column
seqtab.nochim.transpose.md5 <- cbind (seqtab.nochim.transpose, feature = repseq.md5) %>%
  dplyr::select(feature, everything())

# The easiest way to check to makes sure the column was added and also moved is
# to look at the column names of the new tibble
colnames(seqtab.nochim.transpose.md5)

## Export repseq.fas, Feature-Table ============================================
# Here we export the two main outputs of Qiime2, the representative sequences
# (ASV's) fasta file and the feature-table.

# This exports all the ASV's in fasta format, with Feature hash as the sequence
# name. This is analogous to the representative sequence output in Qiime2.
# Change file.out to wherever you want to save your files from this analysis.
write.fasta(
  sequences = as.list(repseq), 
  names = repseq.md5,
  open = "w",
  as.string = FALSE,
  file.out = "data/results/PROJECTNAME_rep-seq.fas"
)

# This exports a feature-table, identical to the feature-table output of Qiime2,
# with rows of feature (ASV) hashs, columns of samples, and values of number
# of reads. It is designed to ignore column 2 (ASV sequences). To export entire
# table, remove bracketed area after seqtab.nochim.transpose.md5.
write.table(
  seqtab.nochim.transpose.md5[,c(1, 3:length(seqtab.nochim.transpose.md5))],
  file="data/results/PROJECTNAME_feature-table.tsv",
  quote = FALSE,
  sep="\t",
  row.names = FALSE
)

## Create and Export Feature-to-Fasta ==========================================
# This creates a fasta file containing all the ASV's for each sample. Each ASV
#  will be labeled with the sample name, ASV hash, and number of reads of that
# ASV in that sample, almost exactly as Matt Kweskin's python script
# "featuretofasta".  The only difference is that Matt's script outputs the ASV's
# for each sample as a separate fasta file, while this outputs a single fasta
# for all the samples in the project.

# Save the ASV sequences from the sequence-list table
# (seqtab.nochim.tall.nozera) as a new list.
repseq.tall <- seqtab.nochim.tall.nozero$ASV
# Look at the top of your new list
head(repseq.tall)

# Convert the sequences into md5 hash's, as we did earlier. md5 hash's are
# consistent across jobs, meaning identical sequences from different projects or
# being converted by different programs will result in the same hash (i.e.
# hash's here will match hash's above)
repseq.tall.md5 <- c()
for (i in seq_along(repseq.tall)) {
  repseq.tall.md5[i] <- digest(
    repseq.tall[i], 
    serialize=FALSE,
    algo="md5"
  )
}
# Examine hashed sequences to make sure conversion worked
head(repseq.tall.md5)

# Attach the ASV hashes as a column (called "feature") to the tall table. The
# table should now have 4 columns, and each row of the "feature" column should
# be a md5 hash of its respective ASV.
seqtab.nochim.tall.nozero.md5 <- cbind(
  seqtab.nochim.tall.nozero,
  feature = repseq.tall.md5
  )
# Check to make sure new table contains "feature" column
colnames(seqtab.nochim.tall.nozero.md5)

# Create a new column in this table that contains "sample", "feature", and
# "count", concatenated. This is the heading for each sequence in the fasta file
# created by Matt Kweskin's script "featuretofasta.py"
seqtab.nochim.tall.nozero.md5.name <- unite(
  seqtab.nochim.tall.nozero.md5,
  sample_feature_count,
  c(sample, feature, count),
  sep = "_",
  remove = FALSE
  )
# Check to make sure the new names are correct, for at least the first 4 rows.
# And check column headings if you like.
head(seqtab.nochim.tall.nozero.md5.name)[1:4,1]
colnames(seqtab.nochim.tall.nozero.md5.name)

# Reorder columns if you want, it doesn't really matter for the next step, but
# it may matter in how you look at the table. You could even remove the "sample"
# "feature" and "count" columns since they are not needed anymore in the table
# (by using [,c(1,3)]), but it would only be for ease of visualizing the table,
# and wouldn't change the next "write.fasta" step. You can also export this
# table with reordered columns as a sequence-list table with an additional
# column of ASV hashes. Follow the directions for exporting the sequence-list
# table.
seqtab.nochim.tall.nozero.md5.name.reorder <- seqtab.nochim.tall.nozero.md5.name[ ,c(1,3,2,5,4)]
# Check reordering
colnames(seqtab.nochim.tall.nozero.md5.name.reorder)

# Create a fasta-formatted file of each row sequence (i.e. ASV), with a heading
# of "sample_feature_count".
write.fasta(
  sequences = as.list(seqtab.nochim.tall.nozero.md5.name$ASV),
  names = seqtab.nochim.tall.nozero.md5.name$sample_feature_count,
  open = "w",
  as.string = FALSE,
  file.out = "data/results/PROJECTNAME_sample-feature-count.fas"
)