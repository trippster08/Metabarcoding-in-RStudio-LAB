#15August22
# 5 IMPORT AND COMBINE SEQUENCE- AND FEATURE-TABLES ############################
# Here we import and combine trimming/denoising results from multiple runs into
  # a single project table for downstream analyses.  The specific procedure used
  # depends upon the format of the information being imported and combined.

## File Housekeeping ===========================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(phyloseq)
library(tidyverse)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)
library(filesstrings)

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd("/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME")
# Make sure all the necessary files are present.
list.files("data/results")

## Import Different File Types =================================================

### Import dada2 sequence-table ------------------------------------------------
# The sequence table (columns of ASVs and rows of sample names) is the format
# that many downstream programs use, so this does not need to be converted, and
# data imported in other formats will be converted into a sequence-table before
# merging.

# Import your sequence-table, preserving the first row as column header and the
# first column as row header. All tables are saved as tab-delimited text files,
# change "sep = '\t'" to your specific separator. We import as a matrix, as
# required by DADA2
seqtab.project.miseqrun1 <- as.matrix(
  read.delim(
    "PROJECTNAME_MISEQRUN1_sequence-table.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
)

# Look at your newly imported miseqrun1 sequence-table for project1.
View(seqtab.project.miseqrun1)

### Import and combine Qiime2 .tsv data and representative sequences -----------
# Qiime2 outputs a feature-table (columns of ASV md5 hashes) and representative
# sequence fasta files (which are ASVs with md5 hashes for names). The data from
# these will be merged into a single sequence-table.

# Import the feature-table. When we use biom convert command in Qiime2 to
# convert the .biom table into a .tsv table, it includes an extra row with
# "# Constructed from biom file". "skip = 1" skips this row when reading in this
# table. I add md5 to the name to denote that this table has uses the md5 hash
# instead of ASV.
feattab.md5.project.miseqrun2 <- read.delim(
  "PROJECTNAME_MISEQRUN2_table-dada2.tsv",
  header = TRUE,
  sep = "\t",
  skip = 1
)

# Import the representative-sequence .fasta, and convert it into a table so it
# can be merged with the feature table. enframe() converts a named vector (in
# this case, repseq.project1.miseqrun2) into a table with two columns, ASV and
# feature (the md5 hash of the ASV, converted from the vector names in
# repseq.project1.miseqrun2)
repseq.project.miseqrun2 <- getSequences(
  "PROJECTNAME_miseqrun2_rep-seqs-dada2.fasta"
) %>%
  enframe(name = "feature", value = "ASV")

# Combine the representative-sequences table with the md5 feature-table.
feattab.repseq.project1.miseqrun2 <- cbind(
  repseq.project1.miseqrun2,
  feattab.md5.project1.miseqrun2
)

# Compare the two columns of md5 hash's to make sure they agree (ensure that
# both tables have the same ASV's in the same order). Check to make sure of the
# name of the feature column from the md5 feature-table. If it is something
# other than "X.OTU.ID", change this in the command below.
feattab.repseq.project1.miseqrun2$md5.agree <- ifelse(
  feattab.repseq.project1.miseqrun2$X.OTU.ID==feattab.repseq.project1.miseqrun2$Feature, "Yes",
  ifelse(
    feattab.repseq.project1.miseqrun2$X.OTU.ID!=feattab.repseq.project1.miseqrun2$Feature, "No"
  )
)

# This counts "yes" and "no" in the column md5.agree. Hopefully there will be
# only "yes". If any "no" occurs, you will have to rearrange your tables before
# binding them together. I am not going to go into that here, because Qiime2
# should export rep-seq and feature-tables in the same order.
table(feattab.repseq.project1.miseqrun2$md5.agree)

# Prune this table to make the feature-table with ASV's instead of md5 hashes
feattab.project1.miseqrun2 <- feattab.repseq.project1.miseqrun2[
  ,
  -c(1, 3, length(feattab.repseq.project1.miseqrun2) - 0)
]

# Transpose this feature-table into a sequence table, and convert the first
# columns/rows into headings.
seqtab.project1.miseqrun2 <- transpose(
  feattab.project1.miseqrun2,
  keep.names = "Sample",
  make.names = 1)

# Convert the first column into rownames to match the format of the DADA2 output
row.names(sectab.project1.miseqrun2) <- sectab.project1.miseqrun2[, 1]

# Remove the 1st column, which is now a duplicate of the row names
seqtab.project1.miseqrun2 <- sectab.project1.miseqrun2[, -1]

# Finally, convert the data.frame into a matrix, as required by DADA2
seqtab.project1.miseqrun2 <- as.matrix(seqtab.project1.miseqrun2)

# Look at your newly miseqrun2 sequence-table for project1.
View(seqtab.project1.miseqrun2)

### Import and combine Qiime2 .biom data and representative sequences ----------

# Install and load the biomformat library
BiocManager::install("biomformat")
library(biomformat)

feattab.md5.project1.miseqrun2.biom <- read_biom("project1_miseqrun2_table-dada2.biom")
feattab.md5.project1.miseqrun2 <- as(biom_data(feattab.md5.project1.miseqrun2.biom), "matrix")

# You now have a md5 feature-table matrix. The procedure from here is similar to
# that above, when importing a md5 feature-table tsv, except here we have a
# matrix with row names and column headings, so we need to convert row names
# to a column.
# Convert row names to the first column
feattab.md5.project1.miseqrun2 <- as_tibble(
  feattab.md5.project1.miseqrun2,
  rownames = "Feature2"
)

# Import the representative-sequence .fasta, and convert it into a table so it
# can be merged with the feature table. enframe() converts a named vector (in
# this case, repseq.project1.miseqrun2) into a table with two columns, ASV and
# Feature (the md5 hash of the ASV, converted from the vector names in
# repseq.project1.miseqrun2)
repseq.project1.miseqrun2 <- getSequences(
  "project1_miseqrun2_rep-seqs-dada2.fasta") %>%
  enframe(name = "Feature", value = "ASV")

# Combine the representative-sequences table with the md5 feature-table.
feattab.repseq.project1.miseqrun2 <- cbind(
  repseq.project1.miseqrun2,
  feattab.md5.project1.miseqrun2
)

# Compare the two columns of md5 hash's to make sure they agree (ensure that
# both tables have the same ASV's in the same order). Check to make sure the
# name of the feature column from the md5 feature-table. If it is something
# other than "X.OTU.ID", change this in the command below.
feattab.repseq.project1.miseqrun2$md5.agree <- ifelse(
  feattab.repseq.project1.miseqrun2$Feature2==feattab.repseq.project1.miseqrun2$Feature, "Yes",
  ifelse(
    feattab.repseq.project1.miseqrun2$Feature2!=feattab.repseq.project1.miseqrun2$Feature, "No"
  )
)

# This counts "yes" and "no" in the column md5.agree. Hopefully there will be
# only "yes". If any "no" occurs, you will have to rearrange your tables before
# binding them together. I am not going to go into that here, because Qiime2
# should export rep-seq and feature-tables in the same order.
table(feattab.repseq.project1.miseqrun2$md5.agree)

# Prune this table to make the feature-table with ASV's instead of md5 hashes.
feattab.project1.miseqrun2 <- feattab.repseq.project1.miseqrun2[
  ,
  -c(1, 3, length(feattab.repseq.project1.miseqrun2) - 0)
]

# Transpose this feature-table into a sequence table, and convert the first
# column/row into headings.
seqtab.project1.miseqrun2 <- transpose(
  feattab.project1.miseqrun2,
  keep.names = "Sample",
  make.names = 1)

# Convert the first column into rownames to match the format of the DADA2 output
row.names(seqtab.project1.miseqrun2) <- seqtab.project1.miseqrun2[, 1]

# Remove the 1st column, which is now a duplicate of the row names
seqtab.project1.miseqrun2 <- seqtab.project1.miseqrun2[, -1]

# Finally, convert the data.frame into a matrix, as required by DADA2
seqtab.project1.miseqrun2 <- as.matrix(seqtab.project1.miseqrun2)

# Look at your newly miseqrun2 sequence-table for project1.
View(seqtab.project1.miseqrun1)


### Import and convert sequence-list tables ------------------------------------
# Often we save our data in this format (columns of sample names, ASVs, read
# counts, and md5 hashes (optional) for each ASV/sample combination). This
# contains all the information we may need for downstream analyses, it is a tidy
# table (each column a variable), and it is easier to concatenate than either
# the feature-table or sequence-table when done outside of R.

# Import a sequence-list-table.
seqlisttab.project1.miseqrun1 <- read.delim(
  "project1_miseqrun1_SeqList_Tall.tsv",
  header = TRUE,
  sep = "\t"
)

# !!!!!OPTIONAL!!!!!
# If your sequence-list-table has a column of md5 hashes, you may want to remove
# them if your other sequence-tables do not contain them.

# Remove md5 hash column. In my case this column is called "feature". Change the
# name to match your column heading
seqlisttab.project1.miseqrun1 <- subset(
  seqlisttab.project1.miseqrun1,
  select = -feature
)
# !!!!!OPTIONAL!!!!!

# Convert this tall (and tidy) table into a wide table, in form like a
# feature-table.  This uses a tidyr command, so it also makes this table into a
# tibble. We will convert to a data.frame later
feattab.project1.miseqrun1.wide <- pivot_wider(
  seqlisttab.project1.miseqrun1,
  names_from = sample,
  values_from = count,
  values_fill = 0
)

# Transpose this feature-table into a sequence table, and convert the first
# column/row into headings.
seqtab.project1.miseqrun1 <- transpose(
  feattab.project1.miseqrun1.wide,
  keep.names = "Sample",
  make.names = 1
)

# Convert the first column into rownames to match the format of the DADA2 output
row.names(seqtab.project1.miseqrun1) <- seqtab.project1.miseqrun1[, 1]

# Remove the 1st column, which is now a duplicate of the row names
seqtab.project1.miseqrun1 <- seqtab.project1.miseqrun1[, -1]

# Finally, convert the data.frame into a matrix, as required by DADA2
seqtab.project1.miseqrun1 <- as.matrix(seqtab.project1.miseqrun1)

# Look at your newly miseqrun2 sequence-table for project1.
View(seqtab.project1.miseqrun1)

## Merge Sequence-Tables =======================================================
# Here we merge the various tables we imported and converted previously. We use
# 'repeats = "sum"' to add together counts when duplicate sample names and ASV's
# are merged. If you don't want them summed, use 'repeats = "error"', and change
# the name of identical samples. While we are merging sequence-tables with ASVs
# as column headings (because these are needed for downstream analyses), this
# also works if you have replaced ASVs with md5 hashes as column headings.

seqtab.project1 <- mergeSequenceTables(
  seqtab.project1.miseqrun1,
  seqtab.project1.miseqrun2,
  repeats = "sum",
  orderBy = "abundance"
)
