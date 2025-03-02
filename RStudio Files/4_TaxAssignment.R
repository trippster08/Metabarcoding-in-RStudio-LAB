# 1June2022
# 8 ASSIGN TAXONOMY ############################################################

## Load Libraries = ============================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(phyloseq)
library(tidyverse)
library(seqinr)
library(ape)
library(DECIPHER)
library(ade4)

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd(
  "/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME"
)

## Assign Taxonomy With DADA2 ==================================================
# Assign taxonomy. tryRC determines whether to also include the reverse
# complement of each sequence. outputBootstraps results in a second table with
# bootstrap support values for each taxonomic assignment. minBoot gives the
# minimum bootstrap required to assign taxonomy. Use whatever table you obtained
# for your earlier analyses (e.g. seqtab.nochim). This command can take a long
# time for large reference databases and large numbers of ASVs. Your reference
# file will have a different path and name then here, please make sure to use
# the correct path and name.
# From the dada2 manual: "assignTaxonomy(...) expects a training fasta file (or
# compressed fasta file) in which the taxonomy corresponding to each sequence is
# encoded in the id line in the following fashion (the second sequence is not
# assigned down to level 6)." Note that all levels above the lowest level must
# have a space, even if there is no name for that level (I've added a third
# example to demonstrate).
#>Level1;Level2;Level3;Level4;Level5;Level6;
#ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC
#>Level1;Level2;Level3;Level4;Level5;
#CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"
#>Level1;Level2;;Level4;Level5;Level6
#CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"

# taxLevels defines what taxonomic rank each of the levels shown in the above
# example represents.

taxonomy <- assignTaxonomy(
  seqtab_nochim,
  "/Users/USERNAME/Dropbox (Smithsonian)/Metabarcoding/Reference_Libraries/REFERENCE.fasta",
  taxLevels = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "species"
  ),
  tryRC = FALSE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

## Add Species With DADA2 ======================================================
# This adds a column to your taxonomy table that contains the species name
# (Genus species) to any ASVs that are exact matches to one of the reference
# sequences (i.e. the ASV is positively identified as a known genetic species).
# tryRC determines whether to also include the reverse complement of each
# sequence. allowMultiple determines whether your sequence can be matched
# to multiple species (this should only occur if your reference database has
# identical sequences).  See dada2 manual for more on allowMultiple. The
# reference database should be in the following format (a unique ID for each
# reference sequence is not necessary):
# >ID <Genus> <species>
# ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCAT
# >ID <Genus> <species>
# CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCAT
# This is an optional step. You can also assign species without assigning
# taxonomy first (just giving an output table containing one column of species
# assignments or "NA") using the command "assignSpecies" (instead of
# "addSpecies"), and keeping everything else identical.
taxonomy_species <- addSpecies(
  taxonomy,
  "/Users/macdonaldk/Dropbox (Smithsonian)/Metabarcoding/Reference_Libraries/REFERENCE_species.fasta",
  tryRC = false
)

## Examine and Manipulate Taxonomy =============================================
# Look at the taxonomic assignments
View(taxonomy$tax)
View(taxonomy$boot)

# You can check to see all the uniqe values exist in each column
unique(taxonomy$tax[, "Phylum"])
unique(taxonomy$tax[, "Class"])
unique(taxonomy$tax[, "Order"])
unique(taxonomy$tax[, "Family"])

### Combine taxonomy and bootstrap tables --------------------------------------
# You can combine the $tax and $boot table, to see simultaneously the taxonomic
# assignment and the bootstrap support for that assignment.

# Convert taxonomy and bootstrap tables into tibbles (with "ASV" as column 1)
taxonomy_tax_tb <- as_tibble(
  taxonomy$tax,
  rownames = "ASV"
)
dim(taxonomy_tax_tb)

taxonomy_boot_tb <- as_tibble(
  taxonomy$boot,
  rownames = "ASV"
)
dim(taxonomy_boot_tb)

# Join the two tables using an inner-join with dbplyr (it shouldn't matter here
# what kind of join you use since the two tables should have the exact same
# number of rows and row headings (actually, now column 1)). I amend bootstrap
# column names with "_boot" (e.g. the bootstrap column for genus would be
# "Genus_boot")
taxonomy_tb <- inner_join(
  taxonomy_tax_tb,
  taxonomy_boot_tb,
  by = "ASV",
  suffix = c("", "_boot")
)
dim(taxonomy_tb)
View(taxonomy_tb)

# Add md5 hash from earlier. The order of ASV's is the same as the sequence-
# table, so there shouldn't be any problem, but you can always redo the md5
# hash conversion here.
taxonomy_tb_md5 <- cbind(
  taxonomy_tb,
  feature = repseq_md5
)
View(taxonomy_tb_md5)

# Rearrange columns so that the md5 hash comes first, then the ASV, then each
# classfication level followed by it's respective bootstrap column.
# fmt: skip
taxonomy_tb_md5 <- taxonomy_tb_md5[ , c(16,1,2,9,3,10,4,11,5,12,6,13,7,14,8,15)]
View(taxonomy_tb_md5)

# Export this table as a .tsv file. I name it with Project Name,
# the reference library used, and taxonomy (vs. speciesID).
write.table(
  taxonomy_tb_md5,
  file = "data/results/PROJECTNAME_REFERENCE_taxonomy.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
