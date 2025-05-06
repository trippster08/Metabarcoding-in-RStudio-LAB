# 8 ASSIGN TAXONOMY ############################################################

## Load Libraries = ============================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(rBLAST)
library(tidyverse)
library(seqinr)
library(R.utils)

## Get Reference Database ======================================================
# If you have your own reference database, you can either enter the path when
# we assign taxonomy below, or move it into the ref/ directory in the main
# project directory

# To use a supplied reference library, download into the ref folder. However, I
# strongly recommend using your own reference.

# Use this link for a full reference database based on the Midori COI database
# (www.reference-midori.info)
reference_url <- "https://www.dropbox.com/scl/fi/7ba1096klgr11n70f8zkx/MIDORI2_UNIQ_NUC_GB260_CO1_DADA2_noInsect.fasta.gz?rlkey=ehr1f9de5x1llr05u89ko07qi&dl=1"
# Use this link for a much reduced reference database. This database only
# contains a single sequence for each Genus in the original Midori pipeline. It
# will not give you great taxonomic assignment, and is really only useful for
# demonstrative purposes (to see if it works) or for a very rough assignment.
reference_url <- "https://www.dropbox.com/scl/fi/1d2ibge46ytbpqoxxry5w/midori_COI_genus_dada2.fasta.gz?rlkey=9lknfo5e03ndzfg86qrml0rg2&dl=1"
# Specify destination file and download
ref_gzip <- sub("\\.fasta\\.gz.*$", ".fasta.gz", basename(reference_url))
download.file(reference_url, paste0("ref/", ref_gzip), mode = "wb")

# Specify decompressed file and unzip
reference_fasta <- sub("\\.gz$", "", paste0("ref/", ref_gzip))
gunzip(paste0("ref/", ref_gzip), destname = reference_fasta, remove = TRUE)

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
# If you used your own database, change the path after "seqtab_nochim" to your
# reference database.

taxonomy <- assignTaxonomy(
  seqtab_nochim,
  reference_fasta,
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
# Save this object, it took a long time to get.
save(taxonomy, file = "data/working/tax_rdp.Rdata")

## Examine and Manipulate Taxonomy =============================================
# Look at the taxonomic assignments
View(taxonomy$tax)
View(taxonomy$boot)

# You can check to see all the uniqe values exist in each column
unique(taxonomy$tax[, "Phylum"])
unique(taxonomy$tax[, "Class"])
unique(taxonomy$tax[, "Order"])
unique(taxonomy$tax[, "Family"])
# You can also get see how many times each unique value exists in that column
# using table
table(taxonomy$tax[, "Phylum"])

### Combine taxonomy and bootstrap tables --------------------------------------
# Join the two tables using an inner-join with dplyr (it shouldn't matter here
# what kind of join you use since the two tables should have the exact same
# number of rows and row headings (actually, now column 1)). I amend bootstrap
# column names with "_boot" (e.g. the bootstrap column for genus would be
# "Genus_boot"). I also add the md5 hash, and rearrange the columnns
taxonomy_rdp <- inner_join(
  as_tibble(taxonomy$tax, rownames = "ASV"),
  as_tibble(taxonomy$boot, rownames = "ASV"),
  by = "ASV",
  suffix = c("", "_boot")
) %>%
  mutate(md5 = repseq_nochim_md5_asv$md5) %>%
  select(
    md5,
    ASV,
    Phylum,
    Phylum_boot,
    Class,
    Class_boot,
    Order,
    Order_boot,
    Family,
    Family_boot,
    Genus,
    Genus_boot,
    species,
    species_boot
  )
# Look at this table
View(taxonomy_rdp)

# Save these objects
save(
  taxonomy,
  taxonomy_rdp,
  file = "data/working/tax_rdp.RData"
)
# Export this table as a .tsv file. I name it with Project Name,
# the reference library used, and taxonomy.

write.table(
  taxonomy_rdp,
  file = paste0("data/results/", project_name, "taxonomy.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Assign Taxonomy With BLAST+ =================================================

# We can also assign taxonomy using BLAST. Here we will use the program rBLAST
# to identify our ASVs. rBLAST allows you to connect directly to the NCBI
# server, or use a locally saved refernce database (in BLAST format)
# One of the reasons I'm using rBLAST is that it has a command to make a
# BLAST-formatted database from a fasta file.

# Make the blast database from the supplied DADA2-formatted Midori reference
#  database, or whatever database you are using.
# The first line is the name of the reference database, change to the database
# you are using. The second line is the name you want to give your BLAST
# database. Currently, this is set up to create a directory in ref/ with the
# name you chose, then will create all the necessary database files in that
# directory, all using the chosen name also. I use the same directory name and
# file name, you do not have to.
makeblastdb(
  reference_url,
  db_name = "ref/midori_COI/midori_COI",
  dbtype = "nucl"
)

# Next we load this database into R in the correct format for rBLAST. Give
# the relative path to the database, and include the name you gave the database
# in the previous step. I.e. db should be the same here as db_name is aboove
midori_coi_db <- blast(db = "ref/midori_COI_genus/midori_COI_genus")

# We need to have our representative sequences (the sequences we are going to
# blast). We have to reformat our representative-sequence table to be a named
# vector. Here is our representative sequence table from DADA2.
View(repseq_nochim_md5_asv)

# Make a DNAStringSet object from our representative sequences
sequences_dna <- DNAStringSet(setNames(
  repseq_nochim_md5_asv$ASV,
  repseq_nochim_md5_asv$md5
))
View(sequences_dna)

# You can also get this from the fasta file we downloaded earlier.
sequences_fasta <- readDNAStringSet("data/results/PROJECTNAME_rep-seq.fas")

# They make the same thing.
head(sequences_dna)
head(sequences_fasta)

# Finally, we blast our representative sequences against the database we created
taxonomy_blast <- predict(
  midori_coi_db,
  sequences_dna,
  outfmt = "6 qseqid sseqid pident",
  BLAST_args = "-perc_identity 85 -max_target_seqs 1 -qcov_hsp_perc 80"
)
View(taxonomy_blast)

# Now lets combine the two taxonomy tables to see how the the two methods
# compare
taxonomy_rdp_blast <- left_join(
  taxonomy_rdp,
  taxonomy_blast,
  join_by(ASV == qseqid)
)


# Export this table as a .tsv file. I name it with Project Name,
# the reference library used, and taxonomy (vs. speciesID).
write.table(
  taxonomy,
  file = "data/results/PROJECTNAME_REFERENCE_taxonomy.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

save(
  midori_coi_db,
  sequences_dna,
  taxonomy_blast,
  taxonomy_rdp_blast,
  file = "tax_rdp_blast.RData"
)
