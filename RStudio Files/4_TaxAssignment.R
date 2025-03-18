# 1June2022
# 8 ASSIGN TAXONOMY ############################################################

## Load Libraries = ============================================================
# Load all R packages you may need, if necessary

library(dada2)
library(digest)
library(rBLAST)
library(tidyverse)
library(seqinr)
library(R.utils)

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.
setwd(
  "/Users/USERNAME/Dropbox (Smithsonian)/Projects_Metabarcoding/PROJECTNAME"
)

## Get Reference Database ======================================================
# If you have your own reference database, you can either enter the path when
# we assign taxonomy below, or move it into the ref/ directory in the main
# project directory

# To use a supplied reference library, download into the ref folder. However, I
# strongly recommend using your own reference.

# Use this link for a full reference database
reference_url <- "https://www.dropbox.com/scl/fi/qagm6nh1chgidlzktaslr/MIDORI2_UNIQ_NUC_GB260_CO1_DADA2_noInsect.fasta.gz?rlkey=ibm0x4k6bspt31u9ekq7gvizd&dl=1"
# Use this link for a much reduced reference database
reference_url <- "https://www.dropbox.com/scl/fi/885yhlno5kwcz0eja5vbe/midori_COI_genus_dada2.fasta.gz?rlkey=pf80nkf9tf3kpwwno2os584o9&dl=1"
# Specify destination file and download
dest_file <- sub("\\.fasta\\.gz.*$", ".fasta.gz", basename(reference_url))
download.file(reference_url, paste0("ref/", dest_file), mode = "wb")

# Specify decompressed file and unzip
decompressed_file <- sub("\\.gz$", "", paste0("ref/", dest_file))
gunzip(paste0("ref/", dest_file), destname = decompressed_file, remove = TRUE)

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
  decompressed_file,
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
