# Metabarcoding Pipeline in RStudio for LAB
1. [Computer and RStudio Preparation](#1---computer-and-rstudio-preparation) </br>
  1.1. [Install and Update Computer Programs](#install-and-update-computer-programs) </br>
  1.2. [Create a New Project](#create-a-new-project) </br>
  1.3. [RStudio Preparation](#rstudio-preparation) </br>
  1.4. [Get Raw Reads](#get-raw-reads) </br>
2. [Cutadapt](#cutadapt) </br>
3. [DADA2](#dada2) </br>
4. [Assign Taxonomy](#assign-taxonomy) </br>
5. [Visualize Results](#assign-taxonomy) </br>
6. [Phyloseq](#6---phyloseq) </br>
4. [Reformat and Export Files](#4---reformat-and-export-files) </br>
5. [Import and Combine Files](#5---import-and-combine-files) </br>



This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on your computer, not on Hydra. It is broken up into sections, each section an `.R` document that can be opened in RStudio. Once in RStudio, each command can be run using the `Run` button, or with `control + return`. The directions for each section are in that section file. You can download this entire pipeline, including the RStudio files using this link: [Metabarcoding Pipeline - RStudio Documents](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/archive/refs/heads/main.zip). I usually download a version of this pipeline for each run I analyse (in case any changes need to be made, and so the primer folder is in the correct place) and save it in the working directory of that run.

However, before running RStudio, you must make sure the necessary programs are installed, and the illumina demultiplexed sequences have been downloaded.

## 1 - Computer and RStudio Preparation
### Install and Update Computer Programs
Make sure you have both R and RStudio already installed and updated on your computer. If you have an SI computer, you can load/update both through the Smithsonian's Self Service Application.

#### Install miniforge
We are going to use the open-source version of miniconda, called miniforge, to install other necessary programs. To install miniforge you need to firts download the correct installer file for your operating system. The following code will do so.
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Next we run this script with:
```
bash Miniforge3-$(uname)-$(uname -m).sh
```
You will get several prompts, including to accept the terms, placement of the program, and whether you want to initialize conda when the terminal is opened. We usually say "yes" to all these. You will need to close and reopen the terminal for conda to be ready to install other programs. 

We installing programs, we will use the command "mamba" instead of the traditional "conda" because it is usually faster at resolving dependencies. You can still use "conda" if you like.

#### Set Channel Preferences
Add bioconda channel and other channels needed for downloading and installing programs. Run this in the order shown, or all at once (this sets priority, with highest priority last).
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
#### Install cutadapt
We install cutadapt and create a cutadapt conda environment (called cutadaptenv) simultaneously. I typically check the anoconda webpage for cutadapt https://anaconda.org/bioconda/cutadapt and specify the version listed. The most up-to-date version is not always installed if not specified. v5.0 is an example, replace with the current version.

First, we need to initialize mamba for our shell.
```
mamba init
```
You will need to close and reopen the terminal.
Next, install cutadapt as an conda environment.
```
mamba create -n cutadapt cutadapt=5.0
```
You may get an error telling you that cutadapt 5.0 does not exist or cannot be found. This typically happens when installing on a Mac with M1/M2 architecture. In this case, you have to use an altered version of this code.
```
CONDA_SUBDIR=osx-64 mamba create -n cutadapt cutadapt
```
Installation usually only has to be done once for your computer. Periodically you may want to update these programs.

#### Update conda and cutadapt
If you have not updated either conda or cutadapt in a while, you may want to do this first. You can update conda from your home directory or anywhere. To update cutadapt, you must first activate the cutadapt enviroment. Sometimes this does not update the environment correctly. If it does not, delete the environment and reinstall you environment as you did the first time.
```
Update conda
```
```
mamba updata --all
```
Activate the cutadapt environment, update it, and check the updated version.
```
mamba activate cutadapt
```
```
mamba update --all
```
```
cutadapt --version
```
If the updated version is not the latest version shown on the webpage: https://anaconda.org/bioconda/cutadapt, delete the current
environment
```
mamba deactivate
```
```
mamba remove -n cutadapt --all
```
Check to make sure its gone
```
mamba env list
```
Reinstall cutadapt with the current version
```
conda create -n cutadapt cutadapt=4.4
```
```
conda activate cutadapt
```
```
cutadapt --version
```
### Create a New Project
The first thing to do is open RStudio and create a new project. If you already havea project directory for this project, select to create it from an "Existing Directory", and chose the directory that you will be using (this folder name will be the name of the project). Or select "New Directory" and name your new directory/project. RStudio will create this directory. Once you have created this project, RStudio will make this directory the current working directory, and you won't need to set your working directory later.

### RStudio Preparation
First, we are goiong to download the entire pipeline into our project directory using the script shown below. Copy the code block below into the Console panel (usually the entire left panel, or the bottom left panel if the Source Editor is open on the top left) of RStudio and run it. This will download the pipeline unzip it, and remove the zipped file.

This is probably the only R code we will be running from the Console. We typically run all the scripts by opening each file in the Source Editor and running from there so we have a record of your analyses, including any changes made and any comments that may be needed along the way.

```{R}
pipeline <- "https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/archive/refs/heads/main.zip"
download.file(pipeline, basename(pipeline))
unzip(basename(pipeline))
file.remove(basename(pipeline))
```

Next we install all the R libraries needed for this pipeline. We also set up our directory structure and find, load, and copy the raw Illumina read files to the directory from which they will be analyzed. 

Open RStudio, and open `1_RstudioPrep.R` in the Source Editor (typically the top left pane). You can run all commands in the source editor using the `Run` button or `control + return`.

[1.3 - Metabarcoding RStudioPrep.R](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/1_RStudioPrep.R)


### Get Raw Reads
Get raw reads. Whether you get them from Illumina basespace downloader or another way, all reads should already be in a directory. Place that direcectory in the main project directory. *NOTE: Remove any "undetermined" read files from the folder containing your raw reads. You do not want to include these reads in your analyses.*

## Cutadapt
We use Cutadapt to remove primer sequences from our raw reads. This section ends with primer-trimmed sequences. There are two versions of Cutadapt in this pipeline. The first version (2a) is for Illumina runs with only a single gene-product. Use the second (2b) if you have more than one gene product in your run. In this case, cutadapt will trim primers, but also sort reads depending upon which gene-specific primers it removed (e.g. it will move reads from which it removed 18S primers into an 18S folder, and reads from which it removed COI primers into a COI folder).

[2a. - Cutadapt-trim](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/2a%20Metabarcoding_Cutadapt_trim.R) </br>
[2b. - Cutadapt-trim and demultiplex genes](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/2b%20Cutadapt_trim_and_demultiplex.R)

## DADA2
Here we use DADA2 to quality-filter and quality-trim reads, estimate error rates and denoise reads, merge paired reads, and remove chimeric sequences. This section ends with a sequence-table, which is a table containing columns of `ASV's` (Amplicon Sequence Variants), rows of `samples`, and cell values equal `# of reads`. There are two versions for this section of the pipeline too. Section 3a is for Illumina runs with a single target gene, where you used Cutadapt 2a. If you used Cutadapt 2b, and had multiple genes in your run, use DADA2 3b.

[3a - DADA2 single gene](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/3a%20Metabarcoding_Dada2.R) </br>
[3b - DADA2 multiple genes](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/3b%20Metabarcoding_Dada2_multiple_genes.R)

## Assign Taxonomy
Here we use an [RDP identifier](https://benjjneb.github.io/dada2/assign.html) through DADA2 and BLAST+ through [rBLAST](https://github.com/mhahsler/rBLAST) to assign taxonomic identities to ASV's. This section requires a reference library.  We will supply you with a reference library based on the [Midori](www.reference-midori.info) reference database, or you can supply your own. Open [4_TaxAssignment.R](RStudio_Files/4_TaxAssignment.R) and follow the directions.

## Visualize Results

Here we primarily use the program [vegan](https://github.com/vegandevs/vegan) to visualize your results. We will explore our results multiple ways. Open [5_VisualizeResults.R](RStudio_Files/5_VisualizeResults.R) and follow the directions. vegan is a very expansive diversity tool and what we do here is only a fraction of it's capabilities. [vegan vignetes](https://vegandevs.r-universe.dev/vegan) is one place to find lots of links to other aspects of the program, although it gets a little into the weeds. Most of the visualization for this pipeline is from an unaffiliated website found here: [Vegan tutorial](https://peat-clark.github.io/BIO381/veganTutorial.html).

## phyloseq
[phyloseq](https://github.com/joey711/phyloseq) is a R library that allows for manipulation, visualization, and analysis of metabarcoding data. This section describes how to set up and load your denoised results from DADA2 into phyloseq, how to perform some preliminary analyses, ana how to visualize a few basic results. Open [6_phyloseq.R](RStudio_Files/6_phyloseq.R) and follow the directions.

## Reformat and Export Files
Here we create and export variants of the sequence-table created in section 3. Many of these are not necessary for your analysis, but may be useful in some cases, as described in the section descriptions. 

One variant is a Sequence-List table (a tidy table containing columns of `sample name`, `ASV`, and   `# of reads`. There is a separate row for each `sample name`/`ASV` combination.) This table is useful if you have sequences for multiple runs, because they can be directly concatenated in a text file and condensed back into a sequence-table for downstream analysis.

Yyou can also create and export your data in a format we refer to as "feature-to-fasta". This creates a fasta file containing all the ASV's/sample combinations found (i.e. each well of the sequence-table or feature-table will have a sequence in the fasta). Each sequence will be labeled with the `sample name`, `ASV hash`, and `# of reads` for that `sample name`/`ASV` combination. This format is useful for making trees (espicially in low-diversity studies, or when sequencing single-organism samples) to look at distribution of ASV's across samples, and to visualize possible "pseudogenes".

[4 - Format and Export Files](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/4%20Metabarcoding_R_Pipeline_RStudio_FormatandExportFiles.R)

## Import and Combine Files
Here we import and combine trimming/denoising results from multiple runs into a single project table for downstream analyses.  The specific procedure used depends upon the format of the information being imported and combined. This section is only needed for importing denoised data for analysis, or if combining data from multiple Illumina runs.

[5 - Import and Combine Files](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/5%20Metabarcoding_R_Pipeline_RStudio_ImportCombine.R)




