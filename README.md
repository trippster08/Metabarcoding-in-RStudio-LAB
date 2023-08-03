# Metabarcoding Pipeline in RStudio for LAB
1. [Computer and RStudio Preparation](#Computer-and-RStudio-Preparation) </br>
  1.1. [Install/Update Computer Programs](#Install/Update-Computer-Programs) </br>
  1.2. [Get Raw Reads](#Get-Raw-Reads) </br>
  1.3. [RStudio Preparation](#RStudio-Preparation) </br>
2. [Cutadapt](#Cutadapt) </br>
  2a. [Cutadapt-trim](#Cutadapt-trim) <br/>
  2b. [Cutadapt-trim and demultiplex genes](#Cutadapt-trim-and-demultiplex-genes) <br/>
3. [DADA2](#DADA2) </br>
  3a. [DADA2 single gene](#DADA2-multiple-genes) <br/>
  3b. [DADA2 multiple genes](#DADA2-multiple-genes) <br/>
4. [Reformat and Export Files](#Reformat-and-Export-Files) </br>
5. [Import and Combine Files](#Import-and-Combine-Files) </br>
6. [Assign Taxonomy](#Assign-Taxonomy) </br>
7. [Phyloseq](#Phyloseq) </br>

This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on your computer, not on Hydra. It is broken up into sections, each section an `.R` document that can be opened in RStudio. Once in RStudio, each command can be run using the `Run` button, or with `control + return`. The directions for each section are in that section file. You can download this entire pipeline, including the RStudio files using this link: [Metabarcoding Pipeline - RStudio Documents](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/archive/refs/heads/multi-gene.zip). I usually download a version of this pipeline for each run I analyse (in case any changes need to be made, and so the primer folder is in the correct place) and save it in the working directory of that run.

However, before running RStudio, you must make sure the necessary programs are installed, and the illumina demultiplexed sequences have been downloaded.

## 1 - Computer and RStudio Preparation
### Install/Update Computer Programs
Make sure you have both R and RStudio already installed and updated on your computer. If you have an SI computer, you can load/update both through the Smithsonian's Self Service Application.

#### Install miniconda
To install miniconda, go to `https://docs.conda.io/en/latest/miniconda.html` and download the Mac OS X 64-bit (bash installer). Open a new command window, go to Downloads/, and run the downloaded shell script to install conda.

Below is just an example, so replace `Miniconda3-latest-MacOSX-x86_64.sh` with the downloaded file name.
```
sh Miniconda3-latest-MacOSX-x86_64.sh
```

Find and enter the folder containing conda (mine is ~/miniconda3) cd ~/miniconda3

#### Install biopython
```
conda install -c conda-forge biopython
```
Add bioconda channel and other channels needed for bioconda. Run this in the order shown (this sets priority, with highest priority last).
```
conda config --add channels defaults
```
```
conda config --add channels bioconda
```
```
conda config --add channels conda-forge
```
#### Install cutadapt
We install cutadapt and create a cutadapt conda environment (called cutadaptenv) simultaneously. I typically check the anoconda webpage for cutadapt https://anaconda.org/bioconda/cutadapt and specify the version listed. The most up-to-date version is not always installed if not specified. v4.4 is an example, replace with the current version.
```
conda create -n cutadapt cutadapt=4.4
```
You may get an error telling you that cutadapt 4.4 does not exist or cannot be found. This typically happens when installing on a Mac with M1/M2 architecture. In this case, you have to use an altered version of this code.
```
CONDA_SUBDIR=osx-64 conda create -n cutadapt cutadapt=4.4
```
Installation usually only has to be done once for your computer. Periodically you may want to update these programs.

#### Update conda and cutadapt
If you have not updated either conda or cutadapt in a while, you may want to do this first. You can update conda from your home directory or anywhere. To update cutadapt, you must first activate the cutadapt enviroment. Sometimes this does not update the environment correctly. If it does not, delete the environment and reinstall you environment as you did the first time.
```
Update conda
```
```
conda update conda
```
Activate the cutadapt environment, update it, and check the updated version.
```
conda activate cutadapt
```
```
conda update cutadapt
```
```
cutadapt --version
```
If the updated version is not the latest version shown on the webpage: https://anaconda.org/bioconda/cutadapt, delete the current
environment
```
conda deactivate
```
```
conda remove -n cutadapt --all
```
Check to make sure its gone
```
conda env list
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
### Get Raw Reads
Get raw reads. If you download them directly using BaseSpace Downloader, it creates a directory that will be your working directory for this pipeline.  Place this directory where you want to run the pipeline.  If you get reads already downloaded in a folder, then either place that folder wherever you want to run this pipeline, and use it as the working directory, or make a new run-specific directory to place these sequences. Below is an example for making a run-specific directory to place your sequences. I perform all my metabarcoding analyses directly in the BaseSpace-downloaded directory within a parent directory called "/Projects_Metabarcoding".

DO NOT USE THIS COMMAND AS-IS, USE YOUR OWN PATH, INCLUDING YOUR USERNAME AND WHATEVER PROJECT NAME YOU WANT TO USE!
```
mkdir -p /Users/USERNAME/Dropbox\ \(Smithsonian\)/Projects_Metabarcoding/PROJECTNAME.
```
Download the folder containing your raw reads into this project directory using basespace.

The rest of this pipeline is run through RStudio.

### RStudio Preparation
The first thing to do after opening RStudio is create a new project, select to creat it from an "Existing Directory", and chose the directory that you will be using as your working directory. Once you have created this project, it will make this directory the current working directory, and you won't need to set your working directory later.

Next we install and load all the R libraries needed for this pipeline. We also set up our directory structure and find, load, and copy the raw Illumina read files to the directory from which they will be analyzed. 

Open RStudio, and open `1_Metabarcoding_R_Pipeline_RstudioPrep.R` in the Source Editor (typically the top left pane). You can run all commands in the source editor using the `Run` button or `control + return`.

[1.3 - Metabarcoding RStudioPrep.R](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/1_3_Metabarcoding_R_Pipeline_RStudioPrep.R)

## 2 - Cutadapt
We use Cutadapt to remove primer sequences from our raw reads. This section ends with primer-trimmed sequences. There are two versions of Cutadapt in this pipeline. The first version (2a) is for Illumina runs with only a single gene-product. Use the second (2b) if you have more than one gene product in your run. In this case, cutadapt will trim primers, but also sort reads dependent upon which gene-specific primers it removed (e.g. it will move reads from which it removed 18S primers into an 18S folder, and reads from which it removed COI primers into a COI folder).

[2a. - Cutadapt-trim](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/2%20Cutadapt_trim.R)
[2b. - Cutadapt-trim and demultiplex genes](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/2%20Cutadapt_demultiplex.R)

## 3 - DADA2
Here we use DADA2 to quality-filter and quality-trim reads, estimate error rates and denoise reads, merge paired reads, and remove chimeric sequences. This section ends with a sequence-table, which is a table containing columns of `ASV's` (Amplicon Sequence Variants), rows of `samples`, and cell values equal `# of reads`. There are two versions for this section of the pipeline too. Section 3a is for Illumina runs with a single target gene, where you used Cutadapt 2a. If you used Cutadapt 2b, and had multiple genes in your run, use DADA2 3b.

[3a - DADA2 single gene](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/3%20Metabarcoding_R_Pipeline_RStudio_Dada2.R)
[3b - DADA2 multiple genes](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/3%20Metabarcoding_R_Pipeline_RStudio_Dada2.R)

## 4 - Reformat and Export Files
Here we creat and export variants of the sequence-table created in section 3. Many of these are not necessary for your analysis, but may be useful in some cases, as described in the section descriptions. 

One variant is a Sequence-List table (a tidy table containing columns of `sample name`, `ASV`, and   `# of reads`. There is a separate row for each `sample name`/`ASV` combination.) This table is useful if you have sequences for multiple runs, because they can be directly concatenated in a text file and condensed back into a sequence-table for downstream analysis.

A second variant is a feature-table. This table contains columns of `sample names` with rows of `ASV's`, and cell values equal to `# of reads`. This is essentially a transposed sequence-table. It is also the output of the metabarcoding program Qiime2, and is included in case you have other analyses in this program that you may want to combine or compare. One aspect of the Qiime2 feature-table is that ASV's are not shown as entire sequences, but as a md5 hash (see section description for for information about md5 encryption), and this section will also add md5 hash information for each ASV. Even if a feature-table is not needed, it is often good to have md5 hash's available for downstreama analyses.

Finally, you can also create and export your data in a format we refer to as "feature-to-fasta". This creates a fasta file containing all the ASV's found in each sample. Each `ASV` will be labeled with the `sample name`, `ASV hash`, and `# of reads` for that `sample name`/`ASV` combination. This format is useful for making trees (espicially in low-diversity studies, or when sequencing single-organism samples) to look at distribution of ASV's across samples.

[4 - Format and Export Files](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/4%20Metabarcoding_R_Pipeline_RStudio_FormatandExportFiles.R)

## - 5 Import and Combine Files
Here we import and combine trimming/denoising results from multiple runs into a single project table for downstream analyses.  The specific procedure used depends upon the format of the information being imported and combined. This section is only needed for importing denoised data for analysis, or if combining data from multiple Illumina runs.

[5 - Import and Combine Files](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/5%20Metabarcoding_R_Pipeline_RStudio_ImportCombine.R)

## - 6 Assign Taxonomy
Here we use DADA2 to assign taxonomic identities to ASV's. This section requires a reference library. LAB has libraries available for both COI and 12S, but you may want to use your own. How to do so is described in the section description.

[6 - Assign Taxonomy](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/6%20Metabarcoding_R_Pipeline_RStudio_TaxAssignment.R)

## - 7 Phyloseq
Phyloseq is a R library that allows for manipulation, visualization, and analysis of metabarcoding data. This section describes how to set up and load your denoised results from DADA2 into Phyloseq, how to perform some preliminary analyses, ana how to visualize a few basic results.

[7 - Phyloseq](https://github.com/trippster08/Metabarcoding-in-RStudio-LAB/blob/main/RStudio%20Files/7%20Metabarcoding_R_Pipeline_RStudio_Phyloseq.R)
