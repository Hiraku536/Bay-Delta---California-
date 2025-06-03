```{r}
library(ShortRead)
library(microbiome)
library(ggplot2)
library(dada2)
library(phyloseq)
library(microViz)
library(dplyr)
library(patchwork) # for arranging groups of plots
library(Biostrings)
library(ggpubr)
library(tidyverse)
```
## Route of the directory
```{r}
pathNUBELSFBioassay <-"/Users/hhuang22/Desktop/SFBioassay"
```
## Checking the files
```{r}
list.files(pathNUBELSFBioassay)
```
## Checking the files
```{r}
fnFs <- sort(list.files(pathNUBELSFBioassay, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathNUBELSFBioassay, pattern = "_R2_001.fastq.gz", full.names = TRUE))

```
## Checking number of FWD files
```{r}
length(fnFs)
```
## Checking number of REV files
```{r}
length(fnRs)
```

### IDENTIFY LINKERS ###

# CYA359F is attached to a CS1 linker. CYA781a,bR are attached to CS2 linkers.

FWD <- "ACACTGACGACATGGTTCTACA"  
REV <- "TACGGTAGCAGAGACTTGGTCT" 
FWD
REV 
### CONFIRM LINKERS & ORIENTATION ###

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Do the sequences look as expected?

### PRE-FILTER SEQUENCES TO REMOVE THOSE WITH Ns ###

# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping
# of short primer sequences difficult. 

fnFs.filtN <- file.path(pathNUBELSFBioassay, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(pathNUBELSFBioassay, "filtN", basename(fnRs))
exists <- file.exists(fnRs.filtN) & file.exists(fnFs.filtN)
filtFs <- fnFs.filtN[exists]
filtRs <- fnRs.filtN[exists]
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = TRUE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

```
## Remove linkers

### LOAD IN CUTADAPT ###

# Cutadapt is normally run in the shell, so this code allows one to run shell commands from RStudio.

cutadapt <- "/Users/hhuang22/miniforge3/envs/QC/bin/cutadapt"
system2(cutadapt, args = "--version")

### RUN CUTADAPT ###

#We now create output filenames for the Cutadapt-ed files, and define the parameters
#we are going to give the Cutadapt command. The critical parameters are the primers,
#and they need to be in the right orientation, i.e. the FWD primer should have been
#matching the forward-reads in its forward orientation, and the REV primer should 
#have been matching the reverse-reads in its forward orientation. 
#Warning: A lot of output will be written to the screen by Cutadapt!

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n 2", # -n 2 required to remove BOTH FWD and REV from reads if both are present in one read.
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Using "-n 2" parameter here because I am assuming that some reads could contain both linkers. This algorithm should also trim off parts of the linkers if only part of the sequences are present.
# Per default, untrimmed reads will pass through into the output file.
# Joel seems to have modified the Cutadapt parameters for his Nubel analysis.

### CONFIRM SUCCESSFUL REMOVAL OF ADAPTERS ###

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                        primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                   fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
## Remove primers

setwd("/Users/hhuang22/Desktop/SFBioassay/cutadapt")
rm(list=ls())

path <- "/Users/hhuang22/Desktop/SFBioassay/cutadapt"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

### IDENTIFY PRIMERS ###

# I am using the CYA359F and CYA781a,bR primers.

FWD <- "GGGGAATYTTCCGCAATGGG"  
REV <- "GACTACWGGGGTATCTAATCCCWTT" 

# Note that there are two reverse primers. They are identical except for A/T substitutions
# at two locations. Thus, I expressed the sequence of the reverse primer as "W" at these sites per IUPAC convention. Cutadapt supports these "wildcard" base calls.

### CONFIRM PRIMERS & ORIENTATION ###

allOrients <- function(primer) {
  
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

### CHECKING FOR PRIMERS ###

# Since reads are already N-filtered, I need to change the code to look at fnFs and fnRs,
# NOT fnFs.filt and fnRs.filt.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                    primerHits, fn = fnRs[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                           fn = fnFs[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

### LOAD IN CUTADAPT ###

cutadapt <- "/Users/hhuang22/miniforge3/envs/QC/bin/cutadapt"
system2(cutadapt, args = "--version") 

### RUN CUTADAPT ###

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n 2", 
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

### CONFIRM SUCCESSFUL REMOVAL OF PRIMERS ###

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                        primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                   fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

setwd("/Users/hhuang22/Desktop/SFBioassay/cutadapt/cutadapt/Working_Reads")


### CLEAR WORKSPACE ###

rm(list=ls())

### ASSIGN PATH ###

path <- "/Users/hhuang22/Desktop/SFBioassay/cutadapt/cutadapt/Working_Reads"

list.files(path)

```{r}
### READ IN FILES ###

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtF50s <- file.path(path, paste0(sample.names, "_F_filt50.fastq.gz"))
filtR50s <- file.path(path, paste0(sample.names, "_R_filt50.fastq.gz"))
names(filtF50s) <- sample.names
names(filtR50s) <- sample.names

### PERFORM FILTERING ###

out <- filterAndTrim(fnFs, filtF50s, fnRs, filtR50s, minLen = 50, rm.phix = TRUE, compress = FALSE, verbose = TRUE)

# Note that the arguments above also remove any reads consistent with the Phi X bacteriophage genome.

### SANITY CHECK ###

# Check at least one FASTQ file that has had the linkers, primers, and short/empty reads removed. Are you detecting any unusual residuals? #

# As needed, rearrange the working directory to prepare for DADA2.


#### END OF WORKFLOW ####

#### Initial Setup ####


### Set the working directory ###

setwd("/Users/hhuang22/Desktop/SFBioassay/cutadapt/cutadapt/Working_Reads")

### Clear workspace ###

rm(list=ls())

### Load required libraries ###

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

### ASSIGN PATH ###

path <- "/Users/hhuang22/Desktop/SFBioassay/cutadapt/cutadapt/Working_Reads"

list.files(path)

### READ IN FILES ###

# First we read in the names of the fastq files, and perform some string manipulation to get lists 
# of the forward and reverse fastq files in matched order

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#### VISUALIZING READ QUALITY PROFILES ####


### VISUALIZING FORWARD READS ###

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnFs[5:8])

#The matrix indicates with which samples the plots will be made (ex. samples 1-4).

#In gray-scale is a heat map of the frequency of each quality score at each base
#position. The mean quality score at each position is shown by the green line, 
#and the quartiles of the quality score distribution by the orange lines. 
#The red line shows the scaled proportion of reads that extend to at least that 
#position (this is more useful for other sequencing technologies, as Illumina 
#reads are typically all the same length, hence the flat red line).

#We generally advise trimming the last few nucleotides to avoid less 
#well-controlled errors that can arise there. 

### VISUALIZING REVERSE READS ###

plotQualityProfile(fnRs[1:4])
plotQualityProfile(fnRs[5:8])

#Use these plots to determine how to set the filtering and trimming parameters.

# For this example data set, I am trimming the forward reads at 230 bp & reverse reads
# at 200 bp. These are the lengths at which average quality score is approximately
# 30.In addition, forward + reverse total length = 430 bp, which provides ample overlap
# greater than the length of the expected contigs (~380 - 400 bp).


#### FILTERING AND TRIMMING #### 


# Assign the filenames for the filtered fastq.gz files

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter forward and reverse reads

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE) 
head(out)


#### LEARN THE ERROR RATES #### 


# The program has to estimate the sequencing error rate to try and 
# distinguish errors from true Amplicon Sequence Variants.

errF <- learnErrors(filtFs, multithread=FALSE, verbose=TRUE)

errR <- learnErrors(filtRs, multithread=FALSE, verbose=TRUE)

### Visualize errors ###

plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

# If the machine learning algorithm is functioning as expected, you expect
# to see the log of the error frequency decreasing as consensus quality score increases.


#### SAMPLE INFERENCE #### 


# Infer the sequence variants in each sample #
# Reads are sorted into ASVs #

# ONLY CHOOSE ONE OF THE METHODS BELOW #

### SAMPLE-BY-SAMPLE INFERENCE (STANDARD) ###

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, verbose=TRUE)

### POOLED SAMPLE INFERENCE (TAKES LONGER, MIGHT CATCH RARE ASVs) ###

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, verbose=TRUE, pool=TRUE)

# There is also an intermediate method that approximates pooling via "pseudo-pooling."

### INSPECTING THE DADA-CLASS OBJECT RETURNED BY DADA ###

dadaFs[]
dadaRs[]


#### MERGE PAIRED READS #### 


### MERGE PAIRED-END READS ###

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

### INSPECT MERGER DATA FRAME (FIRST SAMPLE) ###

head(mergers[[1]])


#### CONSTRUCT SEQUENCE TABLE #### 

### CONSTRUCTING THE TABLE ###

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

### INSPECT SEQUENCE LENGTH DISTRIBUTION ###

table(nchar(getSequences(seqtab)))

# This is the number of ASVs of a given sequence length. For the Nubel amplicons,
# these sequences generally fall into two "bands" - one approximately 380 bp in length, and one approximately
# 400 bp in length. 
# Given the filter/trim parameters from earlier in the pipeline, the expected max contig length is 
# about 430 bp. The expected fragment length as seen in the Nubel paper and by gel imaging is approximately
# 420 bp in length. However, because we trimmed off the adapters and primers, we expect the contigs to be about
# 380 bp in length.  

#Lots of samples shorter/longer than expected? Could be non-specific priming.
#Non-target-length sequences can be removed in an "in silico" version of "cutting out a band"
#Use this code:

seqtab3 <- seqtab[,nchar(colnames(seqtab)) %in% 377:413]


#### REMOVE CHIMERIC SEQUENCES ####


### REMOVE THE CHIMERAS ###

seqtab3.nochim <- removeBimeraDenovo(seqtab3, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab3.nochim)

# It is not uncommon to lose a majority of your ASVs as chimeric. Most reads should remain, however.

### EXAMINING EFFECT OF CHIMERA REMOVAL ###
sum(seqtab3.nochim)/sum(seqtab3)

# The returned result is the proportion of NON-CHIMERIC reads. Per Dr. Callahan, 30% loss is the cutoff at which you should be concerned that something is 
# wrong with your data processing (residual primers/linkers?). 


#### TRACK READS THROUGH PIPELINE #### 


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab3.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track, c(8))

# Make sure there are no overly large drops between any of the steps. #


#### ASSIGN TAXONOMY #### 


### ASSIGN TAXONOMY ###

# Make sure to include PATH info for Silva database
# Note that this database is updated periodically
# Links to different versions are available with the DADA2 tutorial.
silva <- "/Users/hhuang22/Desktop/silva_nr99_v138.2_toSpecies_trainset.fa"
taxa <- assignTaxonomy(seqtab3.nochim, silva, multithread=TRUE)

#Using Silva v. 138.2 from https://zenodo.org/records/14169026 (thanks, Dr. Callahan!).
#Note that this database already assigns to species level

### DISPLAY OUTPUT ###

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#### SAVING & EXPORTING RESULTS ####


### SAVING ABUNDANCE TABLE AS A TXT FILE ###

write.table(t(seqtab3.nochim), "SFBay_Bioassay2023_Chimera_Free_Abundance.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### SAVING TAXONOMY AND ABUNDANCE TABLE AS A CSV FILE ###

write.csv(taxa, "SFBay_Bioassay_2023_Nubel_Taxonomy.csv")
write.csv(seqtab3.nochim, "SFBay_Bio_2023_Abundance.csv")

### SAVING ABUNDANCE AND TAXONOMY TABLES AS R OBJECTS AND AS CSV FILES ###

saveRDS(seqtab3.nochim, "SFBay_Bioassay_2023_Abundance.rds")
saveRDS(taxa, "SFBay_Bioassay_2023_Taxonomy.rds")


### ISOLATING SEQUENCE DATA AS A FASTA FILE ###

library(ShortRead)
uniquesToFasta(seqtab3.nochim, fout='SFBay_Bioassay_2023_Nubel.fna', ids=colnames(seqtab3.nochim))

asv_sfbio2023_seqs_16S <- colnames(seqtab3.nochim)

asv_sfbio2023_headers_16S <- vector(dim(seqtab3.nochim)[2], mode = "character")
for (i in 1:dim(seqtab3.nochim)[2]) {
  asv_sfbio2023_headers_16S[i] <- paste(">ASV", i)
}

asv_sfbio_2023_fasta_16S <- c(rbind(asv_sfbio2023_headers_16S, asv_sfbio2023_seqs_16S))
write(asv_sfbio_2023_fasta_16S, "SFBio2023_ASVs.fa")


#### END OF WORKFLOW ####
### Set the working directory ###

setwd("/Users/hhuang22/Desktop/SFBioassay")

### Clear workspace ###

rm(list=ls())

### Load required libraries ###

library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

### Assign path ###

path <- "/Users/hhuang22/Desktop/SFBioassay"
list.files(path)


#### READING IN REQUIRED DATA ####


### Import the ASV frequency table ###

ASVs_table <- read.table("/Users/hhuang22/Desktop/SFBioassay/SFBay_Bioassay2023_Chimera_Free_Abundance.txt")
ASVs_table1 <- as.data.frame(ASVs_table)
class(ASVs_table1)
head(ASVs_table1) # Reading in ASV frequency table & making a copy as a data frame

rownames(ASVs_table1) <- paste0("ASV", 1:nrow(ASVs_table1))
head(ASVs_table1) # Replacing ASV sequence data with ASV # (ASV1, ASV2, etc.)

ASVs_SFBio2023 <- otu_table(ASVs_table1, taxa_are_rows = TRUE) # Defining the "OTU table" object
# (term recognized by Phyloseq even though these are ASVs, not OTUs)

### Import the taxonomy table ###

Taxa_Cyanos <- read.csv("/Users/hhuang22/Desktop/SFBioassay/SFBay_Bioassay_2023_Nubel_Taxonomy.csv") 
# Reading in taxonomy table

Taxa_Cyanos1 <- data.frame(Taxa_Cyanos) # Making a copy of the taxonomy table as 
# a data frame

rownames(Taxa_Cyanos1) <- paste0("ASV", 1:nrow(Taxa_Cyanos1))
head(Taxa_Cyanos1) # Adding a column for ASV # (ASV1, ASV2, etc.)

Taxa_Cyanos2 = tax_table(as.matrix(Taxa_Cyanos1)) # Converting the taxonomy table into a matrix 
# and defining the resulting matrix as the "taxonomy table" object

#### CREATING A PHYLOSEQ OBJECT AND REMOVING NON-CYANOBACTERIAL ASVs/CHLOROPLASTS ####


### Creating a Phyloseq object ###

Metadata_SFBio2023 <- read.csv("/Users/hhuang22/Desktop/SFBioassay/MetaSFBayBio2023.csv", header=TRUE)
Metadata_SFBio2023_1 <- data.frame(Metadata_SFBio2023, row.names=1)
head(Metadata_SFBio2023_1)
Metadata_SFBio2023_1 = sample_data(Metadata_SFBio2023_1)

SFBayBio2023 <- phyloseq(ASVs_SFBio2023, Taxa_Cyanos2, Metadata_SFBio2023_1)

# Reference https://github.com/joey711/phyloseq/issues/541
str(SFBayBio2023) # Preview the phyloseq object. If you ran the "dummy variable" code, you should have
# n + 1 metadata categories.
SFBayBio2023 # Sanity check - do the OTU (ASV) table and taxonomy tables make sense?

### Gathering crucial information ###

Reads_by_sample <- colSums(otu_table(SFBayBio2023))
head(Reads_by_sample, c(79))

sort(Reads_by_sample, decreasing=TRUE)

SFBayBio2023_NoBlank = subset_samples(SFBayBio2023, sample_names(SFBayBio2023) != "SF079.BN")
SFBayBio2023_NoBlank

Reads_by_sample_Noblank <- colSums(otu_table(SFBayBio2023_NoBlank))
head(Reads_by_sample_Noblank)
min(Reads_by_sample_Noblank)
max(Reads_by_sample_Noblank)
mean(Reads_by_sample_Noblank)

#### RAREFYING READS W/ PHYLOSEQ ####

# Note: per the function used in Phyloseq: "Please note that the authors of phyloseq do not advocate using this
# as a normalization procedure, despite its recent popularity. Our justifications for using alternative approaches to address
# disparities in library sizes have been made available as an article in PLoS Computational Biology. See phyloseq_to_deseq2 for
# a recommended alternative to rarefying directly supported in the phyloseq package, as well as the supplemental materials for 
# the PLoS-CB article and the phyloseq extensions repository on GitHub. Nevertheless, for comparison and demonstration, the rarefying
# procedure is implemented here in good faith and with options we hope are useful."

# This topic is still very contentious in bioinformatics.

SFBayBio2023_phyloseq.rare <- rarefy_even_depth(SFBayBio2023_NoBlank, sample.size = 19001, rngseed = 413314413, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

# Note that default for "replace" is "TRUE," as this is faster and more computationally efficient, but it also means
# that the number of reads in a given ASV could exceed the original number of reads. Here, I decided to try setting it
# to "FALSE" so that the original count is the greatest value possible in any one ASV.
# The calculations ran instantaneously.

sample_sums(SFBayBio2023_phyloseq.rare)

### Removing uncharacterized ASVs at the kingdom/phylum levels ###

SFBayBio2023_1 <- subset_taxa(SFBayBio2023_phyloseq.rare, !is.na(Kingdom) & !Kingdom %in% c("", "Uncharacterized"))
SFBayBio2023_1 #Remove taxa with "N/A" at the Kingdom level
SFBayBio2023_2 <- subset_taxa(SFBayBio2023_1, !is.na(Phylum) & !Phylum %in% c("", "Uncharacterized"))
SFBayBio2023_2 #Remove taxa with "N/A" at the Phylum level

### Checking observed Kingdoms & Phyla ###
table(tax_table(SFBayBio2023_2)[, "Kingdom"])

table(tax_table(SFBayBio2023_2)[, "Phylum"])
### Saving the Phylum frequency table ###

SFBio2023_Phylum_Frequency_Table <- table(tax_table(SFBayBio2023_2)[, "Phylum"])
write.csv(SFBio2023_Phylum_Frequency_Table, "SFBio2023_Phylum_Frequency_Table.csv", row.names = TRUE, quote = FALSE)
### Removing Non-Cyanobacteriota taxa ###

SFBayBio2023_3 <- subset_taxa(SFBayBio2023_2, !is.na(Phylum) & Phylum %in% "Cyanobacteriota")
SFBayBio2023_3

#### SAVING THE "PRUNED" ASV & TAXONOMY TABLES ####


### Saving the pruned ASV frequency table ###

SFBayBio2023_ASVs_pruned <- as(otu_table(SFBayBio2023_3), "matrix")  # Generate new ASV frequency table
write.table(SFBayBio2023_ASVs_pruned, "SFBayBio2023_ASV_Frequency_Pruned.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# Save new ASV frequency table as a txt file

### Saving the pruned taxonomy table ###

SFBayBio2023_taxa_pruned <- as(tax_table(SFBayBio2023_3), "matrix")
write.table(SFBayBio2023_taxa_pruned, "SFBayBio2023_taxa_pruned.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### Saving the pruned taxonomy table as a CSV file ###

write.csv(SFBayBio2023_taxa_pruned, "SFBayBio2023_taxa_pruned.csv")
### Set the working directory ###

setwd("/Users/hhuang22/Desktop/SFBioassay/")

### Clear workspace ###

rm(list=ls())

### Load required libraries ###

library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

### Gathering crucial information ###

Reads_by_sample <- colSums(otu_table(SFBayBio2023_3))
head(Reads_by_sample, c(78))

order(Reads_by_sample)

Reads_by_sample_SFBayBio2023_3 <- colSums(otu_table(SFBayBio2023_3))
head(Reads_by_sample_SFBayBio2023_3)
min(Reads_by_sample_SFBayBio2023_3)
max(Reads_by_sample_SFBayBio2023_3)
mean(Reads_by_sample_SFBayBio2023_3)

### Constructing a simple bar plot w/ R base package ###

barplot(Reads_by_sample_SFBayBio2023_3, height = Reads_by_sample, main = "Cyanobacterial Reads per Sample", xlab = "SF Bay Bioassay Sample",
        ylab = "Total Number of Reads", ylim = c(0, 20000), col = c("darkred", "darkred", "tan4", "tan4", "seagreen", "seagreen",
                                                                    "blue", "blue"))

#### (OPTIONAL) GENERATING ASV ABUNDANCE TABLE FROM RAREFIED SAMPLES ####


write.csv(otu_table(SFBayBio2023_3), "Rarefied_absolute_abundance_SFBayBio2023_otus.csv")
write.csv(tax_table(SFBayBio2023_3), "Rarefied_absolute_abundance_SFBayBio2023_tax.csv")

#### CREATING AN NMDS PLOT ####
library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2)
packageVersion("ggplot2")
library(dplyr)
packageVersion("dplyr")
library(tidyverse)
packageVersion("tidyverse")

phyloseq_proportions_SFBayBio2023_3 <- transform_sample_counts(SFBayBio2023_3, function(otu) otu/sum(otu))
set.seed(413314413)
Ord.nmds.bray <- ordinate(phyloseq_proportions_SFBayBio2023_3, method="NMDS", distance="bray")

NMDS_plot_SFBayBio2023_3 <- plot_ordination(phyloseq_proportions_SFBayBio2023_3, Ord.nmds.bray,
                             color = "Treatment", 
                             title = "Bray NMDS: SF Bay Bioassay 2023")
NMDS_plot_SFBayBio2023_3 + scale_color_brewer(type = "qual", palette = "Set1") + geom_point(size = 2)


#### ANALYZING ALPHA DIVERSITY IN SAMPLES/TREATMENTS ####


### Reporting richness estimates for each sample ###

Richness_estimates_SFBayBio2023 <- estimate_richness(SFBayBio2023_3, split = TRUE, measures = NULL)
write.csv(Richness_estimates_SFBayBio2023, "Richness_estimates_by_sample_SFBayBio2023.csv")

### Plotting richness estimates for each sample ###

plot_richness(SFBayBio2023_3, title = "Alpha Diversity of Cyanobacterial ASVs", nrow = 2)
# Plot has all supported metrics

plot_richness(SFBayBio2023_3, measures = c("Observed", "Shannon", "Simpson"),
              title = "Alpha Diversity of Cyanobacterial ASVs")
# Plot w/ just Shannon, Simpson, and n observed ASVs 

### Plotting richness estimates for each treatment ###

plot_richness(SFBayBio2023_3, x = "Treatment",
              title = "Alpha Diversity of Cyanobacterial ASVs", nrow = 2)

# Plotting richness estimates w/ selected metrics

plot_richness(SFBayBio2023_3, x = "Treatment", measures = c("Observed", "Shannon", "Simpson"),
              title = "Alpha Diversity of Cyanobacterial ASVs")

plot_richness(SFBayBio2023_3, x = "Treatment", color = "Site", 
              title = "Alpha Diversity of Cyanobacterial ASVs", nrow = 2)


#### GENERATING RELATIVE ABUNDANCE PLOTS FOR DISCOVERY BAY BIOASSAY SAMPLES ####
library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2)
packageVersion("ggplot2")
library(dplyr)
packageVersion("dplyr")
library(tidyverse)
packageVersion("tidyverse")

###Creating phyloseq object for 2023 Bioassay with samples collected from Discovery Bay###

SFBio.taxa <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBioTax.csv")
SFBioDB.abund <- read.csv("C:/Users/hira7/OneDrive/Desktop/DBbioabund.csv")
SFBioMetaDB <- read.csv("C:/Users/hira7/OneDrive/Desktop/DBbiometa.csv")

SFBioDB.abund <- SFBioDB.abund %>% 
  tibble::column_to_rownames("ASV")
SFBioST.abund <- SFBioST.abund %>%
  tibble::column_to_rownames("ASV")
SFBio.taxa <- SFBio.taxa %>%
  tibble::column_to_rownames("ASV")
SFBioMetaDB <- SFBioMetaDB %>%
  tibble::column_to_rownames("Samples")
SFBioMetaST <- SFBioMetaST %>%
  tibble::column_to_rownames("Samples")
SFBioDB.abund <- as.matrix(SFBioDB.abund)
SFBioST.abund <- as.matrix(SFBioST.abund)
SFBio.taxa <- as.matrix(SFBio.taxa)
BioOTU.DB = otu_table(SFBioDB.abund, taxa_are_rows=FALSE)

BioOTU.ST = otu_table(SFBioST.abund, taxa_are_rows=FALSE)
BioTAX = tax_table(SFBio.taxa)
BioSamples.DB = sample_data(SFBioMetaDB)
BioSamples.ST = sample_data(SFBioMetaST)

psBioDB <- phyloseq(BioOTU.DB, BioTAX, BioSamples.DB)

# Convert tax table to data frame
BioTAX <- as.data.frame(tax_table(psBioDB), stringsAsFactors = FALSE)

# Identify target ASVs
chloroplast_asvs <- which(BioTAX$Order == "Chloroplast" & is.na(BioTAX$Genus))
length(chloroplast_asvs)

chloroplast_asvs1 <- which(BioTAX$Order == "Chloroplast" & BioTAX$Genus == "Chloroplast")
length(chloroplast_asvs1)

# Update taxonomy
BioTAX$Genus[chloroplast_asvs] <- "Chloroplast"

chloroplast_asvs2 <- which(BioTAX$Order == "Chloroplast" & is.na(BioTAX$Genus))
chloroplast_asvs2

chloroplast_asvs3 <- which(BioTAX$Order == "Chloroplast" & BioTAX$Genus == "Chloroplast")
length(chloroplast_asvs3)

BioTAX$Genus <- str_replace(BioTAX$Genus, "\\s.+$", "") #Keeps only the first word
BioTAX$Genus <- str_replace(BioTAX$Genus, "_.+$", "")  # For "_" separators

# Restore tax table (convert back to matrix)
BioTAX <- as.matrix(BioTAX)
BioTAX1 = tax_table(BioTAX)

# Form the phyloseq object again
psBioDB1 <- phyloseq(BioOTU.DB, BioTAX1, BioSamples.DB)
psBioDB1_cyano <- subset_taxa(psBioDB1, Phylum %in% c("Cyanobacteriota"))
taxa_names(psBioDB1_cyano) <- paste("ASV", 1:ntaxa(psBioDB1_cyano))

# Creating a bar plot based on the top 100 most abundant ASVs from Discovery Bay.

psBioDB1.cyano.rel <- transform_sample_counts(psBioDB1_cyano, function(OTU) OTU/sum(OTU))
top100.BioDB1.cyano <- names(sort(taxa_sums(psBioDB1_cyano), decreasing=TRUE))[1:100]
ps.top100.BioDB1.cyano <- prune_taxa(top100.BioDB1.cyano, psBioDB1.cyano.rel)

ps.top100.BioDB1.cyano

plot.top100.BioDB1.cyano <- plot_bar(ps.top100.BioDB1.cyano, fill="Genus")+theme_gray(base_size=18)+theme(axis.text.x=element_text(angle=270))
plot.top100.BioDB1.cyano

ps.top100.BioDB1.cyano.glom <- tax_glom(ps.top100.BioDB1.cyano, "Genus", NArm=FALSE)

SFDBBioDB1melt <- psmelt(ps.top100.BioDB1.cyano.glom)

plot.top100.BioDB1.cyano.glom <- ggplot(SFDBBioDB1melt)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus), stat="identity", position="stack")+
  theme_gray(base_size=20)+
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))+
  scale_fill_manual(values=c("red", "cornflowerblue","darkolivegreen3", "blue", "orange", "burlywood1", "darkgreen", "deepskyblue", "darkorchid", "cyan", "cornsilk2", "grey"))

plot.top100.BioDB1.cyano.glom

# Creating a bar plot, but this time, group the plots with the number of days post-treatment
# as the x-axis, and have the treatments at the same time point side by side.

ps_SF_top100_BioDB1_cyano_glom_day <- subset_samples(ps.top100.BioDB1.cyano.glom, Day %in% c("1", "3"))

View(ps_SF_top100_BioDB1_cyano_glom_day_avg)

ps_SF_top100_BioDB1_cyano_glom_day_avg <- psmelt(ps_SF_top100_BioDB1_cyano_glom_day) %>%
  group_by(Genus, Treatment, Day) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    .groups='drop'
  ) %>%
  mutate(Day = factor(Day, levels = c("1", "3"))) 

dodge <- position_dodge(width = 0.6)

custom_labels <- function(x) {
  case_when(
    x >= 0.01 ~ format(x, digits = 2), #Regular number (e.g., "0.05")
    x < 0.01 ~ scientific_format(digits = 2)(x), #Scientific if <0.01 (e.g., "1.0e-3")
  )
}

ggplot(ps_SF_top100_BioDB1_cyano_glom_day_avg, aes(x = Day, y = mean_abundance, fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.6, position = dodge) +
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance, 
                    ymax = mean_abundance + sd_abundance),
                width = 0.2,
                linewidth = 0.5,
                position = dodge
  ) +
  facet_wrap(~Genus, scales = "free", ncol = 3, labeller = labeller(Genus = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = custom_labels, breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Days Post-Treatment",
       y = "Relative Abundance",
       title = "Discovery Bay") +
  theme_minimal(base_size=16) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 1, margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold",margin = margin(t = 5, b = 5, l = 5, r = 5)),  # Padding
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        panel.spacing = unit(1.5, "lines"),  # More space between facets
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.8),
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))  # Outer margin

ggsave("SFBay_DB1_Bioassay2023_Top100_Cyano_Genus.png", width = 12, height = 14, units = "in", dpi = 600)

#### GENERATING RELATIVE ABUNDANCE PLOTS FOR STOCKTON BIOASSAY SAMPLES ####

SFBio.taxa <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBioTax.csv")
SFBioST.abund <- read.csv("C:/Users/hira7/OneDrive/Desktop/STbioAbund.csv")
SFBioMetaST <- read.csv("C:/Users/hira7/OneDrive/Desktop/STbiometa.csv")

SFBioST.abund <- SFBioST.abund %>%
  tibble::column_to_rownames("ASV")
SFBio.taxa <- SFBio.taxa %>%
  tibble::column_to_rownames("ASV")
SFBioMetaST <- SFBioMetaST %>%
  tibble::column_to_rownames("Samples")

SFBioST.abund <- as.matrix(SFBioST.abund)
SFBio.taxa <- as.matrix(SFBio.taxa)

BioOTU.ST = otu_table(SFBioST.abund, taxa_are_rows=FALSE)
BioTAX = tax_table(SFBio.taxa)
BioSamples.ST = sample_data(SFBioMetaST)

psBioST <- phyloseq(BioOTU.ST, BioTAX, BioSamples.ST)

# Convert tax table to data frame
BioTAX <- as.data.frame(tax_table(psBioDB), stringsAsFactors = FALSE)

### Modify the tax table so that chloroplast can be visualized at the genus level when
### making relative abundance plots

# Identify target ASVs
chloroplast_asvs <- which(BioTAX$Order == "Chloroplast" & is.na(BioTAX$Genus))
length(chloroplast_asvs)

chloroplast_asvs1 <- which(BioTAX$Order == "Chloroplast" & BioTAX$Genus == "Chloroplast")
length(chloroplast_asvs1)

# Update taxonomy
BioTAX$Genus[chloroplast_asvs] <- "Chloroplast"

chloroplast_asvs2 <- which(BioTAX$Order == "Chloroplast" & is.na(BioTAX$Genus))
chloroplast_asvs2

chloroplast_asvs3 <- which(BioTAX$Order == "Chloroplast" & BioTAX$Genus == "Chloroplast")
length(chloroplast_asvs3)

BioTAX$Genus <- str_replace(BioTAX$Genus, "\\s.+$", "") #Keeps only the first word
BioTAX$Genus <- str_replace(BioTAX$Genus, "_.+$", "")  # For "_" separators

# Restore tax table (convert back to matrix)
BioTAX <- as.matrix(BioTAX)
BioTAX1 = tax_table(BioTAX)

# Form the phyloseq object again
psBioST1 <- phyloseq(BioOTU.ST, BioTAX1, BioSamples.ST)
psBioST1_cyano <- subset_taxa(psBioST1, Phylum %in% c("Cyanobacteriota"))
taxa_names(psBioST1_cyano) <- paste("ASV", 1:ntaxa(psBioST1_cyano))

psBioST1.cyano.rel <- transform_sample_counts(psBioST1_cyano, function(OTU) OTU/sum(OTU))

# Creating a bar plot based on the top 100 most abundant ASVs from Stockton

top100.BioST1.cyano <- names(sort(taxa_sums(psBioST1_cyano), decreasing=TRUE))[1:100]
ps.top100.BioST1.cyano <- prune_taxa(top100.BioST1.cyano, psBioST1.cyano.rel)

ps.top100.BioST1.cyano

plot.top100.BioST1.cyano <- plot_bar(ps.top100.BioST1.cyano, fill="Genus")+theme_gray(base_size=18)+theme(axis.text.x=element_text(angle=270))
plot.top100.BioST1.cyano

ps.top100.BioST1.cyano.glom <- tax_glom(ps.top100.BioST1.cyano, "Genus", NArm=FALSE)

SFBioST1melt <- psmelt(ps.top100.BioST1.cyano.glom)

plot.top100.BioST1.cyano.glom <- ggplot(SFBioST1melt)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus), stat="identity", position="stack")+
  theme_gray(base_size=20)+
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))+
  scale_fill_manual(values=c("aquamarine", "cornflowerblue","darkolivegreen3", "blue", "orange", "chartreuse", "burlywood1", "darkgreen", "darkgoldenrod1", "darkorchid", "cyan", "cornsilk2", "grey"))

plot.top100.BioST1.cyano.glom

# Creating a bar plot, but this time, group the plots with the number of days post-treatment
# as the x-axis, and have the treatments at the same time point side by side.

ps_SF_top100_BioST1_cyano_glom_day <- subset_samples(ps.top100.BioST1.cyano.glom, Day %in% c("1", "3"))

ps_SF_top100_BioST1_cyano_glom_day_avg <- psmelt(ps_SF_top100_BioST1_cyano_glom_day) %>%
  group_by(Genus, Treatment, Day) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    .groups='drop'
  ) %>%
  mutate(Day = factor(Day, levels = c("1", "3"))) 

dodge <- position_dodge(width = 0.6)

custom_labels <- function(x) {
  case_when(
    x >= 0.01 ~ format(x, digits = 2), #Regular number (e.g., "0.05")
    x < 0.01 ~ scientific_format(digits = 2)(x), #Scientific if <0.01 (e.g., "1.0e-3")
  )
}

ggplot(ps_SF_top100_BioST1_cyano_glom_day_avg, aes(x = Day, y = mean_abundance, fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.6, position = dodge) +
  geom_errorbar(aes(ymin = pmax(mean_abundance - sd_abundance, 0),
                    ymax = mean_abundance + sd_abundance),
                width = 0.2,
                linewidth = 0.5,
                position = dodge
  ) +
  facet_wrap(~Genus, scales = "free", ncol = 3, labeller = labeller(Genus = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = custom_labels, breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Day Post-Treatment",
       y = "Relative Abundance",
       title = "Stockton") +
  theme_minimal(base_size=16) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 1, margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold",margin = margin(t = 5, b = 5, l = 5, r = 5)),  # Padding
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        panel.spacing = unit(1.5, "lines"),  # More space between facets
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.8),
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"))  # Outer margin

ggsave("SFBay_ST1_Bioassay2023_Top100_Cyano_Genus.png", width = 12, height = 14, units = "in", dpi = 600)


