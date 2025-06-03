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
pathSFmcyA <-"/Users/hhuang22/Desktop/SFmcyA"
```
## Checking the files
```{r}
list.files(pathSFmcyA)
```
## Checking the files
```{r}
fnFs <- sort(list.files(pathSFmcyA, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathSFmcyA, pattern = "_R2_001.fastq.gz", full.names = TRUE))

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

fnFs.filtN <- file.path(pathSFmcyA, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(pathSFmcyA, "filtN", basename(fnRs))
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

setwd("/Users/hhuang22/Desktop/SFmcyA/cutadapt")
rm(list=ls())

path <- "/Users/hhuang22/Desktop/SFmcyA/cutadapt"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

### IDENTIFY PRIMERS ###

# I am using the mcyA primers.

FWD <- "AAAATTAAAAGCCGTATCAAA"  
REV <- "AAAAGTGTTTTATTAGCGGCTCAT" 

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


setwd("/Users/hhuang22/Desktop/SFTime/cutadapt/Working_Reads")


### CLEAR WORKSPACE ###

rm(list=ls())

### ASSIGN PATH ###

path <- "/Users/hhuang22/Desktop/SFmcyA/cutadapt/Working_Reads"

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

setwd("/Users/hhuang22/Desktop/SFmcyA/cutadapt/Working_Reads")

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

path <- "/Users/hhuang22/Desktop/SFmcyA/cutadapt/Working_Reads"

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

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=TRUE) 
head(out)


#### LEARN THE ERROR RATES #### 


# The program has to estimate the sequencing error rate to try and 
# distinguish errors from true Amplicon Sequence Variants.

errF <- learnErrors(filtFs, multithread=TRUE, verbose=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE, verbose=TRUE)

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

### POOLED SAMPLE INFERENCE (TAKES LONGER, MIGHT CATCH RARE ASVs. We are choosing "pool" for this one) ###

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

#Lots of samples shorter/longer than expected? Could be non-specific priming.
#Non-target-length sequences can be removed in an "in silico" version of "cutting out a band"
#Use this code:

seqtab3 <- seqtab[,nchar(colnames(seqtab)) %in% 240:425]


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

# Make sure to include PATH info for mcyA database
# Note that this database is updated periodically
# Links to different versions are available with the DADA2 tutorial.
mcyA <- "/Users/hhuang22/Desktop/mcyA.txt"
taxa <- assignTaxonomy(seqtab3.nochim, mcyA, multithread=TRUE)

### DISPLAY OUTPUT ###

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#### SAVING & EXPORTING RESULTS ####


### SAVING ABUNDANCE TABLE AS A TXT FILE ###

write.table(t(seqtab3.nochim), "SFBay_mcyA_Chimera_Free_Abundance.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### SAVING TAXONOMY AND ABUNDANCE TABLE AS A CSV FILE ###

write.csv(taxa, "SFBay_mcyA_Taxonomy.csv")
write.csv(seqtab3.nochim, "SFBay_mcyA_Abundance.csv")

### SAVING ABUNDANCE AND TAXONOMY TABLES AS R OBJECTS AND AS CSV FILES ###

saveRDS(seqtab3.nochim, "SFBay_mcyA_Abundance.rds")
saveRDS(taxa, "SFBay_mcyA_Taxonomy.rds")


### ISOLATING SEQUENCE DATA AS A FASTA FILE ###

library(ShortRead)
uniquesToFasta(seqtab3.nochim, fout='SFBay_mcyA.fna', ids=colnames(seqtab3.nochim))

asv_sf_mcyA <- colnames(seqtab3.nochim)

asv_sf_mcyA_headers <- vector(dim(seqtab3.nochim)[2], mode = "character")
for (i in 1:dim(seqtab3.nochim)[2]) {
  asv_sf_mcyA_headers[i] <- paste(">ASV", i)
}

asv_sf_mcyA <- c(rbind(asv_sf_mcyA_headers, asv_sf_mcyA))
write(asv_sf_mcyA, "SFmcyA_ASVs.fa")


#### END OF WORKFLOW ####
### Set the working directory ###

setwd("/Users/hhuang22/Desktop/SFmcyA")

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

path <- "/Users/hhuang22/Desktop/SFmcyA"
list.files(path)


#### READING IN REQUIRED DATA ####


### Import the ASV frequency table ###

ASVs_table <- read.table("/Users/hhuang22/Desktop/SFmcyA/SFBay_mcyA_Chimera_Free_Abundance.txt")
ASVs_table1 <- as.data.frame(ASVs_table)
class(ASVs_table1)
head(ASVs_table1) # Reading in ASV frequency table & making a copy as a data frame

rownames(ASVs_table1) <- paste0("ASV", 1:nrow(ASVs_table1))
head(ASVs_table1) # Replacing ASV sequence data with ASV # (ASV1, ASV2, etc.)

ASVs_SF_mcyA <- otu_table(ASVs_table1, taxa_are_rows = TRUE) # Defining the "OTU table" object
# (term recognized by Phyloseq even though these are ASVs, not OTUs)

### Import the taxonomy table ###

Taxa_Cyanos <- read.csv("/Users/hhuang22/Desktop/SFmcyA/SFBay_mcyA_Taxonomy.csv") 
# Reading in taxonomy table

Taxa_Cyanos1 <- data.frame(Taxa_Cyanos) # Making a copy of the taxonomy table as 
# a data frame

rownames(Taxa_Cyanos1) <- paste0("ASV", 1:nrow(Taxa_Cyanos1))
head(Taxa_Cyanos1) # Adding a column for ASV # (ASV1, ASV2, etc.)

Taxa_Cyanos2 = tax_table(as.matrix(Taxa_Cyanos1)) # Converting the taxonomy table into a matrix 
# and defining the resulting matrix as the "taxonomy table" object

#### CREATING A PHYLOSEQ OBJECT AND REMOVING NON-CYANOBACTERIAL ASVs ####


### Creating a Phyloseq object ###

Metadata_SF_mcyA <- read.csv("/Users/hhuang22/Desktop/SFmcyA/MetaSFBaymcyA.csv", header=TRUE)
Metadata_SF_mcyA_1 <- data.frame(Metadata_SF_mcyA, row.names=1)
head(Metadata_SF_mcyA_1)
Metadata_SF_mcyA_1 = sample_data(Metadata_SF_mcyA_1)

SFmcyA <- phyloseq(ASVs_SF_mcyA, Taxa_Cyanos2, Metadata_SF_mcyA_1)

# Reference https://github.com/joey711/phyloseq/issues/541
str(SFmcyA) # Preview the phyloseq object. If you ran the "dummy variable" code, you should have
# n + 1 metadata categories.
SFmcyA # Sanity check - do the OTU (ASV) table and taxonomy tables make sense?

### Removing uncharacterized ASVs at the kingdom/phylum levels ###

SFmcyA_1 <- subset_taxa(SFmcyA, !is.na(Kingdom) & !Kingdom %in% c("", "Uncharacterized"))
SFmcyA_1 #Remove taxa with "N/A" at the Kingdom level

### Checking observed Kingdoms, Phyla, and Genera ###
table(tax_table(SFmcyA_1)[, "Kingdom"])
table(tax_table(SFmcyA_1)[, "Phylum"])
table(tax_table(SFmcyA_1)[, "Genus"])

### Saving the Genus frequency table ###

SFmcyA_Genus_Frequency_Table <- table(tax_table(SFmcyA_1)[, "Genus"])
write.csv(SFmcyA_Genus_Frequency_Table, "SF_mcyA_Frequency_Table.csv", row.names = TRUE, quote = FALSE)

#### SAVING THE "PRUNED" ASV & TAXONOMY TABLES ####


### Saving the pruned ASV frequency table ###

SFmcyA_ASVs_pruned <- as(otu_table(SFmcyA_1), "matrix")  # Generate new ASV frequency table
write.table(SFmcyA_ASVs_pruned, "SFmcyA_ASV_Frequency_Pruned.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# Save new ASV frequency table as a txt file

### Saving the pruned taxonomy table ###

SFmcyA_taxa_pruned <- as(tax_table(SFmcyA_1), "matrix")
write.table(SFmcyA_taxa_pruned, "SFmcyA_taxa_pruned.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### Saving the pruned taxonomy table as a CSV file ###

write.csv(SFmcyA_taxa_pruned, "SFmcyA_taxa_pruned.csv")
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

#### GENERATING RELATIVE ABUNDANCE PLOTS ####
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

top20_SFmcyA_1 <- names(sort(taxa_sums(SFmcyA_1), decreasing = TRUE))[1:20]
# Create an object consisting of the top 20 most abundant taxa in decreasing order

TransSFmcyA_1 <- transform_sample_counts(SFmcyA_1, function(ASV_Table)ASV_Table/sum(ASV_Table))
# Transform sample counts of ALL taxa into relative abundances and save as a Phyloseq object

top20_SFmcyA_1 <- prune_taxa(top20_SFmcyA_1, TransSFmcyA_1)
# Prune this Phyloseq object to include only the top 20 most abundant taxa we defined above

plot_bar(top20_SFmcyA_1, fill = "Genus") + theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90))+facet_wrap(~Site, scales="free_x")+geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
# Generates the plot of relative abundance by SAMPLE at the GENUS level.
#### END OF WORKFLOW ####
