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
library(tibble)
library(scales)
library(DESeq2)
library(stringr)

```
## Route of the directory
```{r}
pathNUBELSFTime <-"/Users/hhuang22/Desktop/SFTime"
```
## Checking the files
```{r}
list.files(pathNUBELSFTime)
```
## Checking the files
```{r}
fnFs <- sort(list.files(pathNUBELSFTime, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathNUBELSFTime, pattern = "_R2_001.fastq.gz", full.names = TRUE))

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

fnFs.filtN <- file.path(pathNUBELSFTime, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(pathNUBELSFTime, "filtN", basename(fnRs))
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

setwd("/Users/hhuang22/Desktop/SFTime/cutadapt")
rm(list=ls())

path <- "/Users/hhuang22/Desktop/SFTime/cutadapt"
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

 
setwd("/Users/hhuang22/Desktop/SFTime/cutadapt/Working_Reads")


### CLEAR WORKSPACE ###

rm(list=ls())

### ASSIGN PATH ###

path <- "/Users/hhuang22/Desktop/SFTime/cutadapt/Working_Reads"

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

setwd("/Users/hhuang22/Desktop/SFTime/cutadapt/Working_Reads")

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

path <- "/Users/hhuang22/Desktop/SFTime/cutadapt/Working_Reads"

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

seqtab3 <- seqtab[,nchar(colnames(seqtab)) %in% 377:404]


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

write.table(t(seqtab3.nochim), "SFBay_Time_Series_Chimera_Free_Abundance.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

### SAVING TAXONOMY AND ABUNDANCE TABLE AS A CSV FILE ###

write.csv(taxa, "SFBay_Time_Series_Nubel_Taxonomy.csv")
write.csv(seqtab3.nochim, "SFBay_Time_Series_Abundance.csv")

### SAVING ABUNDANCE AND TAXONOMY TABLES AS R OBJECTS AND AS CSV FILES ###

saveRDS(seqtab3.nochim, "SFBay_Time_Series_Abundance.rds")
saveRDS(taxa, "SFBay_Time_Series_Taxonomy.rds")


### ISOLATING SEQUENCE DATA AS A FASTA FILE ###

library(ShortRead)
uniquesToFasta(seqtab3.nochim, fout='SFBay_Time_Series_Nubel.fna', ids=colnames(seqtab3.nochim))

asv_sftime_seqs_16S <- colnames(seqtab3.nochim)

asv_sftime_headers_16S <- vector(dim(seqtab3.nochim)[2], mode = "character")
for (i in 1:dim(seqtab3.nochim)[2]) {
  asv_sftime_headers_16S[i] <- paste(">ASV", i)
}

asv_sftime_fasta_16S <- c(rbind(asv_sftime_headers_16S, asv_sftime_seqs_16S))
write(asv_sftime_fasta_16S, "SFTimeSeries_ASVs.fa")


#### END OF WORKFLOW ####
### Set the working directory ###

setwd("/Users/hhuang22/Desktop/SFTime")

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

path <- "/Users/hhuang22/Desktop/SFTime"
list.files(path)


#### READING IN REQUIRED DATA ####


### Import the ASV frequency table ###

ASVs_table <- read.table("/Users/hhuang22/Desktop/SFTime/SFBay_Time_Series_Chimera_Free_Abundance.txt")
ASVs_table1 <- as.data.frame(ASVs_table)
class(ASVs_table1)
head(ASVs_table1) # Reading in ASV frequency table & making a copy as a data frame

rownames(ASVs_table1) <- paste0("ASV", 1:nrow(ASVs_table1))
head(ASVs_table1) # Replacing ASV sequence data with ASV # (ASV1, ASV2, etc.)

ASVs_SFTime <- otu_table(ASVs_table1, taxa_are_rows = TRUE) # Defining the "OTU table" object
# (term recognized by Phyloseq even though these are ASVs, not OTUs)

### Import the taxonomy table ###

Taxa_Cyanos <- read.csv("/Users/hhuang22/Desktop/SFTime/SFBay_Time_Series_Nubel_Taxonomy.csv") 
# Reading in taxonomy table

Taxa_Cyanos1 <- data.frame(Taxa_Cyanos) # Making a copy of the taxonomy table as 
# a data frame

rownames(Taxa_Cyanos1) <- paste0("ASV", 1:nrow(Taxa_Cyanos1))
head(Taxa_Cyanos1) # Adding a column for ASV # (ASV1, ASV2, etc.)

Taxa_Cyanos2 = tax_table(as.matrix(Taxa_Cyanos1)) # Converting the taxonomy table into a matrix 
# and defining the resulting matrix as the "taxonomy table" object

#### CREATING A PHYLOSEQ OBJECT AND REMOVING NON-CYANOBACTERIAL ASVs ####

### Creating phyloseq object for cyanobacterial 16S PCR of SF Bay Time Series###

### The tables were altered outside of R such that the files related to Disocvery Bay and Stockton
### are being processed separately before being reuploaded to R to make phyloseq objects.

### We will start with creating the phyloseq object and relative abundance plots
### for Discovery Bay samples.

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tibble)
library(scales)
library(DESeq2)
library(stringr)

SFBay.taxa <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayTimeAbunDB.csv")
SFBayDB.abund <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayAbundDB.csv")
SFBayMetaDB <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayTimeMetaDB.csv")

SFBayDB.abund <- SFBayDB.abund %>% 
  tibble::column_to_rownames("ASV")
SFBay.taxa <- SFBay.taxa %>%
  tibble::column_to_rownames("ASV")
SFBayMetaDB <- SFBayMetaDB %>%
  tibble::column_to_rownames("Samples")

SFBayDB.abund <- as.matrix(SFBayDB.abund)
SFBay.taxa <- as.matrix(SFBay.taxa)
OTU.DB = otu_table(SFBayDB.abund, taxa_are_rows=FALSE)

Samples.DB = sample_data(SFBayMetaDB)
SFBay.taxa <- as.matrix(SFBay.taxa)
TAX = tax_table(SFBay.taxa)

psDB <- phyloseq(OTU.DB, TAX, Samples.DB)


### Identify ASVs with Order == "Chloroplast" AND Genus == NA,
### so that we can replace those NA under Genus with "Chloroplasts 
### for these ASVs because we will be doing a lot of the analysis 
### at the Genus level, and we might want to know which ones might actually 
### be eukaryotes

# Convert tax table to data frame

TAX <- as.data.frame(tax_table(psBioDB), stringsAsFactors = FALSE)

# Identify target ASVs
chloroplast_asvs <- which(TAX$Order == "Chloroplast" & is.na(TAX$Genus))
length(chloroplast_asvs)

chloroplast_asvs1 <- which(TAX$Order == "Chloroplast" & TAX$Genus == "Chloroplast")
length(chloroplast_asvs1)

# Update taxonomy
TAX$Genus[chloroplast_asvs] <- "Chloroplast"

chloroplast_asvs2 <- which(TAX$Order == "Chloroplast" & is.na(TAX$Genus))
chloroplast_asvs2

chloroplast_asvs3 <- which(TAX$Order == "Chloroplast" & TAX$Genus == "Chloroplast")
length(chloroplast_asvs3)

TAX$Genus <- str_replace(TAX$Genus, "\\s.+$", "") #Keeps only the first word
TAX$Genus <- str_replace(TAX$Genus, "_.+$", "")  # For "_" separators

# Restore tax table (convert back to matrix)
TAX <- as.matrix(TAX)
TAX1 = tax_table(TAX)

# Form the phyloseq object again
psDB1 <- phyloseq(OTU.DB, TAX1, SFBayMetaDB)
psDB1_cyano <- subset_taxa(psDB1, Phylum %in% c("Cyanobacteriota"))

taxa_names(psDB1_cyano) <- paste("ASV", 1:ntaxa(psDB1_cyano))

psDB1.cyano.rel <- transform_sample_counts(psDB1_cyano, function(OTU) OTU/sum(OTU))

### Create bar plot based on the top 100 most abundant ASVs from Discovery Bay samples

top100.DB1.cyano <- names(sort(taxa_sums(psDB1_cyano), decreasing=TRUE))[1:100]
ps.top100.DB1.cyano <- prune_taxa(top100.DB1.cyano, psDB1.cyano.rel)

ps.top100.DB1.cyano

plot.top100.DB1.cyano <- plot_bar(ps.top100.DB.cyano, x="Date", fill="Genus")+theme_gray(base_size=18)+theme(axis.text.x=element_text(angle=270))
plot.top100.DB1.cyano

ps.top100.DB1.cyano.glom <- tax_glom(ps.top100.DB.cyano, "Genus", NArm=FALSE)

### Creating a bar plot in which the information has collapsed into the genus
 
SFBayDB1melt <- psmelt(ps.top100.DB1.cyano.glom)

plot.top100.DB1.cyano.glom <- ggplot(SFBayDBmelt)+
  geom_bar(aes(x=Date, y=Abundance, fill=Genus), stat="identity", position="stack")+
  theme_gray(base_size=20)+
  theme(axis.text.x=element_text(angle=270))+
  scale_fill_manual(values=c("red", "darkolivegreen3", "blue", "orange", "cyan", "darkgoldenrod1", "darkseagreen", "deeppink", "purple", "cornsilk2", "deepskyblue", "coral", "darkgreen"))

plot.top100.DB1.cyano.glom

### Processing data from Stockton

SFBay.taxa <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayTimeAbunDB.csv")
SFBayST.abund1 <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayAbundST.csv")
SFBayMetaST <- read.csv("C:/Users/hira7/OneDrive/Desktop/SFBayTimeMetaST.csv")

SFBayST.abund <- SFBayST.abund %>%
  tibble::column_to_rownames("ASV")
SFBay.taxa <- SFBay.taxa %>%
  tibble::column_to_rownames("ASV")
SFBayMetaST <- SFBayMetaST %>%
  tibble::column_to_rownames("Samples")

SFBayST.abund <- as.matrix(SFBayST.abund)
SFBay.taxa <- as.matrix(SFBay.taxa)

OTU.ST = otu_table(SFBayST.abund, taxa_are_rows=FALSE)
Samples.ST = sample_data(SFBayMetaST)
TAX = tax_table(SFBay.taxa)


psST <- phyloseq(OTU.ST, TAX, Samples.ST)

# Moving the chloroplast from "Order" to "Genus", so that we can visualize how many
# ASVs are eukaryotic at the genus level.

# Convert tax table to data frame

TAX <- as.data.frame(tax_table(psST), stringsAsFactors = FALSE)

# Identify target ASVs
chloroplast_asvs <- which(TAX$Order == "Chloroplast" & is.na(TAX$Genus))
length(chloroplast_asvs)

chloroplast_asvs1 <- which(TAX$Order == "Chloroplast" & TAX$Genus == "Chloroplast")
length(chloroplast_asvs1)

# Update taxonomy
TAX$Genus[chloroplast_asvs] <- "Chloroplast"

chloroplast_asvs2 <- which(TAX$Order == "Chloroplast" & is.na(TAX$Genus))
chloroplast_asvs2

chloroplast_asvs3 <- which(TAX$Order == "Chloroplast" & TAX$Genus == "Chloroplast")
length(chloroplast_asvs3)

TAX$Genus <- str_replace(TAX$Genus, "\\s.+$", "") #Keeps only the first word
TAX$Genus <- str_replace(TAX$Genus, "_.+$", "")  # For "_" separators

# Restore tax table (convert back to matrix)
TAX <- as.matrix(TAX)
TAX1 = tax_table(TAX)

# Form the phyloseq object again
psST1 <- phyloseq(OTU.ST, TAX1, SFBayMetaST)
psST1_cyano <- subset_taxa(psST1, Phylum %in% c("Cyanobacteriota"))

taxa_names(psST1_cyano) <- paste("ASV", 1:ntaxa(psST1_cyano))

psST1.cyano.rel <- transform_sample_counts(psST1_cyano, function(OTU) OTU/sum(OTU))

top100.ST1.cyano <- names(sort(taxa_sums(psST1_cyano), decreasing=TRUE))[1:100]
ps.top100.ST1.cyano <- prune_taxa(top100.ST1.cyano, psST1.cyano.rel)
ps.top100.ST1.cyano

write.csv2(otu_table(ps.top100.ST1.cyano),"C://Users/hira7/OneDrive/Desktop/ST100cyano2.xlsx")

# Plotting the distribution of the top 100 ASVs from Stockton's samples.

plot.top100.ST1.cyano <- plot_bar(ps.top100.ST1.cyano, x="Date", fill="Genus")+theme_gray(base_size=18)+theme(axis.text.x=element_text(angle=270))
plot.top100.ST1.cyano

ps.top100.ST1.cyano.glom <- tax_glom(ps.top100.ST1.cyano, "Genus", NArm=FALSE)

SFBayST1melt <- psmelt(ps.top100.ST1.cyano.glom)
plot.top100.ST1.cyano.glom <- ggplot(SFBayST1melt)+
  geom_bar(aes(x=Date, y=Abundance, fill=Genus), stat="identity", position="stack")+
  theme_gray(base_size=20)+
  theme(axis.text.x=element_text(angle=270))+
  scale_fill_manual(values=c("darkolivegreen3", "blue", "cyan", "darkseagreen", "aquamarine", "deeppink", "purple", "cornsilk2","black"))

plot.top100.ST1.cyano.glom
