####CO1 Script#####

# SCRIPT prepared by Matthieu Leray
# April 22nd 2020
#
# SCRIPT modified by Zoé Chamot
# May 30th 2022
#
# Last modified by Luisa Meister
# Feburary 2025
# READ PROCESSING



# Steps followed

# (1) Filter raw reads & Infer Amplicon Sequence Variants (ASVs)
# (2) Assign taxonomy
# (3) Create Phyloseq object
# (4) Remove contaminants, control and outlier samples
# (5) Cluster ASVs in OTUs at 3%
# (6) OTU curation with LULU
# (7) Taxonomic assignments with BLAST
# (8) Create curated Phyloseq object


# Installing packages
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("speedyseq") # funktioniert nicht
BiocManager::install("decontam")
BiocManager::install("metagMisc")# funktioniert nicht
BiocManager::install("DECIPHER")
BiocManager::install("lulu") # funktioniert nicht
BiocManager::install("Biostrings")


# Install speedyseq, metagMisc, and lulu (install.packages or github)
devtools::install_github("vmikk/metagMisc")
devtools::install_github("tobiasgf/lulu")
remotes::install_github("mikemc/speedyseq")


#(0)  Load packages ####

library(BiocGenerics)
library(Biostrings)
library(XVector)
library(GenomeInfoDb)
library(usethis)
library(Rcpp)
library(dada2)
library(ggplot2)
library(phyloseq)
library(speedyseq)
library(tibble)
library(dplyr)
library(decontam) # Removal of contaminants
library(metagMisc)
library(DECIPHER) # Align sequences
library(lulu) # OTU curation
library(ape) #needed for phagorn
library(phangorn) # phylogenetic tree
library(ShortRead)
library(RcppParallel)


#_________________________________________________________________________
#
#### (1) Sequencing Data Processing ####
#_________________________________________________________________________

#creating file path to data
path <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/Fastq_ROHR27_rerun/unzip/trimmed"
list.files(path)

##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
#sample.names <- sapply(strsplit(fnFs, "-"), '[', 1)
sample.names <- tools::file_path_sans_ext(basename(fnFs))
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

################################################################
# Remove empty sample files
################################################################

# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.exists <- fnFs[file.size(fnFs) > 50]
length(fnFs.exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.exists <- fnRs[file.size(fnRs) > 50]
length(fnRs.exists)
file.size(fnFs.exists)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs <- fnFs.exists
fnRs <- fnRs.exists
length(fnFs)
length(fnRs)
file.size(fnFs)

# Update your samples names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
length(sample.names)

#plotting quality profiles
qprofile_fw <- print(plotQualityProfile(fnFs, aggregate = TRUE)
                      + ggtitle("Forward"))
# -> cut at 270 or 280
qprofile_rev <- print(plotQualityProfile(fnRs, aggregate = TRUE)
                       + ggtitle("Reverse"))

#placing filtered files in a new filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering and trimming, here truncation at 220 (Fwd) and 180 (Rev) bp,
#2expected errors max (N discarded automatically)

bp <- SnowParam(workers = 8)
register(bp)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 230),
                      maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = TRUE)
head(out)


#pathfiltered <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/ITS2/Fastq_ROHR31/cutadapt/filtered"
#list.files(pathfiltered)

##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
#filtFs <- sort(list.files (pathfiltered, pattern = "_F_filt.fastq"))
#filtRs <- sort(list.files (pathfiltered, pattern = "_R_filt.fastq"))
#sample.names <- sapply(strsplit(filtFs, "_"), '[', 1)
#fnFs <- file.path(path, fnFs)
#fnRs <- file.path(path, fnRs)
#to just consider files that pass the filtering
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs [exists]
filtRs <- filtRs [exists]

# Learning error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plotting errors
plotErrors(errF, nominal = FALSE)
plotErrors(errR, nominal = FALSE)

# Dereplicating reads
#sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sam.names <- tools::file_path_sans_ext(basename(filtFs))
derepFs <- derepFastq(filtFs)
names(derepFs) <- sam.names
derepRs <- derepFastq(filtRs)
names(derepRs) <- sam.names

# Save copies of dereplicated objects
saveRDS(derepFs, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/derepFs.rds")
saveRDS(derepRs, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/derepRs.rds")


#____________________________________________________________________________
#
#### (2) Infering Amplification Sequence Variants with dada ####
#____________________________________________________________________________

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaFs [[1]]
dadaFs [[5]]
dadaFs [[100]]
dadaFs [[200]]
#etc.
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
dadaRs [[1]]
dadaRs [[2]]
dadaRs [[200]]
#etc.

# Save copies of dada objects
saveRDS(dadaFs, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/dadaFs_CO1.rds")
saveRDS(dadaRs, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/dadaRs_CO1.rds")


##Merging paired ends
merge <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab <- makeSequenceTable(merge)
dim(seqtab)

#CO1 [1]   253 70614

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merge[[1]])
head(merge[[6]])


table(nchar(getSequences(seqtab)))

#CO1
# 270   271   272   273   274   275   276   277   278   279   280   281   282   283   284   285   286   287   288   289
#113     6     8    13    10     6    34    12     3     8    16     8     8     7     6     3    13     8    12    28
#290   291   292   293   294   295   296   297   298   299   300   301   302   303   304   305   306   307   308   309
#11     8    21    11     7    47    14    21   245    26    78   623    27    53  1643    55   105  1626    59    98
#310   311   312   313   314   315   316   317   318   319   320   321   322   323   324   325   326   327   328   329
#1709   537  1966 56149   197    86  1648    26    37   251    22    16    66    75     6    11     1     5    27     6
#330   331   332   333   334   335   336   337   338   339   340   341   342   343   344   346   347   348   349   350
#6    19     2     5    13     1     3    21    44     3    63     4     1     9     4     3    22    18   686    18
#351   352   353   354   355   356   357   358   359   363   364   366   368   373   374   375   376   379   380   381
#42  1417     9     4    59     1     4    46     3     3     9     1     2     3     3     2     1     2     1     1
#382   384   385   386   387   388   391   392   397   402   403   405   406   407   409   413   415   417   418   419
#9     1    10     1     1     6     1     1     1     2     1     1     2     2     2     1     3     1     1     1
#420   421   424   425   427   428   430   433   436   438   439   440   441   442   443   445   446   448   449   451
#1     2     9     1     2     1     1     1     1     1     4     2     4     4     1     1     2     3     1     1
#452   453   457   458   465   466   469   470   471   472   481   483   485   487
#1     1     1     1     1     4     1     7     1     1     1     1     1     2

# Remove non specific amplifications -> keep only sequences that are 313bp
Final_data <- seqtab[, nchar(colnames(seqtab)) %in% 310:316]
Final_data2seqs <- getSequences(Final_data)
writeFasta(Final_data2seqs, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data2seqs_new.fasta")

# Saving files to use in the next part of the workflow
saveRDS(Final_data, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data_new.rds")


# Identifying and removing chimeras # set multithread to 10 next time
seqtab.nochim <- removeBimeraDenovo(Final_data,
                                    method = "consensus",
                                    multithread = 10)

dim(seqtab.nochim)
#CO1 [1]   253 68618 new 253 60448

# Saving table without chimeras and downloading as .csv file

saveRDS(seqtab.nochim, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/seqtab_nochim_new.rds")
write.csv(seqtab.nochim, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/nochim_seq_new.csv")

# Generate .txt document with the changes generated in each step

getN <- function(x) sum(getUniques(x))
track <-    cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                  sapply(merge, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sam.names
write.table(track, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/tracked_COI_new.txt")
saveRDS(seqtab.nochim, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/seqtab.nochim_new.rds")

################################################################
#### (3) Assign taxonomy to create Phyloseq #####
###############################################################
# This is just a quick assignment to create a taxonomy table, and then a phyloseq oject that can be used to detect contaminants
# To speed things up, we use a version of MIDORI with a few sequences

#assigning taxonomy with MIDORI
#reference datasets formatted for DADA2 can be found here: http://www.reference-midori.info/download.php#
set.seed(119)

#prepare small Midori dataset

#library(ShortRead)
#Midoripath <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/MiFish/MIDORI2_UNIQ_NUC_GB261_srRNA_DADA2.fasta/MIDORI2_UNIQ_NUC_GB261_srRNA_DADA2.fasta"
#Midori_data <- read.FASTA(Midoripath)
#head(Midoripath)
#Midorisubset <- Midori_data[1:200]
#print(Midorisubset)
#write.FASTA(Midorisubset, file="C:/Users/luisa/Documents/eDNA/Data_Analysis/R/MiFish/MIDORI2_UNIQ_NUC_GB261_srRNA_DADA2.fasta/Midorisubset.FASTA")

#use small Midori dataset

Midori_path <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/MiFish/MIDORI2_UNIQ_NUC_GB261_srRNA_DADA2.fasta/Midorisubset.FASTA"
Midori_data <- read.FASTA(Midori_path)


#go on
Final_data.nochim.tax <- assignTaxonomy(seqtab.nochim,
                                        "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/MiFish/MIDORI2_UNIQ_NUC_GB261_srRNA_DADA2.fasta/Midorisubset.FASTA",
                                        multithread = 10)


saveRDS(Final_data.nochim.tax, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data.nochim.tax_new.rds")


# Inspect taxonomic assignments
taxa.print <- Final_data.nochim.tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_dada_tax <- tax_table(Final_data.nochim.tax)
write.csv(tax_dada_tax, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data.nochim.tax_new.csv")

tax_dada_otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
rownames(tax_dada_otu) <- sub("^(([^_]*_){2}[^_]*)_.*", "\\1", rownames(tax_dada_otu))

write.csv(tax_dada_otu, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data.nochim.tax.otu_new.csv")


################################################################
# (4) Create Phyloseq object ####
################################################################

# Read RDS object
Final_data.nochim.tax <- readRDS("C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data.nochim.tax_new.rds")
Final_data.nochim <- readRDS("C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data.nochim_new.rds")

dim(Final_data.nochim.tax)
#COI 68618     6

dim(seqtab.nochim)
#COI 253 68618 new 60448     6

##Creating phyloseq object
##Opening and extracting sample data from a .csv file
#creating path to the .csv
# metadatasheet has to have column with sample or control

Final_data_sam <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Metadatasheet_CO1_NEW1.csv"
#reading the sample data sheet
Final_data_metadata <- read.csv(Final_data_sam, header = TRUE, sep = ";", row.names = 1)
head(Final_data_metadata)

##Creating a unique phyloseq object with sample-by-sequence feature table,
# the sample metadata and the sequence taxonomies
Final_data_obj <- phyloseq(otu_table(tax_dada_otu, taxa_are_rows = FALSE),
                           sample_data(Final_data_metadata), tax_table(tax_dada_tax))

Final_data_obj
#CO1
#otu_table()   OTU Table:          [ 60448 taxa and 253 samples ]:
 # sample_data() Sample Data:        [ 253 samples by 19 sample variables ]:
  #tax_table()   Taxonomy Table:     [ 60448 taxa by 6 taxonomic ranks ]:
################################################################
# (5) Remove contaminants, control and outlier samples ####
################################################################

# Identify Contaminants - Prevalence ##### R package decontam
sample_data(Final_data_obj)$is.neg <- sample_data(Final_data_obj)$SAMPLE_OR_CONTROL == "Control"
contamdf.prev <- isContaminant(Final_data_obj, method = "prevalence", neg = "is.neg")
contaminants <- table(contamdf.prev$contaminant)
write.csv(contamdf.prev, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data_contamdf.prev_new.csv")
print(contaminants)
subset_contaminants <- contamdf.prev[contamdf.prev$contaminant == TRUE, ]
write.csv(subset_contaminants, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data_contaminants_only_new.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
Final_data_obj.pa <- transform_sample_counts(Final_data_obj, function(abund) 1 * (abund > 0))

Final_data_obj.pa.neg <- prune_samples(sample_data(Final_data_obj.pa)$SAMPLE_OR_CONTROL == "Control", Final_data_obj.pa)
Final_data_obj.pa.pos <- prune_samples(sample_data(Final_data_obj.pa)$SAMPLE_OR_CONTROL == "Sample", Final_data_obj.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(Final_data_obj.pa.pos), pa.neg = taxa_sums(Final_data_obj.pa.neg),
                    contaminant = contamdf.prev$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


#Co1 27.01.2025
badTaxa <- c("TCTAGCAGCAACGCTAGGTCATAGAGGCTCTTCTGTCGACTTGGCTATTTTCTCTTTACACCTAGCCGGAGCCTCATCTATCTTAGGAGCTGTTAACTTCATCTCAACCGTTATCAACATACGCTCAGAAGGTATGTATATAGACCGCATCCCTCTATTTGTTTGATCAATCTTTATCACAGCAATTTTATTACTACTATCTCTACCCGTGTTAGCCGGAGCGATCACTATATTACTTACAGACCGAAATTTAAACACCTCTTTCTTCGACCCTTCAGGAGGGGGAGATCCTATTCTATACCAACACCTGTTT",
"TCTAGCAAGTAATATTGCCCATGCAGGCCCTTCCGTAGATTTAGCAATCTTTTCTCTTCATCTAGCAGGTGTTTCATCTATCTTGGGTGCAATCAATTTTATATCAACAGTCATCAATATACGATCAAAAGGGCTACGGATAGAACGAGTCCCCTTATTCGTTTGGTCAGTGCTAGTAACTGCAACTCTACTCCTTTTATCCCTCCCTGTACTAGCCGCAGCAATTACCATATTGCTGACGGACCGAAACCTCAACACCACTTTCTTCGACCCAAGAGGGGGAGGTGACCCTATCCTTTACCAACACTTGTTC",
"ATTAGCCGCAAATATTTCCCACAGTGGACCATCAGTTGATCTAGCTATTTTTTCTCTACATTTAGCAGGAGTTTCCTCAATTCTAGGTTCTATTAACTTTATTACCACAGTACGCAATATACGTCCAAAAGGAATAACTCCTGAACGAATCCCTCTATTTGTATGAGCAGTAGTACTAACAACCATTTTATTACTACTATCTCTACCAGTCCTAGCAGGAGGAATTACAATATTATTAACTGATCGTAATCTAAATACATCATTCTTCGATCCAGGAGGAGGAGGAGACCCAATTCTATACCAACATTTATTT",
"TTATCCGGAAGGACTGCACACAGAGGCCCCCCAGTTGATATAGCTATCTTTTCTCTTCATTTGGCCGGTGCAAGGTCTATTTTAGGAGCAGCAAATTTTATTGCTACACTAATTAACATACGACTTAAAAAGCAGAGGTCAGACAAAATACCTCTATTCCCATGATCTGTGATTCTAACCGCCATTCTTCTACTTCTCTCTCTACCTGTTTTAGCTGGGGCCATCACTATGTTATTGACAGATCGTAACCTAAATACAACTTTCTTTGATAGTGGTGGAGGGGGAGATCCAATTCTATATCAGCACCTTTTC",
"TCTATCTAATAACACTGCCCACTCTGGGTTTTCTGTTGACCTAACCATCTTCTCCCTCCATTTAGCTGGGGTCAGCTCTATTTTAGGCGCAATCAATTTCATAACCACTATTCTCAACATACGAACATCGGGAATAACCCTAGAACGAATACCCCTATTTGTTTGATCTATTTTCATTACTACCCTTCTACTACTCTTATCCCTACCAGTGTTGGCAGGAGCAATTACCATATTATTAACAGATCGCAATTTTAATACTAGTTTTTTTGACCCTACAGGAGGGGGAGACCCTATTTTATATCAACACCTATTC",
"ATTATCAAGCATACAATCACACTCAGGTGGTTCTGTTGATCTTGCAATTTTCAGTTTACACTTAGCTGGAGTATCATCCTTATTAGGTGCTATTAATTTTATTACTACAGTTTTAAATATGAGAACTAACGGAATGAGTCTACACAAATTACCTTTATTTGTTTGGGCTATCTTCGTAACAGCTATTCTATTATTATTATCATTACCAGTATTAGCTGGTGCTATTACAATGCTTTTAACTGATAGAAACTTTAATACTAGTTTCTACGATCCAGCAGGTGGAGGTGATCCAATACTTTATCAACATTTATTC",
"TTTATCGTCAAAAGCTGGACAACCTGGTCCTGCAATGGATTTAGCGATTTTAAGTTTACATTTGGCTGGTGCTTCATCAATTTTAGGAGCAGTTAATTTTATCACAACAATTCTTAATATGAGAACTCCTGGAATGACTTTACACAAGATGCCATTATTTGCTTGGTCAATTTTAGTAACAGCTTTTTTATTGTTATTGTCTTTACCAGTTCTAGCAGGAGCGATAACAATGCTATTAACTGATAGGAACTTTGGAACCGCCTTTTTTGAAGCTCAAACTGGAGGAGACCCAGTGTTGTTTCAACATTTATTC",
"GCTAGCCTCAGCAGTCGCCCATAGTGGGGCATCTGTCGATTTAGCAATCTTTTCTCTTCATCTAGCAGGAGCTTCGTCTATTCTGGGGGCTATTAATTTTATTTCCACAGTTATCAATATACGAACCGCGGGAATATTTATAGATCGAGTCCCCCTATTCGTATGATCCGTCTTCATTACTGCGATTCTTCTACTCCTATCTTTGCCAGTCCTAGCCGGAGCTATTACCATACTACTGACTGACCGAAACGTTAACACCTCCTTTTTCGACCCTATGGGCGGAGGAGACCCGATTCTCTATCAACACCTGTTT",
"ACTAAGTACATCATTAATGAGTTTATCATCTTTAGGTATATCATTAGTTATTTATGGTTTAGTATTATTAGGTATATCATCAACTTTAACATCTTTTAACTTCTTTGTTACTTTTGTTTATATGAAATCTTATGGTATGACTTTATCTTCTATGTCAGTTTATGTATGGTCTATTAATATTACAACATCTATGCTATTATTAGTTTTACCAATACTAACTGGAGCTCTTACTATGTTAACTTCAGATATTCATTTTAATACAAGTATCTTTGATTCATTATTTGGTGGTGATCCTGTGTTCTATCAACATTTATTT",
"TCTTGCAGGAAATGTTGCCCATGCCGGAGGTGCAGTTGATGCTGCTATCTTTTCTCTTCATTTAGCCGGAGCTTCTTCTATTTTAGGTGCTGTTAATTTTATTACCACTGTTATTAATATACGAACCCCAAGAATAACTATAGACCGTGTCCCACTCTTTGTATGATCTATTTTTATTACCGTAATTTTATTACTCCTGTCCCTACCTGTTTTAGCAGGGGCTATTACAATGCTTTTGACTGACCGAAATTTAAATACCTCATTCTTTGATCCCAGAGGGGGTGGAGATCCTATTTTATACCAACATTTATTT",
"ATTAAGTTCTATACAAAGTCATTCAGGGGGTGCTGTTGATTTAGCTATTTTTAGTTTACATATAGCAGGAGCTTCATCAATTTTAGGAGCTGTAAATTTTATCTCAACTATTTTAAATATGCGTAATCCTGGCCAAAGTATGTATAGGATGCCTTTATTTGTTTGATCAATTTTTATTACAGCTTTGTTATTATTATTAGCTGTACCTGTTTTAGCGGGAGCTATTACAATGCTTTTAACAGATAGAAATTTTAATACATCTTTTTTTGATCCTGCAGGAGGTGGAGATCCTATATTATATCAACATTTGTTT",
"TTTAGCAGGTATAACAGCTCATTCAGGAGGATCTGTAGATTTAGCAATTTTCAGTCTTCATTTAGCTGGAGCGTCTTCAATATTAGGAGCAATTAATTTTATTTGTACAATATGTAACATGCGTACTGAAAGTTTGCCTTTCCATAAATTGCCTTTATTCGTTTGGGCTGTTCTTATTACTGCAGTATTATTATTATTATCGTTACCTGTTTTAGCTGGAGCTATTACTATGTTATTAACAGATAGAAATTTTAATACTACCTTTTTTGATCCAGCAGGGGGAGGAGATCCTGTACTTTACCAACATTTATTT")


goodTaxa <- setdiff(taxa_names(Final_data_obj), badTaxa)
Final_data_obj <- prune_taxa(goodTaxa, Final_data_obj)

# Remove all control samples ## did not do this to see the controls
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_27266")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_27322")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-1")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-2")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-3")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-4")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-5")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-6")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-7")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-8")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-9")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-10")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-11")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-12")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_00-13")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC1.1")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC1.2")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC1.3")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC2.1")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC2.2")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC3.1")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC3.2")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RRR_PlateC3.3")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RCT_01900")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RCT_01910")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RCT_C1")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RCT_C2")
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "CO1_RCT_C3")

# Remove sample with no reads
Final_data_obj <- subset_samples(Final_data_obj, sample_names(Final_data_obj) != "...")

# removings ASVs that are not present in any sample (if any)
Final_data_obj <- prune_taxa(taxa_sums(Final_data_obj) > 0, Final_data_obj)
Final_data_obj

#26.02. CO1 2025
#otu_table()   OTU Table:          [ 60388 taxa and 232 samples ]:
#  sample_data() Sample Data:        [ 232 samples by 20 sample variables ]:
 # tax_table()   Taxonomy Table:     [ 60388 taxa by 6 taxonomic ranks ]:

#new26.02 with controls :
#otu_table()   OTU Table:          [ 60436 taxa and 253 samples ]:
 # sample_data() Sample Data:        [ 253 samples by 19 sample variables ]:
  #tax_table()   Taxonomy Table:     [ 60436 taxa by 6 taxonomic ranks ]:

saveRDS(Final_data_obj, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/CO1_Control-free_new.rds")
saveRDS(Final_data_obj, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/CO1_with_Controls.rds")

# Look at distribution of library size
df <- as.data.frame(sample_data(Final_data_obj)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(Final_data_obj)
df <- df[order(df$LibrarySize), ]
df$Index <- seq_len(nrow(df))
ggplot(data = df, aes(x = Index, y = LibrarySize, color = SAMPLE_OR_CONTROL)) + geom_point()
sample_data(Final_data_obj)


################################################################
### (6) Clusterinng ASVs into 97% OTUs ####
################################################################
setwd("C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files")

nproc <- 10

COI <- readRDS("CO1_Control-free_new.rds")
sequences <- taxa_names(COI)
sample_names <- sample_names(COI)

# Install decipher

dna <- Biostrings::DNAStringSet(sequences)
set.seed(123) # initialize the random number generator
clusters <- DECIPHER::Clusterize(dna, method = "overlap",
                                 cutoff = 0.03,
                                 penalizeGapLetterMatches = NA,
                                 includeTerminalGaps = TRUE,
                                 processors = nproc)
set.seed(NULL)

# penalize gap-to-letter mismatches once per insertion or deletion,
# which treats runs of gaps (i.e., indels) as equivalent to a single mismatch
# the calculation of distance will use the entire (global) alignment

# Note: function "merge_taxa_vec" was in the package mikemc/speedyseq

# Merge clusters with taxonomy
COI.otu <- merge_taxa_vec(
  COI,
  group = clusters$cluster,
  tax_adjust = 0)

# Checking objects

COI
#otu_table()   OTU Table:          [ 60388 taxa and 232 samples ]:
#sample_data() Sample Data:        [ 232 samples by 20 sample variables ]:
 # tax_table()   Taxonomy Table:     [ 60388 taxa by 6 taxonomic ranks ]:

COI.otu #97%
#otu_table()   OTU Table:          [ 33944 taxa and 232 samples ]:
#  sample_data() Sample Data:        [ 232 samples by 20 sample variables ]:
#  tax_table()   Taxonomy Table:     [ 33944 taxa by 6 taxonomic ranks ]:


#otu_table()   OTU Table:          [ 32045 taxa and 213 samples ]:
#sample_data() Sample Data:        [ 213 samples by 20 sample variables ]:
#  tax_table()   Taxonomy Table:     [ 32045 taxa by 6 taxonomic ranks ]:

# Take a look at files
tax_dada_tax <- tax_table(COI.otu)
write.csv(tax_dada_tax, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/otu_tax_new.csv")

tax_dada_otu <- otu_table(COI.otu, taxa_are_rows = FALSE)
write.csv(tax_dada_otu, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/otu_matrix_new.csv")

################################################################
### (7) OTU curation ####
################################################################
# File preparation
# A - OTU table with samples as columns and OTUs as rows


tax_dada_otu <- t(tax_dada_otu)
tax_dada_otu.df <- phyloseq_to_df(tax_dada_otu, addtax = FALSE, addtot = FALSE, addmaxrank = FALSE,
                                  sorting = "abundance") # with package metagMisc
tax_dada_otu.df  <- data.frame(tax_dada_otu.df, row.names = 1)

# B - Fasta file to prepare the match list
library(dplyr)
library(tibble)

#prepare OTU table and Sequences for matchlist both with the same OTU-IDs because matchlist does it with OTU IDs and this is why lulu doesnt work because the names are not the same
OTU_table_OTU <- tax_dada_otu.df %>%
  rownames_to_column("Seq") %>%                     # Convert row names into a column named "Seq"
  mutate(OTU = paste0("OTU", row_number())) %>%   # Generate OTU IDs based on row number#
  select(OTU, everything())  # Ensure OTU is the first column

OTU_Sequences_match <- OTU_table_OTU %>%
   select(OTU, Seq)                                  # Select only the OTU and Seq columns

# OTU als rownames setzen und die "Seq"-Spalte explizit entfernen
OTU_table_lulu <- OTU_table_OTU %>%
  column_to_rownames("OTU") %>%
  select(-Seq)  # Entfernt die Spalte mit dem Namen "Seq"


fa <- character(2 * nrow(OTU_Sequences_match))
fa[c(TRUE, FALSE)] <- paste0(">", OTU_Sequences_match$OTU)
fa[c(FALSE, TRUE)] <- as.character(OTU_Sequences_match$Seq)

writeLines(fa, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/otus97COI_new.fasta")

# C - Create match list with vsearch
# Place input file "otus97.fasta" in "/Applications/vsearch/bin" folder of vsearch and run:
# vsearch --usearch_global otus97.fasta --db otus97.fasta --self --id .84 --iddef 1 --userout match_list97.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
matchlist_name <- read.table("C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/match_list97_new.txt")
names(matchlist_name)[names(matchlist_name) == "V1"] <- "OTUid"
names(matchlist_name)[names(matchlist_name) == "V2"] <- "hit"
names(matchlist_name)[names(matchlist_name) == "V3"] <- "match"
matchlist_name$OTUid <- as.character(matchlist_name$OTUid)
matchlist_name$hit <- as.character(matchlist_name$hit)


# Run OTU curation
curated_result <- lulu(OTU_table_lulu, matchlist_name)
## look at cureated lulu settings
#curated_result$minimum_match
#curated_result$minimum_relative_cooccurence
#curated_result$discarded_count

# Curated OTU table
curated_table <- curated_result$curated_table
write.csv(curated_table, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated_table97_new.csv")

# keep only the OTUs that passed lulu for taxonomic assignment

# Extract the OTUs that stayed
retained_otus <- rownames(curated_result$curated_table)

# Filter the original OTU + Sequence table to include only retained OTUs for taxonomic assignment
curated_otu_sequences <- OTU_Sequences_match %>%
  filter(OTU %in% retained_otus)
write.csv(curated_otu_sequences, file <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated97_seq_otu_new.csv")


# Curation Results
#  32588 OTUS were obtained. -> less then before

# Prepare fasta file for taxonomic assignment using the curated sequences.

curated_otu_sequences <- read.csv("C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated97_seq_otu_new.csv")
#remove numbers
curated_seqs <- curated_otu_sequences[, -1]
fa <- character(2 * nrow(curated_otu_sequences))
fa[c(TRUE, FALSE)] <- paste0(">", curated_otu_sequences$OTU)
fa[c(FALSE, TRUE)] <- as.character(curated_otu_sequences$Seq)
writeLines(fa, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated_97OTUseqs_COI_new.fasta")


#____________________________________________________________________________
#
## Download latest MIDORI CO1 unique database from http://www.reference-midori.info/download.php#
# Unzip files and place Final_data.nochim.tax.byOTU.fasta in folder

# Run BLASTn
# cd in directory MIDORI_UNIQ_NUC_SP_GB249_CO1_BLAST
# Place curated_97OTUseqs_COI_new.fasta in directory MIDORI_UNIQ_NUC_SP_GB249_CO1_BLAST
# Run following Blast command

# For top hit only
# blastn -db MIDORI2_UNIQ_NUC_GB260_CO1_BLAST -query curated_OTUseqs.fasta -num_alignments 1 -evalue 0.01 -word_size 11 -outfmt "6 qseqid sallseqid evalue bitscore length nident pident" -out curated_OTU-BLAST.out -num_threads 8

# For top ten hits
# blastn -db MIDORI2_UNIQ_NUC_GB261_CO1_BLAST -query curated_97OTUseqs_COI_new.fasta -evalue 0.01 -word_size 11 -culling_limit 100 -outfmt "6 qseqid sallseqid evalue bitscore length nident pident qcovs" -out curated_OTU_BLAST_top10_COI.fasta -num_threads 10 -max_target_seqs 10
#blastn -db MIDORI2_UNIQ_NUC_GB261_CO1_BLAST -query newOTUScuratedforBlastshort.txt -evalue 0.01 -word_size 11 -culling_limit 100 -outfmt "6 qseqid sallseqid evalue bitscore length nident pident qcovs" -out trialcurated_OTU_BLAST_top10_COI.fasta -num_threads 10 -max_target_seqs 10
#blastn -db MIDORI2_UNIQ_NUC_GB261_CO1_BLAST -query newOTUScuratedforBlastshort.txt -evalue 0.01 -word_size 11 -culling_limit 100 -outfmt "6 qseqid sallseqid evalue bitscore length nident pident qcovs" -out trialcurated_OTU_BLAST_top5_COI.fasta -num_threads 10 -max_target_seqs 5
# "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit"
# the lower (5 - and default) number can produce odd results when there are several species with similar high scores.



################################################################
### (9) Curation LCA with galaxy tool ####
################################################################
# Edit output of blastn search so that it works as input in "galaxy-tool-lca"
# Import blast output

#out <- "C:/Users/luisa/Documents/eDNA/R/MIDORI2_UNIQ_NUC_GB260_CO1_BLAST/curated_OTU-BLAST_top10.out"
out <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/BLASTcurated_97OTUseqs_COI_new.out"
out <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/trialcurated_OTU_BLAST_top5_COI.fasta"
#reading the sample data sheet
Blast.out <- data <- read.delim(out, header = FALSE, sep = "\t", row.names = NULL) #added fill=true, um das problem mit fehlenden elementen in einer zeile zu lösen


###original
library(stringr)
Blast.out$V2 <- gsub("root_1;", "", Blast.out$V2)
Blast.out$V2 <- gsub(";", " / ", Blast.out$V2)
Blast.out <- subset(Blast.out, select = -V5)
Blast.out <- subset(Blast.out, select = -V6)
Blast.out[c("#Subject accession", "V2b")] <- str_split_fixed(Blast.out$V2, "###", 2)
Blast.out[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] <- str_split_fixed(Blast.out$V2b, " / ", 7)
Blast.out <- subset(Blast.out, select = -Kingdom)
Blast.out <- subset(Blast.out, select = -Phylum)
Blast.out <- subset(Blast.out, select = -Class)
Blast.out <- subset(Blast.out, select = -Order)
Blast.out <- subset(Blast.out, select = -Family)
Blast.out <- subset(Blast.out, select = -Genus)
Blast.out$Species2 <- gsub("^.*\\_", "", Blast.out$Species)
Blast.out <- subset(Blast.out, select = -V2)

colnames(Blast.out)[colnames(Blast.out) == "V1"] <- "#Query ID"
colnames(Blast.out)[colnames(Blast.out) == "V2b"] <- "#Taxonomy"
colnames(Blast.out)[colnames(Blast.out) == "V3"] <- "#evalue"
colnames(Blast.out)[colnames(Blast.out) == "V4"] <- "#bitscore"
colnames(Blast.out)[colnames(Blast.out) == "V7"] <- "#Identity percentage"
colnames(Blast.out)[colnames(Blast.out) == "V8"] <- "#Coverage"
colnames(Blast.out)[colnames(Blast.out) == "Species"] <- "#Subject"
colnames(Blast.out)[colnames(Blast.out) == "Species2"] <- "#Subject Taxonomy ID"
Blast.out$Source <- "Midori GB261"
colnames(Blast.out)[colnames(Blast.out) == "Source"] <- "#Source"

col_order <- c("#Query ID", "#Subject", "#Subject accession", "#Subject Taxonomy ID", "#Identity percentage", "#Coverage", "#evalue", "#bitscore", "#Source", "#Taxonomy")
Blast.out <- Blast.out[, col_order]
colnames(Blast.out)

# Export as tab delimited text ignoring row names
write.table(Blast.out, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/OTU-LCA_top10_new.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Blast.out, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/trialLCA_top5_new.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# speichern in galaxy-tools ordner nicht im CO1 1. analysis, damit in wsl darauf zugegriffen werden kann
# Summarize Lowest Common Ancestor with galaxy-tool-lca (https://github.com/naturalis/galaxy-tool-lca)

# galaxy-tool must be run with wsl -> navigate to projects (cd), galaxy-tools-lca then copy the lca file in there with cd /mnt/c/Users/l...
#python3 lca.py -i /mnt/c/Users/luisa/Documents/eDNA/R/tools/galaxy-tool-lca/OTU-LCA_top10.txt -o /mnt/c/Users/luisa/Documents/eDNA/R/tools/galaxy-tool-lca/curated_OTU-BLAST_top10.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 97 -tcov 90 -flh unknown

#different settings
#----> this one -->#python3 lca.py -i trialLCA_top10_new.txt -o trialLCA_top10.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 97 -tcov 90 -flh unknown
#---> this one -->#python3 lca.py -i trialLCA_top5_new.txt -o trialLCA_top5.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 98 -tcov 100 -flh unknown
#----> 98 % instead of 97 and best_hit and cov 100 -->#python3 lca.py -i OTU-LCA_top10.txt -o curated_OTU_BLAST_top1098.LCA.txt -b 8 -id 80 -cov 80 -t best_hit -tid 98 -tcov 100 -flh unknown
#----> 98 % instead of 97 and best_hits_range and cov 100 -->#python3 lca.py -i OTU-LCA_top10.txt -o curated_OTU_BLAST_top1098_Range.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 98 -tcov 100 -flh unknown
#----> 99 % instead of 97 and best_hits_range and cov 100 -->#python3 lca.py -i OTU-LCA_top10.txt -o curated_OTU_BLAST_top1098_Range.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 99 -tcov 100 -flh unknown
#----> 99 % instead of 97 and best_hit and cov 100 -->#python3 lca.py -i OTU-LCA_top10.txt -o curated_OTU_BLAST_top1098_Best.LCA.txt -b 8 -id 80 -cov 80 -t best_hit -tid 99 -tcov 100 -flh unknown
#------>98 % instead of 97 and best_hits_range and cov 90 -->#python3 lca.py -i OTU-LCA_top10.txt -o curated_OTU_BLAST_top1098_Range.LCA.txt -b 8 -id 80 -cov 80 -t best_hits_range -tid 98 -tcov 90 -flh unknown
##python3 lca.py -i OTU-LCA_top10_new.txt -o curated_OTU_BLAST_top9890_Range.LCA_new.fasta -b 8 -id 80 -cov 80 -t best_hits_range -tid 97 -tcov 90 -flh unknown
################################################################
# (8) Create curated Phyloseq object ####
################################################################

# Select OTUs in "tax_dada_otu" that remain after curation

# Transpose curated_table
#curated_table = t(curated_table)
#asv_sequences_for_blast_t = t(asv_sequences_for_blast)

# read metadata
Final_data_sam <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Metadatasheet_CO1_NEW1.csv"
Final_data_metadata <- read.csv(Final_data_sam, header = TRUE, sep = ";", row.names = 1)

# Curated taxonomy table

tax_dada <- "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated_BLAST_LCA9790_new.fasta"
tax_dada_tax <- read.csv(tax_dada, header = TRUE, sep = "\t")
tax_dada_tax <- tax_dada_tax[!duplicated(tax_dada_tax$X.Query), ]
tax_dada_tax2 <- tax_dada_tax[, -1]
rownames(tax_dada_tax2) <- tax_dada_tax[, 1]
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.lca.rank"] <- "lca.rank"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.lca.taxon"] <- "lca.taxon"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.kingdom"] <- "Kingdom"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.phylum"] <- "Phylum"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.class"] <- "Class"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.order"] <- "Order"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.family"] <- "Family"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.genus"] <- "Genus"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.species"] <- "Species"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.method"] <- "method"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.identity"] <- "identity"
colnames(tax_dada_tax2)[colnames(tax_dada_tax2) == "X.coverage"] <- "coverage"
write.csv(tax_dada_tax2, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/curated_BLAST_LCA9790_new.csv")



#otu_df_ID[!rownames(otu_df_ID)%in%rownames(tax_dada_tax2),]

#length(unique(tax_dada_tax$X))
#n_occur <- data.frame(table(tax_dada_tax$X))
#n_occur[n_occur$Freq > 1,]
#duplicates <- tax_dada_tax[tax_dada_tax$X %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(duplicates, file="C:/Users/luisa/Documents/eDNA/R/COI 1. Anlysis/duplicates.csv")


library(tibble)

# Convert rownames to a column
curated_OTU_table_merge <- rownames_to_column(curated_table, var = "OTU")
# Convert rownames to a column
tax_dada_tax2_merge <- rownames_to_column(tax_dada_tax2, var = "Query")

# Merge Curated OTU table and taxonomy # changed y= from X to X.Query because there is only a column name called like this in tax.data.tax
curated_table_withtax <- merge(curated_OTU_table_merge, tax_dada_tax2_merge, by.x = "OTU", by.y = "Query", all = FALSE)
write.csv(curated_table_withtax, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/final_curated_table_withtax_new.csv")


Final_data_obj_curated <- phyloseq(otu_table(curated_table, taxa_are_rows = TRUE),
                                   sample_data(Final_data_metadata), tax_table(as.matrix(tax_dada_tax2)))


Final_data_obj_curated

tax_final <- tax_table(Final_data_obj_curated)
otu_final <- otu_table(Final_data_obj_curated)
sam_final <- sample_data(Final_data_obj_curated)


write.csv(tax_final, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/final_tax_new.csv")
write.csv(otu_final, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/final_otu_new.csv")
write.csv(sam_final, file = "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/final_sam_new.csv")


# Save Phyloseq object
saveRDS(Final_data_obj_curated, "C:/Users/luisa/Documents/eDNA/Data_Analysis/R/COI/files/Final_data_obj_curated9097_new.rds")
