## Sequencing data processing using DADA2

## PACKAGES AND DIRECTORIES =======================
#BiocManager::install("dada2", version = "3.19")
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)

if(Sys.info()["user"] == "junyinglim"){
  main.dir <- "~/Dropbox/Projects/NUSHS-eDNA/"
} else {
  main.dir <- "" # specify your main path here
}

data.dir <- file.path(main.dir, "data_raw")
output.dir <- file.path(main.dir, "data_processed")
fig.dir <- file.path(main.dir, "figures")
if(!dir.exists(output.dir)){dir.create(output.dir)}
if(!dir.exists(fig.dir)){dir.create(fig.dir)}

## DEFINE SAMPLES =======================
# Define sample ids (TO DO)
sample_ids <- c("A1_1", "A1_2")

# Define input parameter 
if(Sys.info()[["user"]] == "junyinglim"){
  multithread <- TRUE  
} else {
  multithread <- FALSE
}

run_ids <- file.path(data.dir, sample_ids[1]) %>%
  list.files(pattern = "\\.fq\\.gz") %>%
  gsub(., pattern = "_1.fq.gz|_2.fq.gz", replacement = "") %>%
  gsub(., pattern = paste0(sample_ids[1], "_"), replacement = "") %>%
  #gsub(., pattern = "_L1|_L2", replacement = "") %>%
  unique()

run_id <- run_ids[2] # RERUN THE CODE FOR ALL RUN IDS

# Get file names
fname <- list.files(data.dir, recursive = TRUE, pattern = "\\.fq\\.gz", full.names = TRUE)

# Get file names for this unique run only
run_id_fnames <- fname[grepl(pattern = run_id, x= fname)]

# Get files names for the forward reads for this unique run only
run_id_for_fnames <- run_id_fnames[grepl(run_id_fnames, pattern = "_1.fq.gz")]

# Get files names for the reverse reads for this unique run only
run_id_rev_fnames <- run_id_fnames[grepl(run_id_fnames, pattern = "_2.fq.gz")]

## INPUT PRIMER SEQUENCE ===============
# Checking that the right primers are present + correct orientation
FWD <- 'GGWACWGGWTGAACWGTWTAYCCYCC' # forward primer
REV <- 'TANACYTCNGGRTGNCCRAARAAYCA' # reverse primer

allOrients <- function(primer) {  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna,
               Complement = Biostrings::complement(dna),
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

## TRIM AMBIGUOUS BASES ===============
# Define subdirectory
trim.dir <- file.path(output.dir, "1_Trimmed")

# Define file names
fnFs.filtN <- file.path(trim.dir, paste(sample_ids, run_id, "1_trim.fq.gz", sep = "_")) 
fnRs.filtN <- file.path(trim.dir, paste(sample_ids, run_id, "2_trim.fq.gz", sep = "_"))

# Filter and trim
filterAndTrim(fwd = run_id_for_fnames,
              filt = fnFs.filtN,
              rev = run_id_rev_fnames,
              filt.rev = fnRs.filtN,
              maxN = 0,
              multithread = multithread)

#Count the number of times the primers appear in the forward and reverse read while accounting for all possible primer orientations
primerHits <- function(primer, fn) { # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer,
                         sread(readFastq(fn)),
                         fixed = FALSE)
  return(sum(nhits > 0))
}
before_RemovingPrimers <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), #Only processes the first sample since library preparation is assumed to be the same
                                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
                                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
print(before_RemovingPrimers)
# The reverse complement is low because the insert size is larger than 250 bp

## TRIM PRIMERS ====================
# Define cutadapt path
if(Sys.info()[["user"]] == "junyinglim"){
  cutadapt <- "/Users/junyinglim/my_virtual_env/bin/cutadapt"
  system2(cutadapt, args = "--version") #Run shell commands from R to check that cutadapt is working; version 4.6 as of Feb 2024; on HPC, put cutadapt down as a string (i.e., "cutadapt")
} else {
  ## TO DO
}

# Define output directory for cuts
cut.dir <- file.path(output.dir, "2_Cutadapt")
if(!dir.exists(cut.dir)) dir.create(cut.dir) # if folder doesn't exist, make one

# Define output file names
fnFs.cut <- file.path(cut.dir, paste(sample_ids, run_id, "1_cutadapt.fq.gz", sep = "_"))
fnRs.cut <- file.path(cut.dir, paste(sample_ids, run_id, "2_cutadapt.fq.gz", sep = "_"))

# Define reverse complement for primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Defining search parameters for forward reads
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Defining search parameters for reverse reads
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  system2(cutadapt, # HPC, change to string
          args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                   "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                   fnFs.filtN[i], fnRs.filtN[i], # input files
                   "--minimum-length", 1)) ###### minimum length of 1
}

# Checlking to see if primers have been removed properly (should be 0)
after_RemovingPrimers <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
                               FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
                               REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
                               REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
print(after_RemovingPrimers)
print(before_RemovingPrimers)

## FILTER AND TRIM BASED ON READ QUALITY ================
# Inspect read quality profiles
FWD_QualityScores <- plotQualityProfile(fnFs.cut)
ggsave(paste(run_id, "FWD_QualityScores.pdf", sep = "_"),
       FWD_QualityScores,
       path = fig.dir,
       width = 5,
       height = 7)

REV_QualityScores <- plotQualityProfile(fnRs.cut)
ggsave(paste(run_id, "REV_QualityScores.pdf", sep = "_"),
       REV_QualityScores,
       path = fig.dir,
       width = 5,
       height = 7)

# Define output directory for this step
filter.dir <- file.path(output.dir, "3_Filtered")
if(!dir.exists(filter.dir)){dir.create(filter.dir)}

# Define output files for this step
filtFs <- file.path(filter.dir,
                    paste(sample_ids,
                          run_id, 
                          "1_filtered.fastq.gz",
                          sep = "_"))
filtRs <- file.path(filter.dir,
                    paste(sample_ids,
                          run_id, 
                          "2_filtered.fastq.gz",
                          sep = "_"))

out <- filterAndTrim(fnFs.cut,
                     filtFs,
                     fnRs.cut,
                     filtRs,
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2, 
                     minLen = 50, # min length default == 20 bp
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = multithread)
# on windows, set multithread = FALSE
print(out) # no of reads that have been filtered and trimmed

# Learn the error rates
errF <- learnErrors(filtFs, multithread = multithread)
errR <- learnErrors(filtRs, multithread = multithread)

For_Error_Rates_Plot <- plotErrors(errF, nominalQ = TRUE)
ggsave(paste(run_id, "For_Error_Rates.pdf", sep = "_"),
       For_Error_Rates_Plot,
       path = fig.dir,
       width = 5,
       height = 7)

Rev_Error_Rates_Plot <- plotErrors(errF, nominalQ = TRUE)
ggsave(paste(run_id, "Rev_Error_Rates.pdf", sep = "_"),
       Rev_Error_Rates_Plot,
       path = fig.dir,
       width = 5,
       height = 7)

## SAMPLE INTERFENCE ====================
# Sample inference
dadaFs <- dada(filtFs, err = errF, multithread = multithread)
dadaRs <- dada(filtRs, err = errR, multithread = multithread)

# Merged paired reads
mergers <- mergePairs(dadaF = dadaFs,
                      derepF = filtFs,
                      dadaR = dadaRs,
                      derepR = filtRs,
                      minOverlap = 50,
                      verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # rows are samples, columns are unique sequences
table(nchar(getSequences(seqtab))) # most are 313 bp long

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method="consensus",
                                    multithread=multithread,
                                    verbose=TRUE)

# Number of unique sequences should be same or lower
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# Proportion of unique reads that are not chimeras
sum(seqtab.nochim)/sum(seqtab) 
table(nchar(getSequences(seqtab.nochim)))

# Export results
asv.dir <- file.path(output.dir, "4_Final")
if(!dir.exists(asv.dir)) {dir.create(asv.dir)}
unique_seq_in_sample <- colnames(seqtab.nochim)
write.csv(data.frame("sample.id" = run_id,
                     "sample.seq.id" = 1:length(unique_seq_in_sample),
                     "seq" = unique_seq_in_sample),
          file = file.path(asv.dir, paste0(run_id, "_sample_seq.csv")),
          row.names = FALSE)
colnames(seqtab.nochim) <- 1:length(unique_seq_in_sample)
rownames(seqtab.nochim) <- sample_ids
write.csv(seqtab.nochim,
          file = file.path(asv.dir, paste0(run_id,"_sample_asv_table.csv")),
          row.names = TRUE)

# Tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_ids
head(track)