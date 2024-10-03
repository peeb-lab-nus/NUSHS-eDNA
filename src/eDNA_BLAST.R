# The following code obtains the taxonomy assignments for the ASVs on the HPC.

if(Sys.info()["user"] == "junyinglim"){
  setwd("")
} else {
  setwd("/Users/jasmine/Documents/Work/eDNA")
}

#################### Split up fasta file for faster processing on the HPC ####################
library(Biostrings)
library(seqinr)

subset_Sequences <- function(input_file, output_file, start, end) {
  # Read the FASTA file
  sequences <- read.fasta(input_file)
  
  # Extract the desired range of sequences
  selected_sequences <- sequences[start:end]
  
  # Write the selected sequences to a new FASTA file
  write.fasta(sequences = selected_sequences, names = names(selected_sequences),
              file.out = output_file)
}
edna_fasta <- read.fasta("/Users/jasmine/Documents/Work/eDNA/master_seq.fasta")
edna <- "/Users/jasmine/Documents/Work/eDNA/master_seq.fasta"
edna_sequences <- subset_Sequences(edna, "/Users/jasmine/Documents/Work/eDNA/master_seq_3001_3430.fasta", 3001, 3430)


#################### Run BLAST on the HPC ####################


#################### Combine blast results from split files + obtain NCBI accession No ####################
library(stringr)
library(dplyr)
library(tidyr)
library(rgbif)

edna_0001_1000 <- read.table("0001_eDNA_GenBankAssignments.txt")
edna_0001_1000$V1 <- gsub("Query_", "", edna_0001_1000$V1)

edna_1001_2000 <- read.table("1001_eDNA_GenBankAssignments.txt")
edna_1001_2000$V1 <- gsub("Query_", "", edna_1001_2000$V1)
edna_1001_2000$V1 <- as.integer(edna_1001_2000$V1) + 1000

edna_2001_3000 <- read.table("2001_eDNA_GenBankAssignments.txt")
edna_2001_3000$V1 <- gsub("Query_", "", edna_2001_3000$V1)
edna_2001_3000$V1 <- as.integer(edna_2001_3000$V1) + 2000

edna_3001_3430 <- read.table("3001_eDNA_GenBankAssignments.txt")
edna_3001_3430$V1 <- gsub("Query_", "", edna_3001_3430$V1)
edna_3001_3430$V1 <- as.integer(edna_3001_3430$V1) + 3000

edna_genbank <- rbind(edna_0001_1000, edna_1001_2000, 
                      edna_2001_3000, edna_3001_3430)
edna_accessionno <- edna_genbank |>
  select(V2) |>
  rename("ID" = "V2")

write.csv(edna_genbank, "eDNA_GenBank_BLASTOutput.csv", row.names = FALSE)
write.csv(edna_accessionno, "eDNA_GenBank_AccessionNo.csv", row.names = FALSE)



#################### Obtain taxonomic information from accession number ####################
#This was run on the HPC
library(dplyr)
library(taxonomizr)

NCBI_AccessionNo <- function(AccessionNoFile){
  ID <- read.csv(AccessionNoFile)
  ID <- as.character(ID$ID)
  sqFile <- "taxdump/accessionTaxa.sql" #load the SQL database
  read.nodes("taxdump/nodes.dmp") #read nodes from downloaded nodes file so that function can fetch taxon info from database
  taxaID <- accessionToTaxa(ID, sqFile) #get taxonomic information ID from accession number 
  TaxaLevels <- getTaxonomy (taxaID, "taxdump/accessionTaxa.sql")
  return(TaxaLevels)
}


#################### Combine taxonomic information with  ####################
library(dplyr)

edna_genbank_blast <- read.csv("eDNA_GenBank_BLASTOutput.csv")
edna_genbank_accessionno_results <- read.csv("eDNA_AccessionNo_Results.csv")

edna_combined <- cbind(edna_genbank_blast, edna_genbank_accessionno_results)
write.csv(edna_combined, "eDNA_GenBank.csv", row.names = FALSE)









