## RECONCILING SEQUENCE TABLES

## PACKAGES AND DIRECTORIES =======================
#BiocManager::install("dada2", version = "3.19")
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr)
library(tidyverse)

if(Sys.info()["user"] == "junyinglim"){
  main.dir <- "~/Dropbox/Projects/NUSHS-eDNA/"
} else {
  main.dir <- "" # specify your main path here
}

output.dir <- file.path(main.dir, "data_processed")
seq.dir <- file.path(output.dir, "4_Final")
clean.dir <- file.path(output.dir, "5_Clean")
if(!dir.exists(clean.dir)){dir.create(clean.dir)}

sample_id <- c("FKDL240103905-1A_HFKHCDRX5_L2",
               "FKDL240103905-1A_HFM7VDRX5_L2",
               "FKDL240103905-1A_HFV3YDRX5_L1")

## IMPORT DATA =======================
# Import sequence tables
df_list <- list.files(x, pattern = "sample_seq", full.names = TRUE) %>%
  lapply(FUN = function(x) read.csv(x))

# Generate master seq id from the union of all individual datasets
uniq_seq_vect <- lapply(df_list, FUN = function(x) { x$seq}) %>%
  unlist() %>%
  unique()
seq_master_db <- data.frame(master.seq.id = 1:length(uniq_seq_vect),
                            seq = uniq_seq_vect)

# Write a fasta file
fasta_write_list <- list()
counter <- 1
for(i in 1:nrow(seq_master_db)){
  fasta_write_list[[counter]] <- paste0("> master.seq.id.", seq_master_db$master.seq.id[i])
  fasta_write_list[[counter + 1]] <- seq_master_db$seq[i]
  counter <- counter + 2
}
writeLines(text = unlist(fasta_write_list),
           con = file.path(clean.dir, "master.seq.fasta"),
           sep = "\n")

# # Example of how match works
# match(c("D", "E"), c("A", "B", "C", "D", "E"))

# THis code tells you what is the master seq id
df_list[[1]]$master.seq.id <- seq_master_db$master.seq.id[match(df_list[[1]]$seq, seq_master_db$seq)]
df_list[[2]]$master.seq.id <- seq_master_db$master.seq.id[match(df_list[[2]]$seq, seq_master_db$seq)]
df_list[[3]]$master.seq.id <- seq_master_db$master.seq.id[match(df_list[[3]]$seq, seq_master_db$seq)]

# Create a new reference table
seq_df_combined <- do.call("rbind", df_list) %>%
  .[,c("sample.seq.id", "sample.id", "master.seq.id")] %>%
  unique()

# Import seq tables 
table_list <- list.files(x, pattern = "asv", full.names = TRUE) %>%
  lapply(FUN = function(x) read.csv(x)) %>%
  lapply(FUN = function(x) reshape2::melt(data = x))
  
for(i in 1:length(sample_id)){
  table_list[[i]] <- table_list[[i]] %>%
    mutate(sample.id = sample_id[i]) %>%
    dplyr::rename(sample.seq.id = variable) %>%
    mutate(sample.seq.id = as.numeric(gsub(sample.seq.id, pattern = "X", replacement = ""))) %>%
    dplyr::rename(sample = X)
}

seq_table_combined <- do.call("rbind",table_list) %>%
  left_join(seq_df_combined, by = c("sample.seq.id","sample.id")) %>%
  mutate(sample_main = stringr::str_split_fixed(string = sample, n = 2, pattern = "_")[,1])

write.csv(seq_table_combined,
          file = file.path(clean.dir, "seq_table_combined.csv"),
          row.names = FALSE)


test <- plyr::ddply(.data = seq_table_combined,
      .variables = c("sample_main", "master.seq.id"),
      .fun = summarise,
      present = ifelse(sum(value) > 0, 1, 0),
      .progress = "text")
head(test)
