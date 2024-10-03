
seq_table_combined <- read.csv("~/Dropbox/Projects/NUSHS-eDNA/data_processed/5_Clean/seq_table_combined.csv")
test <- plyr::ddply(.data = seq_table_combined,
                    .variables = c("sample_main", "master.seq.id"),
                    .fun = summarise,
                    present = ifelse(sum(value) > 0, 1, 0),
                    .progress = "text")

contamination_ids <- subset(test, sample_main %in% c("EXTC", "PCR") & present == 1)$master.seq.id


test2 <- subset(test, (!sample_main %in% c("EXTC", "PCR")) & (!master.seq.id %in% contamination_ids))

nrow(test)
nrow(test2)

unique(test2$sample_main)
contamination_ids %in% unique(test2$master.seq.id)

taxon_assignments <- read.csv("~/Dropbox/Projects/NUSHS-eDNA/data_processed/6_TaxonAssignment/eDNA_consensusassignments.csv")
head(taxon_assignments)

# Merge taxon assignments with my cleaned-up file
head(test2)
test3 <- dplyr::left_join(x = test2,
                 y = taxon_assignments,
                 by = c("master.seq.id" = "ASV"))
test3$Phylum[is.na(test3$Phylum)] <- "Unassigned"
test3$Order[is.na(test3$Order)] <- "Unassigned"
test3$Family[is.na(test3$Family)] <- "Unassigned"

head(test3)

phylum_counts<- ddply(.data = subset(test3, present == 1),
                      .variables = .(sample_main, Phylum),
                      .fun = summarise,
                      n_seq = length(master.seq.id))

library(ggplot2)
ggplot(data = phylum_counts) +
  geom_bar(aes(x = sample_main, y = n_seq, fill = Phylum), stat = "identity")

subset(test3, grepl(sample_main, pattern = "D"))
test3

arthropod_counts<- ddply(.data = subset(test3, present == 1 & Phylum == "Arthropoda"),
                      .variables = .(sample_main),
                      .fun = summarise,
                      n_seq = length(master.seq.id))
median(arthropod_counts$n_seq)

arthropod_order_counts <- ddply(subset(test3, Phylum == "Arthropoda" & present == 1),
      .variables = .(sample_main, Order),
      .fun = summarise,
      n_seq = length(master.seq.id))

ggplot(data = arthropod_order_counts) +
  geom_bar(aes(x = sample_main, y = n_seq, fill = Order), stat = "identity")

## Convert this into a 
pollinator_dataset <- subset(test3, Phylum == "Arthropoda" &
         Order %in% c("Coleoptera",
                      "Hymenoptera",
                      "Diptera",
                      "Lepidoptera",
                      "Entomobryomorpha"))
head(pollinator_dataset)
pollinator_matrix<- acast(sample_main ~ master.seq.id, fill = 0, value.var = "present",
                          data = pollinator_dataset)

#install.packages("vegan")
library(vegan)
pollinator_dissim <- vegan::vegdist(x = pollinator_matrix,
                                    method = "jaccard")
pollinator_dissim_df <- reshape2::melt(as.matrix(pollinator_dissim))

ggplot(data = pollinator_dissim_df) +
  geom_tile(aes(y = Var1, x = Var2, fill = value))
