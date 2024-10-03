
library(dplyr)

if(Sys.info()["user"] == "junyinglim"){
  setwd("")
} else {
  setwd("/Users/jasmine/Documents/Work/eDNA")
}



edna_taxa <- read.csv("eDNA_GenBank.csv")

#################### GenBank Sequences Results ####################
genbankresults <- edna_taxa |> #only selected for similarity score and e-value, feel free to change
  select(V1, V3, V11, superkingdom, phylum, class, order, family, genus, species) |>
  dplyr::rename(ASV = V1,
                Similarity = V3,
                E_value = V11,
                Superkingdom = superkingdom,
                Phylum = phylum,
                Class = class,
                Order = order,
                Family = family,
                Genus = genus,
                Species = species)


#################### Removal of Non-Target Taxa ####################
cleantaxa <- genbankresults |>
  filter(Superkingdom == "Eukaryota") |> #can filter by phylum accordingly
  select(-Superkingdom)

#################### Consensus Taxonomy Assignments ####################
top3hits <- cleantaxa |>
  filter(E_value < 0.001) |>
  select(-E_value) |>
  group_by(ASV) |> 
  slice_head(n = 3) |>
  mutate(No_of_Hits = n()) 

consensus_assignments <- top3hits |>
  mutate(Species = ifelse(No_of_Hits == 3, Species, NA),
         Species = ifelse((No_of_Hits == 3 & all(Similarity >= 99) & n_distinct(Species) == 1), Species, NA),
         Genus = ifelse((all(Similarity >= 95) & n_distinct(Genus) == 1), Genus, NA),
         Family = ifelse((all(Similarity >= 95) & n_distinct(Family) == 1), Family, NA),
         Order = ifelse((all(Similarity >= 90) & n_distinct(Order) == 1), Order, NA),
         Class = ifelse((all(Similarity >= 90) & n_distinct(Class) == 1), Class, NA),
         Phylum = ifelse((all(Similarity >= 90) & n_distinct(Phylum) == 1), Phylum, NA)) |>
  select(-No_of_Hits, -Similarity) |>
  distinct() 

write.csv(consensus_assignments, "edna_consensus.csv", row.names = FALSE)

















