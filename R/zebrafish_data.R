
library(ape)
library(tidyverse)
library(phy)

zebra_tree <- ape::read.tree(here::here("Data","Zebrafish",
                                        "microbiome_data","alignedsequences.tre"))

tax_final <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/microbiome_data/tax_final.rds")

sequence_table <- readRDS("~/Desktop/phylo_gee_project/Data/Zebrafish/microbiome_data/sequence_table.rds")

tax_table <- tax_table(tax_final)
otu_table <- otu_table(sequence_table, taxa_are_rows = F)
zebra_phyloseq <- phyloseq(tax_table, otu_table, zebra_tree)

pcap_np_metadata <- read_excel("Data/Zebrafish/metadata/pcap_np_metadata.xlsx")
meta1 <- sample_data(pcap_np_metadata)
zebra_phyloseq <- phyloseq(tax_table, otu_table, meta1)
sample_names(meta1)

gaulke_metabolomics_metadata <- read_excel("Data/Zebrafish/metadata/gaulke_metabolomics_metadata.xlsx")
meta2 <- sample_data(gaulke_metabolomics_metadata)
zebra_phyloseq <- phyloseq(tax_table, otu_table, meta2)
sample_names(meta2)

pcap_metabolomics_metadata <- read_excel("Data/Zebrafish/metadata/pcap_metabolomics_metadata.xlsx")


library(ape)
# tests: 
tips <- zebra_tree$tip.label[1:100]
pruned.tree<-drop.tip(zebra_tree,zebra_tree$tip.label[-match(tips, zebra_tree$tip.label)])
write.tree(pruned.tree)

# Try using the distTips function 
# from adephylo 
#install.packages("adephylo")
library(adephylo)
zebra_dist <- distTips(pruned.tree, method = "patristic")
zebra_dist

zebra_dist_full <- distTips(zebra_tree, method = "patristic")


length(zebra_tree$tip.label)
full_ps

tax_table(full_ps) %>% View()

## Paired data
