## Microbiome dataset (from glmmtree)
library(tidyverse)
load("Data/data1.rda")
# This data is simulated 

# This has several included objects 
# z.tr: 778 OTU proportions for 100 samples in training set
# y.tr: Continous outcome for 100 samples in training set
# z.te: 778 OTU proportions for 200 samples in testing set
# y.te: Continous outcome for 200 samples in training set
# D: Patristic distance matrix among 778 OTUs

# Use the train data 
# Their z is our Y 
mb_data_otu <- z.tr
# Their Y is our X
mb_data_sample <- data.frame(id = 1:200,y = y.tr)


# Convert into phyloseq for easier filtering
library(phyloseq)
otu <- otu_table(mb_data_otu, taxa_are_rows = F)
sam <- sample_data(mb_data_sample)
ps <- phyloseq(otu,sam)
ps

# Helper function for prevalence filtering 
filterTaxaByPrevalence <- function(ps, percentSamplesPresentIn){
  #define threshold
  prevalenceThreshold <- percentSamplesPresentIn * phyloseq::nsamples(ps)
  #apply a function to all samples that determines if a taxa is fully absent from a particular sample
  toKeep <- apply(data.frame(otu_table(ps)), 2, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  #this is placed into toKeep which is a TRUE/FALSE vector of whether or not a that taxa is present in enough samples 
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

# 99% prevalence filter to really filter data down
psf_99 <- filterTaxaByPrevalence(ps, .99)
psf_99
# Arbitrary select only the first 50 of the existing samples. 
psf <- prune_samples(as.numeric(substr(sample_names(psf_99),
                                       3,length(sample_names(psf_99)))) <= 50, psf_99)
psf

# Filter so we have an arbitrarily small dataset to practice building the code with: 
# 50 samples each with 43 observations of OTUs. 

# Convert into the 'long' format where every row is a otu observation for the given 
# sample. Each sample has 43 rows associated with it. 
mb_data_otu <- otu_table(psf)
colnames(mb_data_otu) <- paste0("OTU",1:ntaxa(psf))
mb_data <- cbind(X = y.tr[1:50],
                 sampleID = 1:phyloseq::nsamples(psf),
                 data.frame(mb_data_otu)) %>% 
  pivot_longer(-c(X,sampleID), names_to = "OTU_ID", values_to = "Y") %>% 
  separate(OTU_ID,3, into = c("chr","OTU_ID")) %>% 
  mutate(OTU_ID = as.integer(OTU_ID))



Y_data_small <- mb_data$Y
id_small <- mb_data$sampleID
X_data_small <- y.tr[1:50]

# Test existing geem model using independence cor structure. 
# geem_ind_mb <- geem(formula = Y ~ X, id = sampleID,
#                     data = mb_data) 
# geem_ind_mb

# Corresponding Distance matrix on filtered data 
D_filtered <- D[taxa_names(psf),taxa_names(psf)]




# Proper format for the entire dataset 
mb_data_otu_full <- otu_table(psf_99)
colnames(mb_data_otu_full) <- paste0("OTU",1:ntaxa(psf_99))
mb_data_full <- cbind(X = y.tr[1:50],
                 sampleID = 1:nsamples(psf_99),
                 data.frame(mb_data_otu_full)) %>% 
  pivot_longer(-c(X,sampleID), names_to = "OTU_ID", values_to = "Y") %>% 
  separate(OTU_ID,3, into = c("chr","OTU_ID")) %>% 
  mutate(OTU_ID = as.integer(OTU_ID))

# All samples, but using 99% prevalence. 
Y_data_all_samples <- mb_data_full$Y
id_all_samples <- mb_data_full$sampleID
X_data_all_samples <- y.tr


# Dont want to pre-transform ?

# Do a less strict filter - maybe 30%?

