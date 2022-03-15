
make_long_data <- function(X,ps_y){
  cbind(X = X,
        sampleID = 1:phyloseq::nsamples(ps_y),
        data.frame(otu_table(ps_y))) %>% 
    pivot_longer(-c(X,sampleID), names_to = "OTU_ID", values_to = "Y") %>% 
    separate(OTU_ID,3, into = c("chr","OTU_ID"))  %>% 
    dplyr::select(sampleID, OTU_ID, Y, X)
}
