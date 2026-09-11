###################################################################################
##########################  Reduce Labels   #######################################
###################################################################################

######################### Packages
#BiocManager::install("rrvgo")
#BiocManager::install("AnnotationDbi", update = FALSE, force = TRUE)
library(tidyverse)
library(rrvgo)
library(GOSemSim)
library(dplyr)

######################### Importing dataset

data <- read_rds("/Users/myrtomitletton/Documents/BioHackathon2026/gsea_results.rds")

#I recommend that we filter the dataset first, otherwise Generating the new terms takes a lot of time
#p.adj=0.05
#up=T
#library(clusterProfiler)
#data <- data %>% filter(p.adjust < p.adj)
#data <- if (up) filter(data, NES > 0) else filter(data, NES < 0)

######################### Parameters
orgdb = "org.Hs.eg.db"
ont = "all"
method = 'Wang'

######################### Commands
if (ont == "all"){
  ont <- c("CC", "MF", "BP")}

#### Generating the new terms
      reducedTerms <- list()
      reducedTerms <- lapply(ont, function(ONT){
        simMatrix <- calculateSimMatrix(data$ID, orgdb = orgdb, ont = ONT, method = method)
        scores <- setNames(-log10(data$p.adjust), data$ID)
        reducedTerms[[ONT]] <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = orgdb)
      })
      names(reducedTerms) <- ont #columns: go, term, cluster, parent, parentTerm, size, score
      
      # Combine all reducedTerms into one data frame
      all_terms <- do.call(rbind, lapply(ont, function(o) {
        df <- reducedTerms[[o]]
        df$source_ont <- o # Optional: keep track of source if needed
        return(df)
      }))

#### Transferring the labels to data table 
      # Ensure column names match for merging
      # Assuming 'go' is the ID column in reducedTerms
      colnames(all_terms)[colnames(all_terms) == "go"] <- 'ID'
      
      # This keeps all rows in data@result and fills parentTerm where matches exist
      data@result$parentTerm <- NA
      
      data@result <- data@result %>%
        select(-parentTerm) %>% # Remove old column to avoid suffix issues
        left_join(all_terms %>% select(ID, parentTerm), by = 'ID')
      
      
      
      
      
    