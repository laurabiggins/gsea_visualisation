###################################################################################
##########################  Reduce Labels   #######################################
###################################################################################

######################### Packages
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#install.packages("tidyverse")
library(tidyverse)
#BiocManager::install("rrvgo")
library(rrvgo)
#BiocManager::install("GOSemSim")
library(GOSemSim)
#install.packages("dplyr")ß
library(dplyr)

######################### Function
ReduceLabels <- function(data, orgdb = "org.Hs.eg.db", ont = "all", method = 'Wang', treemapPlot=F, folder.out=''){ #up = TRUE, GO_ID='ID', p.adj=0.05
  
  #method = c("Resnik", "Lin", "Rel", "Jiang", "Wang")
  #ont = c('CC', 'BP', 'MF', 'all')
  
  if (ont == "all"){
    ont <- c("CC", "MF", "BP")
    }
  
  #Generating the new terms
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
  
  #Creating optional plot if specified by the user with `treemapPlot=F` and `folder.out=''`
  if (treemapPlot){
    if (!dir.exists(folder.out)) {
      dir.create(folder.out, recursive = TRUE)
    }
    pdf(file.path(folder.out, "ReducedTermsTreemapPlot.pdf"), width = 10, height = 8)
    treemapPlot(all_terms)
    dev.off()
  }
  
  #Transferring the labels to data
  # Ensure column names match for merging
  # Assuming 'go' is the ID column in reducedTerms
  colnames(all_terms)[colnames(all_terms) == "go"] <- 'ID'
  
  # This keeps all rows in data@result and fills parentTerm where matches exist
  data@result$parentTerm <- NA
  
  data@result <- data@result %>%
    select(-parentTerm) %>% # Remove old column to avoid suffix issues
    left_join(all_terms %>% select(ID, parentTerm), by = 'ID')
  
    return(list(data=data, reduced_terms=all_terms))
}