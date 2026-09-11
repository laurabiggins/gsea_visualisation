###################################################################################
##########################  Reduce Labels   #######################################
###################################################################################

######################### Packages
library(rrvgo)
library(GOSemSim)

######################### Function

gsea_results2 <- read_rds("/Users/myrtomitletton/Documents/BioHackathon2026/gsea_results.rds")

ReduceLabels <- function(data, up = TRUE, GO_ID='ID', 
                            orgdb = "org.Hs.eg.db", ont = "all", p.adj=0.05, method = 'Wang', treemapPlot=F) {
  
  #method = c("Resnik", "Lin", "Rel", "Jiang", "Wang")
  #ont = c('CC', 'BP', 'MF', 'all')
  
  #data <- data %>% filter(p.adjust < p.adj)
  #data <- if (up) filter(data, NES > 0) else filter(data, NES < 0)
  
  
  if (ont == "all"){
    ont <- c("CC", "MF", "BP")}
  
  reducedTerms <- list()
  reducedTerms <- lapply(ont, function(ONT){
    simMatrix <- calculateSimMatrix(data$GO_ID, orgdb = orgdb, ont = ONT, method = method)
    scores <- setNames(-log10(data$p.adjust), data$GO_ID)
    reducedTerms[[ONT]] <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = orgdb)
  })
  names(reducedTerms) <- ont #columns: go, term, cluster, parent, parentTerm, size, score
  # Combine all reducedTerms into one data frame
  all_terms <- do.call(rbind, lapply(ont, function(o) {
    df <- reducedTerms[[o]]
    df$source_ont <- o # Optional: keep track of source if needed
    return(df)
  }))
  
  #Transferring the labels to data
  # Ensure column names match for merging
  # Assuming 'go' is the ID column in reducedTerms
  colnames(all_terms)[colnames(all_terms) == "go"] <- GO_ID
  
  # Perform a left join
  # This keeps all rows in data@result and fills parentTerm where matches exist
  # If you can use tidyverse
  data@result$parentTerm <- NA
  
  data@result <- data@result %>%
    select(-parentTerm) %>% # Remove old column to avoid suffix issues
    left_join(all_terms %>% select(GO_ID, parentTerm), by = GO_ID)
  
  # If you need base R only:
  # merged <- merge(data@result, all_terms[, c(GO_ID, "parentTerm")], by = GO_ID, all.x = TRUE)
  # data@result <- merged
  # Note: merge() might reorder rows, so match by ID afterwards if row order matters strictly.
  
}