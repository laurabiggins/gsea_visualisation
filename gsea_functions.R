#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#install.packages("tidyverse")
#install.packages('ggraph')
#install.packages('tidygraph')

library(readr)
library(clusterProfiler)

gsea_results <- read_rds("data/gsea_results.rds")

gsea_results@result$p.adjust

gsea_results@result <- gsea_results@result[gsea_results@result$p.adjust < 0.05, ]

simplified_views_file <- "data/gsea_simplify_cache.rds"

simplify_cutoffs <- c(0.5, 0.6, 0.7, 0.8, 0.9)
simplified_views <- simplify_cutoffs |>
  set_names(as.character(simplify_cutoffs)) |>
  map(function(cutoff) {
    message("Computing simplify() for cutoff = ", cutoff)
    clusterProfiler::simplify(
      gsea_results,
      cutoff = cutoff,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    )
  })

simplified_views[["1"]] <- gsea_results

write_rds(simplified_views, simplified_views_file)
