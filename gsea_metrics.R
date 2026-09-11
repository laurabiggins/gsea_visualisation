library(clusterProfiler)
library(tidyverse)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("clusterProfiler")

#install.packages("enrichit")
# BiocManager::install("org.Hs.eg.db")
#install.packages("enrichplot")
#install.packages("ggraph")
#install.packages("tidygraph")
#BiocManager::install("pathview")

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggraph)
library(tidygraph)
library(pathview)


setwd("~/OneDrive/Bureau/Perso/Fac/Cambridge/Academic PDN/Year 2026-2027/gsea_project_biohackathon")
gsea_results <- read_rds("data/gsea_results.rds")


col_gradient <- c(
  "#2166AC",
  "#FFFFFF",
  "#B2182B"
)


make_gsea_graph <- function(
    data,
    min_in_cluster = 2,
    gene_size = FALSE,
    metrics = "NES",
    col_gradient = NULL,
    show = "ALL"
) {
  
  data |>
    filter(p.adjust < 0.05) -> data
  
  
  # Check show argument
  if (!(show %in% c("ALL", "DOWN", "UP"))) {
    stop("show should be either 'ALL', 'DOWN' or 'UP'.")
  }
  
  
  # Select up/down pathways based on NES
  if (show == "UP") {
    
    data |>
      filter(NES > 0) -> data
    
  } else if (show == "DOWN") {
    
    data |>
      filter(NES < 0) -> data
  }
  
  
  # Check metric
  if (!is.null(metrics)) {
    
    if (length(metrics) != 1) {
      stop("Please select only one metric.")
    }
    
    if (!(metrics %in% colnames(data@result))) {
      stop(
        paste0(
          "'", metrics,
          "' is not present in the GSEA results."
        )
      )
    }
  }
  
  
  # Pairwise term similarity
  data |>
    pairwise_termsim() -> data_termsim
  
  data_termsim@termsim |>
    as_tibble(rownames = "from") |>
    pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "weight"
    ) |>
    filter(from != to) |>
    filter(weight > 0.3) -> data_pairs
  
  
  # Graph
  tbl_graph(
    edges = data_pairs,
    directed = FALSE
  ) -> data_graph
  
  data_graph |>
    activate(nodes) |>
    mutate(component = group_components()) |>
    group_by(component) |>
    mutate(component_size = n()) |>
    ungroup() |>
    filter(component_size >= min_in_cluster) -> data_graph
  
  
  # Add gene set size
  if (gene_size) {
    
    data@result |>
      mutate(
        size = lengths(strsplit(core_enrichment, "/"))
      ) |>
      select(Description, size) -> node_size
    
    data_graph <- data_graph |>
      activate(nodes) |>
      left_join(
        node_size,
        by = c("name" = "Description")
      )
  }
  
  
  # Add selected metric
  if (!is.null(metrics)) {
    
    data@result |>
      select(
        Description,
        all_of(metrics)
      ) -> node_metric
    
    data_graph <- data_graph |>
      activate(nodes) |>
      left_join(
        node_metric,
        by = c("name" = "Description")
      )
  }
  
  
  set.seed(17434)
  
  plot <- ggraph(
    data_graph,
    layout = "fr"
  ) +
    geom_edge_link(alpha = 0.2)
  
  
  # Nodes
  if (gene_size && !is.null(metrics)) {
    
    plot <- plot +
      geom_node_point(
        aes(
          fill = .data[[metrics]],
          size = size
        ),
        pch = 21
      )
    
  } else if (gene_size) {
    
    plot <- plot +
      geom_node_point(
        aes(
          fill = as.factor(component),
          size = size
        ),
        pch = 21
      )
    
  } else if (!is.null(metrics)) {
    
    plot <- plot +
      geom_node_point(
        aes(
          fill = .data[[metrics]]
        ),
        pch = 21,
        size = 4
      )
    
  } else {
    
    plot <- plot +
      geom_node_point(
        aes(
          fill = as.factor(component)
        ),
        pch = 21,
        size = 4)
  }
  # Metric colour gradient
  if (!is.null(metrics) && !is.null(col_gradient)) {
    plot <- plot +
      scale_fill_gradientn(colours = col_gradient)
  }
  plot <- plot +
    theme_graph() +
    theme(
      legend.position = "right"
    )
  # Gene set size scale
  if (gene_size) {
    plot <- plot +
      scale_size_continuous(
        range = c(1, 10)
      ) +
      labs(size = "Gene set size")
  }
  # Metric legend title
  if (!is.null(metrics)) {
    plot <- plot +
      labs(fill = metrics)
  } else {
    plot <- plot +
      guides(fill = "none")
  }
  return(list(plot,data_graph))
}

gsea_results |>
  as_tibble() |>
  arrange(p.adjust) |> 
  filter(p.adjust<0.05)

gsea_graph <- make_gsea_graph(gsea_results, gene_size = TRUE, show = "ALL", metrics = "NES", col_gradient = col_gradient)
gsea_graph[[1]]

