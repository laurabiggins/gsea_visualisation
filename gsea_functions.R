library(clusterProfiler)
library(tidyverse)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("clusterProfiler")

# install.packages("enrichit")
#BiocManager::install("org.Hs.eg.db")
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


gsea_results <- read_rds("data/gsea_results.rds")


make_gsea_graph <- function(data, up=TRUE, min_in_cluster = 2) {
  
  data |>
    filter(p.adjust<0.05) -> data
  
  if (up) {
    data |>
      filter(NES>0) -> data
  } else {
    data |>
      filter(NES<0) -> data
  }
  
  data |>
    pairwise_termsim() -> data_termsim
  
  data_termsim@termsim |>
    as_tibble(rownames="from") |>
    pivot_longer(
      cols=-from,
      names_to="to",
      values_to="weight"
    ) |>
    filter(from != to) |>
    filter(weight > 0.3) -> data_pairs
  
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
  
  set.seed(17434)
  ggraph(
    data_graph,
    layout = "fr"
  ) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes(fill=as.factor(component)),pch=21, size=4,show.legend = TRUE) +
    geom_node_text(aes(label=name), repel=TRUE, size=3, colour="black") +
    theme_graph() +
    theme(legend.position = "none") -> plot
  
  return(list(plot,data_graph))
}

gsea_results |>
  as_tibble() |>
  arrange(p.adjust) |> 
  filter(p.adjust<0.05)

gsea_up_graph <- make_gsea_graph(gsea_results, up=TRUE)
gsea_up_graph[[1]]

gsea_down_graph <- make_gsea_graph(gsea_results, up=FALSE)
gsea_down_graph[[2]]
