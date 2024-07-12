# Load the igraph package
library(igraph)
#plot_tree
plot_dp_comp_exposure <- function(hdpsample){
  # hdpsample: HDPSample object
  # Install igraph package if not already installed
  if (!requireNamespace("igraph", quietly = TRUE)) {
    install.packages("igraph")
  }



  # Your parent vector
  parent_vector <- c(0, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3)

  # Create edge list
  edges <- c()
  for (i in 2:length(parent_vector)) {
    edges <- c(edges, parent_vector[i], i)
  }

  # Create the graph from the edge list
  g <- graph(edges, directed = TRUE)

  # Plot the graph
  plot(g, vertex.label = V(g)$name, main = "Tree Structure from Parent Vector")
}
