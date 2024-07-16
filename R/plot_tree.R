# Load the igraph package
library(igraph)
#plot_tree
plot_tree <- function(hdpsample){
  # hdpsample: HDPSample object
  # Install igraph package if not already installed

 #rpart.plot
  # Your parent vector
  parent_vector <- hdpsample@ppindex

  plot_tree_from_parent_vector(parent_vector)

}

plot_tree_from_parent_vector <- function(parent_vector){
  # Create edge list
  edges <- c()
  for (i in 2:length(parent_vector)) {
    edges <- c(edges, parent_vector[i], i)
  }

  # Create the graph from the edge list
  g <- graph(edges, directed = TRUE)


  # Plot the tree
  plot(g, layout = layout_as_tree(g, root = 1),
       vertex.label = 1:length(ppindex),
       vertex.size = 30,
       vertex.color = "skyblue",
       edge.arrow.size = 0.5,
       main = "Tree Plot from Parent Index Vector")
}

