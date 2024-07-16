#create tree
library(ape)
library(dplyr)

# Function to convert edge matrix to parent array
edge_matrix_to_parent_array <- function(edge_matrix) {
  if (is.matrix(edge_matrix)){
    edge_matrix <- as.data.frame(edge_matrix)
    edge_matrix <- edge_matrix %>%
      rename(from = V1, to = V2, name = V3)
  }

  all_nodes <- unique(c(edge_matrix$from, edge_matrix$to))
  num_nodes <- length(all_nodes)
  # Initialize parent array with NA values
  parent_array <- rep(NA, num_nodes)
  names_child  <- rep(NA, num_nodes)

  #Populate the parent array with the parent of each node
  for (i in 1:nrow(edge_matrix)) {
    parent <- as.integer(edge_matrix[i, 1])
    child <- as.integer(edge_matrix[i, 2])
    name_i <- edge_matrix[i, 3]
    parent_array[child] <- parent
    names_child[child] <- name_i
  }

  if(sum(is.na(parent_array)) != 1){
    browser()
    stop("There are more than one NA in parent_array")
  }

  parent_array[is.na(parent_array)] <- 0
  names_child[is.na(parent_array)] <- "root"
  names(parent_array) <- names_child
  #print(parent_array)

  parent_array_conv <- convert_tree_representation(parent_array)
  return(parent_array_conv)
}

# Example usage
#parent_array <- c(7, 7, 8, 8, 5, -1, 5, 6, 6)
#parent_array <- parent_array + 1
#new_parent_array_ordered <- convert_tree_representation(parent_array)
#new_parent_array_ordered
convert_tree_representation <- function(parent_array) {
  #artifact from python
  parent_array <- parent_array -1

  # Step 1: Create a list of children for each node to easily traverse the tree
  children <- vector("list", length(parent_array))
  for (i in 1:length(parent_array)) {
    parent <- parent_array[i]
    if (parent != -1) {  # Skip the root node
      # Adjust for R's 1-based indexing
      if (is.null(children[[parent + 1]])) {
        children[[parent + 1]] <- c(i - 1)  # Store original 0-based index
      } else {
        children[[parent + 1]] <- c(children[[parent + 1]], i - 1)
      }
    }
  }

  # Step 2: BFS to create a mapping from old indices to new indices
  old_to_new <- list()
  new_to_old <- list()
  queue <- list(which(parent_array == -1) - 1)  # Start with the root node, adjust for 0-based indexing
  new_index <- 0

  while (length(queue) > 0) {
    old_index <- queue[[1]]
    queue <- queue[-1]  # Remove the first element
    old_to_new[[as.character(old_index)]] <- new_index
    new_to_old[[as.character(new_index)]] <- old_index
    new_index <- new_index + 1
    if (!is.null(children[[old_index + 1]])) {
      for (child in children[[old_index + 1]]) {
        queue <- c(queue, list(child))
      }
    }
  }

  # Step 3: Create a new parent array using the mapping
  new_parent_array <- sapply(parent_array, function(parent) {
    if (parent == -1) {
      return(-1)
    } else {
      return(old_to_new[[as.character(parent)]])
    }
  })

  #artifact from python
  new_parent_array <- new_parent_array + 1

  # Order the new_parent_array array
  new_parent_array_ordered <- sort(new_parent_array)
  return(new_parent_array_ordered)
}



parentArrayToEdgeMatrix <- function(parentArray) {
  edgeMatrix <- matrix(nrow=0, ncol=2)  # Initialize empty edge matrix
  for (child in 1:length(parentArray)) {
    parent <- parentArray[child]
    if (parent != 0) {  # If the node has a parent
      edgeMatrix <- rbind(edgeMatrix, c(parent, child))  # Add edge to matrix
    }
  }
  return(edgeMatrix)
}

create_hdp_tree <- function(ape_tree, tricl_mut_counts){

  #TODO sanity check that the tree and the mutation counts match
  if (length(ape_tree[["tip.label"]]) != nrow(tricl_mut_counts)){
    stop("length(ape_tree[[\"tip.label\"]]) != nrow(tricl_mut_counts)")
  }

  ape_tree_edges <- ape_tree[["edge"]]
  ape_tree_edges_df <- ape_tree_edges %>%
    as.data.frame %>%
    rename(from = V1, to = V2)  %>%
    arrange(to)

  ape_tree_edges_df$to_name <- c(ape_tree$tip.label,
                                 paste0("n", seq_len(ape_tree$Nnode-1)))

  ppindex <- edge_matrix_to_parent_array(ape_tree_edges_df)

  cpindex <- 1 + ppindex #TODO
  hh <- rep(1, ncol(tricl_mut_counts)) #TODO
  alphaa <- rep(1, 1 + max(ppindex)) #TODO
  alphab <- rep(1, 1 + max(ppindex)) #TODO

  hdp_mut_tree <- hdp_init(ppindex = ppindex, # index of parental nodes
                           cpindex = cpindex, # index of concentration param
                           hh=hh, # prior is uniform over the 96 mutation categories
                           alphaa=alphaa, # shape hyperparams for five different CPs
                           alphab=alphab) # rate hyperparams for five different CPs

  names_ppindex <- names(ppindex)
  tricl_mut_counts_cell <- rownames(tricl_mut_counts)#tricl_mut_counts[,'Cell']
  dpindex <- match(tricl_mut_counts_cell,
                   names_ppindex)

  #cnsignatures_input_filtered_cin_quant <- as.data.frame(cnsignatures_input_filtered_cin_quant)
  # add data to leaf nodes (one per cancer sample, in row order of mut_count)
  hdp_mut_tree <- hdp_setdata(hdp_mut_tree,
                              dpindex = dpindex, # index of nodes to add data to
                              data = tricl_mut_counts) # input data (mutation counts, sample rows match up with specified dpindex)
  return(hdp_mut_tree)
}
