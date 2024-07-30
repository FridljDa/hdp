#Daniel_experimental
extract_signature_matrix <- function(mut_example_multi, chains_selection = seq_along(chains(mut_example_multi)), plot_bars = FALSE){

  mut_example_multi_copy <- mut_example_multi

  #if(chains_selection != seq_along(chains(mut_example_multi))){
    ##---remove chains to be ignored----
    mut_example_multi_copy@numcomp <- integer()
    mut_example_multi_copy@prop.ex <- double()
    mut_example_multi_copy@comp_cos_merge <- double()
    mut_example_multi_copy@comp_categ_counts <- list()
    mut_example_multi_copy@comp_dp_counts <- list()
    mut_example_multi_copy@comp_categ_distn <- list()
    mut_example_multi_copy@comp_dp_distn <- list()

    chains_full = seq_along(chains(mut_example_multi_copy))
    chains_selection = 1

    #
    chain_diff <- setdiff(chains_full, chains_selection)
    for(chain_diff_i in chain_diff){
      mut_example_multi_copy@chains[[chain_diff_i]] <- NULL
    }
  #}

  ##----- Extract the components of the model----
  mut_example_multi_copy <- hdp_extract_components(mut_example_multi_copy)
  categ_counts <- comp_categ_counts(mut_example_multi_copy)

  categ_counts_means <- lapply(categ_counts, colMeans)
  categ_counts_means_bound <- do.call(rbind, categ_counts_means)

  #colnames(categ_counts_means_bound) <- colnames(vcf_tricl_mut_counts)
  #rownames(categ_counts_means_bound) <- paste0("Signature.tree.", seq_len(nrow(categ_counts_means_bound)))

  categ_counts_means_bound <- categ_counts_means_bound/rowSums(categ_counts_means_bound)

  if(plot_bars){
    par(mfrow=c(1,3))
    for (i in 1:3){
      barplot(categ_counts_means_bound[i,], xpd=F, las=1)
    }
  }

  return(categ_counts_means_bound)
}

extract_signature_matrix_chains <- function(mut_example_multi){
  signatures_hdp_chains <- seq_along(chains(mut_example_multi)) %>%
    lapply(function(i){
      signatures_hdp_chain_i <- extract_signature_matrix(mut_example_multi, chains_selection = i, plot_bars = FALSE)

      colnames(signatures_hdp_chain_i) <- colnames(mut_count)
      rownames(signatures_hdp_chain_i) <- paste0("HDP.Chain", i,".Signature.", seq_len(nrow(signatures_hdp_chain_i)))

      return(signatures_hdp_chain_i)
    })


  signatures_hdp_chains <- do.call(rbind, signatures_hdp_chains)
  return(signatures_hdp_chains)
}

plot_cosine_sim <- function(signatures, filename, width = 800, height = 800){
  signatures_matrix <- as.matrix(signatures)

  cosine_similarity_matrix <- lsa::cosine(t(signatures_matrix))

  png(filename = filename, width = width, height = height)

  # Step 4: Plot the cosine similarity matrix
  heatmap(cosine_similarity_matrix, Rowv=NA, Colv="Rowv", col = cm.colors(256), scale="column", margins=c(5,10))

  dev.off()

  # Step 4: Plot the cosine similarity matrix
  heatmap(cosine_similarity_matrix, Rowv=NA, Colv="Rowv", col = cm.colors(256), scale="column", margins=c(5,10))
}
