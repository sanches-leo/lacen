#' Bootstrap Stability
#'
#' Repeats network construction multiple times with subsets of genes to assess gene stability.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param numberOfIterations Number of bootstrap iterations to perform.
#' @param maxBlockSize Maximum number of genes processed simultaneously.
#' @param pathModGroupsPlot Filename to save the module stability plot as PNG.
#' @param pathStabilityPlot Filename to save the stability ratio plot as PNG.
#' @param nThreads Number of threads to use for parallel processing.
#' @param ... Additional parameters for internal functions.
#'
#' @return A 'lacen' S3 object with updated `bootstrapStability` and `bootstrap` data.
#' @export
#' @import foreach doParallel WGCNA
lacenBootstrap <- function(lacenObject,
                           numberOfIterations = 100,
                           maxBlockSize = 50000,
                           pathModGroupsPlot = "3_modgroups.png",
                           pathStabilityPlot = "4_stability_bootstrap.png",
                           nThreads = 4,
                           ...) {
  UseMethod("lacenBootstrap")
}

#' @export
lacenBootstrap.lacen <- function(lacenObject,
                                  numberOfIterations = 100,
                                  maxBlockSize = 50000,
                                  pathModGroupsPlot = "3_modgroups.png",
                                  pathStabilityPlot = "4_stability_bootstrap.png",
                                  nThreads = 4,
                                  ...) {
  # Extract necessary data
  expr_data <- lacenObject$datExpr
  soft_power <- lacenObject$indicePower
  
  # Capture additional arguments
  additional_args <- list(...)
  if (length(additional_args) > 0) {
    for (arg_name in names(additional_args)) {
      assign(arg_name, additional_args[[arg_name]])
    }
  }
  
  # Determine if parallel processing should be used
  if (!exists("parallel")) {
    if (nThreads < 8) {
      WGCNA::disableWGCNAThreads()
      use_parallel <- FALSE
    } else {
      WGCNA::enableWGCNAThreads(4)
      use_parallel <- TRUE
      num_parallel <- floor(nThreads / 4)
    }
  }
  
  # Initialize cutBootstrap if not set
  if (!exists("cutBootstrap")) {
    cutBootstrap <- FALSE
  }
  
  # Internal function to perform bootstrap iterations
  perform_bootstrap <- function(iterations, expr_data, power, max_block, use_parallel, num_parallel) {
    if (power > 0) {
      status_message <- "Starting hierarchical clustering and module detection."
      message(status_message)
    }
    
    bootstrap_results <- matrix(NA, nrow = iterations + 1, ncol = ncol(expr_data))
    colnames(bootstrap_results) <- colnames(expr_data)
    
    # Initial module detection without bootstrap
    initial_network <- WGCNA::blockwiseModules(
      expr_data, 
      power = power, 
      maxBlockSize = max_block,
      TOMType = "unsigned", 
      minModuleSize = 30,
      reassignThreshold = 0, 
      mergeCutHeight = 0.3,
      numericLabels = TRUE, 
      pamRespectsDendro = FALSE,
      verbose = 0
    )
    bootstrap_results[1, ] <- initial_network$colors
    
    # Assign each gene to a bootstrap iteration for exclusion
    set.seed(123)  # For reproducibility
    exclusion_assignments <- sample(rep(1:iterations, length.out = ncol(expr_data)))
    
    if (use_parallel) {
      doParallel::registerDoParallel(num_parallel)
      
      bootstrap_parallel <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
        # Exclude specified genes for this iteration
        genes_to_include <- exclusion_assignments != i
        subset_expr <- expr_data[, genes_to_include, drop = FALSE]
        
        # Perform module detection on the subset
        subset_network <- WGCNA::blockwiseModules(
          subset_expr, 
          power = power, 
          maxBlockSize = max_block,
          TOMType = "unsigned", 
          minModuleSize = 30,
          reassignThreshold = 0, 
          mergeCutHeight = 0.3,
          numericLabels = TRUE, 
          pamRespectsDendro = FALSE,
          verbose = 0
        )
        
        # Assign modules, marking excluded genes as NA
        modules_subset <- subset_network$colors
        modules_full <- rep(NA, ncol(expr_data))
        modules_full[genes_to_include] <- modules_subset
        return(modules_full)
      }
      
      bootstrap_results[-1, ] <- bootstrap_parallel
    } else {
      # Sequential processing without parallel
      for (i in 1:iterations) {
        genes_to_include <- exclusion_assignments != i
        subset_expr <- expr_data[, genes_to_include, drop = FALSE]
        subset_network <- WGCNA::blockwiseModules(
          subset_expr, 
          power = power, 
          maxBlockSize = max_block,
          TOMType = "unsigned", 
          minModuleSize = 30,
          reassignThreshold = 0, 
          mergeCutHeight = 0.3,
          numericLabels = TRUE, 
          pamRespectsDendro = FALSE,
          verbose = 0
        )
        modules_subset <- subset_network$colors
        modules_full <- rep(NA, ncol(expr_data))
        modules_full[genes_to_include] <- modules_subset
        bootstrap_results[i + 1, ] <- modules_full
      }
    }
    
    return(as.data.frame(bootstrap_results))
  }
  
  # Internal function to assess module stability
  assess_module_stability <- function(bootstrap_df, mod_groups_plot_path) {
    # Helper function to find the most similar module
    find_most_similar <- function(vector, list_of_vectors) {
      similarity_scores <- sapply(list_of_vectors, function(vec) sum(vector %in% vec))
      most_similar <- names(which.max(similarity_scores))
      return(most_similar)
    }
    
    message("Assessing module stability across bootstrap iterations.")
    
    # Initialize list to store module groups
    module_groups <- list()
    initial_modules <- as.character(bootstrap_df[1, ])
    unique_modules <- unique(initial_modules)
    unique_modules <- unique_modules[unique_modules != "0"]
    
    for (mod in unique_modules) {
      module_groups[[mod]] <- mod  # Initialize with original module
      for (i in 2:nrow(bootstrap_df)) {
        current_modules <- unique(as.character(bootstrap_df[i, ]))
        current_modules <- current_modules[current_modules != "0"]
        # Find the most similar module in the current iteration
        similar_mod <- find_most_similar(module_groups[[mod]][i - 1], current_modules)
        module_groups[[mod]] <- c(module_groups[[mod]], similar_mod)
      }
    }
    
    # Calculate stability ratios
    stability_ratios <- sapply(colnames(bootstrap_df), function(gene) {
      initial_mod <- as.character(bootstrap_df[1, gene])
      if (initial_mod == "0") {
        return(NA)
      }
      assigned_mod <- module_groups[[initial_mod]]
      return(mean(bootstrap_df[, gene] == assigned_mod))
    })
    
    # Plot stability ratios
    grDevices::png(filename = pathModGroupsPlot, width = 1200, height = 1200)
    boxplot(stability_ratios ~ as.factor(bootstrap_df[1, ]),
            main = "Module Stability",
            xlab = "Module",
            ylab = "Stability Ratio",
            ylim = c(0, 1),
            cex.lab = 2,
            cex.axis = 2,
            cex.main = 2,
            cex.sub = 2)
    abline(h = 1, lty = 3)
    grDevices::dev.off()
    
    return(stability_ratios)
  }
  
  # Internal function to filter genes based on stability ratios
  filter_genes_by_stability <- function(mod_groups, bootstrap_df) {
    gene_stability <- sapply(colnames(bootstrap_df), function(gene) {
      stability <- sum(bootstrap_df[, gene] == mod_groups[[bootstrap_df[1, gene]]][1:nrow(bootstrap_df)])
      return(stability / nrow(bootstrap_df))
    })
    return(as.data.frame(gene_stability))
  }
  
  # Internal function to plot stability ratios
  plot_stability_ratios <- function(stability_ratios, cutoff, stability_plot_path) {
    sequences <- seq(0, 1, by = 0.1)
    cumulative_counts <- sapply(sequences, function(x) sum(stability_ratios < x, na.rm = TRUE))
    
    grDevices::png(filename = stability_plot_path, width = 1200, height = 1200)
    plot(sequences, cumulative_counts, type = "b",
         main = "Cumulative Gene Distribution by Module Stability",
         xlab = "Stability Proportion",
         ylab = "Cumulative Number of Genes",
         cex.lab = 2,
         cex.axis = 2,
         cex.main = 2,
         cex.sub = 2)
    text(sequences, cumulative_counts, labels = cumulative_counts, pos = 3, cex = 1)
    if (!isFALSE(cutoff) && is.numeric(cutoff) && cutoff >= 0 && cutoff <= 1) {
      abline(v = cutoff, col = "red", lty = 2)
    }
    grDevices::dev.off()
  }
  
  # Perform bootstrap iterations
  bootstrap_df <- perform_bootstrap(
    iterations = numberOfIterations,
    expr_data = expr_data,
    power = soft_power,
    max_block = maxBlockSize,
    use_parallel = use_parallel,
    num_parallel = num_parallel
  )
  
  # Assess module stability
  module_groups <- assess_module_stability(
    bootstrap_df = bootstrap_df,
    mod_groups_plot_path = pathModGroupsPlot
  )
  
  # Calculate stability ratios
  bootstrap_stability <- filter_genes_by_stability(
    mod_groups = module_groups,
    bootstrap_df = bootstrap_df
  )
  
  # Plot stability ratios
  plot_stability_ratios(
    stability_ratios = bootstrap_stability$sigs,
    cutoff = cutBootstrap,
    stability_plot_path = pathStabilityPlot
  )
  
  # Update lacenObject with bootstrap results
  lacenObject$bootstrapStability <- bootstrap_stability
  lacenObject$bootstrap <- bootstrap_df
  
  return(lacenObject)
}