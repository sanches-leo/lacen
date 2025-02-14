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
                           pathModGroupsPlot = "3_modgroupsnew.png",
                           pathStabilityPlot = "4_stability_bootstrap.png",
                           nThreads = 4,
                           ...) {
  # Create a generic function to allow method dispatch based on lacenObject class
  UseMethod("lacenBootstrap")
}

#' @export
lacenBootstrap.lacen <- function(lacenObject,
                                 numberOfIterations = 100,
                                 maxBlockSize = 50000,
                                 pathModGroupsPlot = "3_modgroups.png",
                                 pathStabilityPlot = "4_stability_bootstrap.png",
                                 nThreads = 4,
                                 ...){
  
  # Ensure that WGCNA uses its own correlation function instead of base cor()
  cor <- WGCNA::cor
  
  # Extract the expression data and soft-thresholding power from the lacenObject
  expr_data <- lacenObject$datExpr
  soft_power <- lacenObject$indicePower
  
  # Capture any additional arguments passed via ... and assign them as separate variables
  additional_args <- list(...)
  if (length(additional_args) > 0) {
    for (arg_name in names(additional_args)) {
      assign(arg_name, additional_args[[arg_name]])
    }
  }
  
  # Determine whether to use parallel processing based on nThreads.
  # Also, configure the WGCNA threading accordingly.
  if (!exists("parallel")) {
    if (nThreads < 8) {
      # For fewer threads, disable WGCNA internal multithreading
      WGCNA::disableWGCNAThreads()
      use_parallel <- FALSE
    } else {
      # For sufficient threads, enable WGCNA multithreading using 4 threads
      WGCNA::enableWGCNAThreads(4)
      use_parallel <- TRUE
      # Set the number of parallel workers based on the available threads
      num_parallel <- floor(nThreads / 4)
    }
  }
  
  # If cutBootstrap is not provided in the additional arguments, set its default value
  if (!exists("cutBootstrap")) {
    cutBootstrap <- FALSE
  }
  
  # ------------------------------------------------------------------------------
  # Nested function: makeBootstrap
  # This function performs the bootstrapping iterations, running network construction
  # on subsets of genes (by leaving out a subset for each iteration), either in parallel
  # or sequentially.
  # ------------------------------------------------------------------------------
  makeBootstrap <- function(numberOfIterations,
                            datExpr,
                            indicePower,
                            maxBlockSize,
                            parallel,
                            nparallel,
                            WGCNAThreads){
    
    # Set WGCNA threading depending on the input parameter WGCNAThreads.
    # If threading is effectively disabled (FALSE or 1), then disable it.
    if(WGCNAThreads == FALSE | WGCNAThreads == 1){
      WGCNA::disableWGCNAThreads()
    } else {
      WGCNA::enableWGCNAThreads(WGCNAThreads)
    }
    
    # If parallel processing is enabled, use foreach for iterations
    if(isTRUE(parallel)){
      # Register parallel backend with the specified number of workers
      doParallel::registerDoParallel(nparallel)
      
      # Initialize a data frame to store bootstrap results:
      # rows represent each bootstrap iteration, columns correspond to genes.
      bootstrap = as.data.frame(matrix(data = NA, nrow = numberOfIterations, ncol = ncol(datExpr)))
      colnames(bootstrap) = colnames(datExpr)
      
      # Run blockwiseModules on the full dataset to obtain the reference module assignment
      net_no_bootstrap <- WGCNA::blockwiseModules(datExpr, power = indicePower, maxBlockSize = maxBlockSize,
                                                  TOMType = "unsigned", minModuleSize = 30,
                                                  reassignThreshold = 0, mergeCutHeight = 0.3,
                                                  numericLabels = TRUE, pamRespectsDendro = FALSE,
                                                  verbose = 0)
      # Save the module assignment from the full dataset as the first line
      first_line <- net_no_bootstrap$colors
      
      # Create a vector that randomly assigns each gene an iteration number
      # in which it will be left out from the network construction.
      this_gene_wont_be_in_iteratiom_number = sample(rep(1:numberOfIterations, 
                                                         ceiling(dim(datExpr)[[2]]/numberOfIterations)), 
                                                     dim(datExpr)[[2]], replace = FALSE)
      
      # Use foreach to iterate over the bootstrap iterations in parallel.
      boot_foreach <- foreach::foreach(i = 1:(numberOfIterations/1)) %dopar% {
        # For each iteration, exclude genes assigned to iteration i.
        net = WGCNA::blockwiseModules(datExpr[,this_gene_wont_be_in_iteratiom_number != i], 
                                      power = indicePower, maxBlockSize = maxBlockSize,
                                      TOMType = "unsigned", minModuleSize = 30,
                                      reassignThreshold = 0, mergeCutHeight = 0.3,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      verbose = 0)
        # Initialize an empty vector to store module assignments for all genes
        next_line = c()
        counter = 1
        # Reconstruct a full vector for all genes:
        # - For genes included in the current bootstrap, assign the computed module color.
        # - For genes left out, assign NA.
        for(j in 1:dim(datExpr)[[2]]){
          if(this_gene_wont_be_in_iteratiom_number[j] != i){
            next_line = c(next_line, net$colors[counter])
            counter = counter + 1
          }
          else{
            next_line = c(next_line, NA)
          }
        }
        # Return the reconstructed module assignment vector for iteration i
        return(next_line)
      }
      # Combine the results from parallel iterations and transpose for proper formatting
      boot_foreach <- t(as.data.frame(boot_foreach))
      # Prepend the full dataset module assignment (first_line) as the first row
      boot_foreach <- rbind(first_line, boot_foreach)
      row.names(boot_foreach) <- 0:(nrow(boot_foreach) - 1)
      return(boot_foreach)
    } else {
      # Sequential processing if parallel is not enabled
      
      # Initialize bootstrap data frame similar to the parallel case
      bootstrap = as.data.frame(matrix(data = NA, nrow = numberOfIterations, ncol = ncol(datExpr)))
      colnames(bootstrap) = colnames(datExpr)
      
      # Run network construction on the full dataset
      net_no_bootstrap <- WGCNA::blockwiseModules(datExpr, power = indicePower, maxBlockSize = maxBlockSize,
                                                  TOMType = "unsigned", minModuleSize = 30,
                                                  reassignThreshold = 0, mergeCutHeight = 0.3,
                                                  numericLabels = TRUE, pamRespectsDendro = FALSE,
                                                  verbose = 0)
      bootstrap[1,] <- net_no_bootstrap$colors
      
      # Create the assignment vector indicating which genes to leave out in each iteration
      this_gene_wont_be_in_iteratiom_number = sample(rep(1:numberOfIterations, 
                                                         ceiling(dim(datExpr)[[2]]/numberOfIterations)), 
                                                     dim(datExpr)[[2]], replace = FALSE)
      # Loop through each bootstrap iteration
      for(i in 1:(numberOfIterations/1)){
        # Exclude genes assigned to the current iteration i
        net = WGCNA::blockwiseModules(datExpr[,this_gene_wont_be_in_iteratiom_number != i], 
                                      power = indicePower, maxBlockSize = maxBlockSize,
                                      TOMType = "unsigned", minModuleSize = 30,
                                      reassignThreshold = 0, mergeCutHeight = 0.3,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      verbose = 0)
        next_line = c()
        counter = 1
        # Reconstruct the complete vector of module assignments for all genes
        for(j in 1:dim(datExpr)[[2]]){
          if(this_gene_wont_be_in_iteratiom_number[j] != i){
            next_line = c(next_line, net$colors[counter])
            counter = counter + 1
          }
          else{
            next_line = c(next_line, NA)
          }
        }
        # Save the iteration result (note: row 1 is already used by the full dataset)
        bootstrap[i+1,] <- next_line
      }
      return(bootstrap)
    }
  }
  
  # ------------------------------------------------------------------------------
  # Nested function: moduleStability
  # This function calculates the stability of the modules across the bootstrap iterations.
  # It constructs a list for each bootstrap iteration mapping each module to its genes,
  # compares the modules across iterations, and produces a boxplot to visualize stability.
  # ------------------------------------------------------------------------------
  moduleStability <- function(bootstrap,
                              pathModGroupsPlot){
    # Internal helper function to find the most similar module (vector) in a list
    # compared to the given vector, based on the number of shared genes.
    most_similar_vector <- function(vector, list_of_vectors){
      sim <- 0
      most <- NULL
      for(vec in names(list_of_vectors)){
        if(sum(vector %in% list_of_vectors[[vec]]) > sim){
          sim <- sum(vector %in% list_of_vectors[[vec]])
          most <- vec
        }
      }
      return(most)
    }
    print("Creating Module Stability Matrix (1/2)")
    
    # Build a list of module compositions for each bootstrap repetition.
    bootlist <- list()
    counter <- 0
    totcounter <- nrow(bootstrap)
    for(rep in 1:nrow(bootstrap)){
      replist <- list()
      # For each unique module label in the current iteration,
      # record the indices of the genes assigned to that module.
      for(mod in unique(unlist(bootstrap[rep,]))){
        if(!is.na(mod)){
          replist[[as.character(mod)]] <- which(bootstrap[rep,] == mod)
        }
      }
      bootlist[[as.character(rep)]] <- replist
      counter <- counter + 1
    }
    
    # Create a mapping (modGroups) from each module in the reference (first iteration)
    # to a vector of module labels (one per bootstrap iteration) that are the most similar.
    modGroups <- list()
    for(mod in names(bootlist[[1]])){
      modlist <- mod  # start with the module label from the first iteration
      for(rep in bootlist[2:length(bootlist)]){
        # Find the module label in the current iteration that is most similar to the reference module.
        modlist <- c(modlist, most_similar_vector(bootlist[[1]][[mod]], rep))
        modGroups[[as.character(mod)]] <- modlist
      }
    }
    
    # Calculate stability ratios for each module across bootstrap iterations.
    boxplotlist <- list()
    print("Calculating Module Stability (2/2)")
    counter <- 0
    totcounter <- length(modGroups) * nrow(bootstrap)
    for(mod in names(modGroups)){
      boot <- list()
      for(row in 1:nrow(bootstrap)){
        # For each bootstrap iteration, find indices where the module label equals the expected label
        boot[[row]] <- which(unlist(bootstrap[row,]) == modGroups[[mod]][row])
        counter <- counter + 1
      }
      # Calculate the proportion of overlapping gene indices between the full dataset (first iteration)
      # and each bootstrap iteration (ignoring the first iteration itself).
      boot <- unlist(lapply(boot[2:nrow(bootstrap)], function(x) sum(boot[[1]] %in% x)/length(boot[[1]])))
      boxplotlist[[as.character(mod)]] <- boot
    }
    
    # Plot the stability boxplots and save to a PNG file.
    grDevices::png(filename = pathModGroupsPlot, width = 1200, height = 1200)
    graphics::boxplot(boxplotlist[sort.int(as.numeric(names(boxplotlist)))],
                      main="Module Stability",
                      xlab="Module Number",
                      ylab="Transcripts in the same module",
                      ylim=c(0,1.2),
                      graphics::par(yaxs='i'),
                      graphics::par(mfrow=c(1,1)),
                      cex.lab = 2,
                      cex.axis = 2,
                      cex.main = 2,
                      cex.sub = 2)
    graphics::abline( h = 1, lty = 3)
    grDevices::dev.off()
    return(modGroups)
  }
  
  # ------------------------------------------------------------------------------
  # Nested function: filterStability
  # This function computes a stability ratio for each gene, defined as the proportion
  # of bootstrap iterations in which the gene was assigned to the module indicated by modGroups.
  # ------------------------------------------------------------------------------
  filterStability <- function(modGroups,
                              bootstrap){
    names <- colnames(bootstrap)
    sigs <- c() # To store the stability ratio for each gene
    for(i in 1:ncol(bootstrap)){
      gene <- bootstrap[-c(1), i]    # All bootstrap iterations (excluding the first full-data iteration)
      module <- bootstrap[1, i]       # The reference module assignment from the full dataset
      # Calculate the stability ratio: proportion of iterations where the gene's module
      # equals the expected module label from modGroups.
      sig <- sum(gene == modGroups[[as.character(module)]][-c(1)], na.rm = TRUE) / length(gene)
      sigs <- c(sigs, sig)
    }
    names(sigs) <- colnames(bootstrap)
    return(as.data.frame(sigs))
  }
  
  # ------------------------------------------------------------------------------
  # Nested function: stabilityRatioPlot
  # This function creates a plot of the cumulative distribution of genes by their module
  # stability proportion, optionally marking a threshold (cutBootstrap).
  # ------------------------------------------------------------------------------
  stabilityRatioPlot <- function(bootstrapStability,
                                 cutBootstrap = FALSE,
                                 pathStabilityPlot = pathStabilityPlot){
    sigs <- bootstrapStability
    # Define thresholds from 0 to 1 in increments of 0.1
    seqs <- seq(0, 1, 0.1)
    # For each threshold, count the cumulative number of genes with a stability ratio below that threshold
    cumulativeSig <- unlist(lapply(seqs, function(x) sum(sigs < x)))
    names(cumulativeSig) <- as.character(seqs)
    grDevices::png(filename = pathStabilityPlot, width = 1200, height = 1200)
    plot(x = seqs,
         y = cumulativeSig,
         main = "Cumulative genes distribuction by module stability proportion",
         ylab = "Cumulative number of genes",
         xlab = "Proportions of genes maintained in module",
         cex.lab = 2,
         cex.axis = 2,
         cex.main = 2,
         cex.sub = 2)
    graphics::text(seqs, cumulativeSig, labels = cumulativeSig, cex = 1, pos = 3)
    # If cutBootstrap is provided and is numeric, draw a vertical red line at that threshold.
    if(!isFALSE(cutBootstrap) & !isTRUE(cutBootstrap) & is.numeric(as.numeric(cutBootstrap))){
      graphics::abline(v = cutBootstrap, col = "red")
    }
    grDevices::dev.off()
  }
  
  # ------------------------------------------------------------------------------
  # Main processing: Run bootstrap analysis and update the lacenObject
  # ------------------------------------------------------------------------------
  
  # Generate the bootstrap table containing module assignments per iteration.
  bootstrap <- makeBootstrap(numberOfIterations = numberOfIterations,
                             datExpr = expr_data,
                             indicePower = soft_power,
                             maxBlockSize = maxBlockSize,
                             parallel = parallel,
                             nparallel = nparallel,
                             WGCNAThreads = WGCNAThreads)
  
  # Compute the module stability groups and generate a boxplot saved to pathModGroupsPlot.
  modGroups <- moduleStability(bootstrap = bootstrap,
                               pathModGroupsPlot = pathModGroupsPlot)
  
  # Calculate the stability ratio for each gene based on the modGroups and bootstrap table.
  bootstrapStability <- filterStability(modGroups = modGroups,
                                        bootstrap = bootstrap)
  
  # Generate the stability ratio plot and save it to pathStabilityPlot.
  stabilityRatioPlot(bootstrapStability = bootstrapStability,
                     cutBootstrap = cutBootstrap,
                     pathStabilityPlot = pathStabilityPlot)
  
  # Update the lacenObject with the new bootstrap stability information and bootstrap table.
  lacenObject$bootstrapStability <- bootstrapStability
  lacenObject$bootstrap <- bootstrap
  
  return(lacenObject)
}