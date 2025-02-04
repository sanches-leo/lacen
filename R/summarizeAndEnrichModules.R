#' Summarize and Enrich Modules
#'
#' Constructs the co-expression network, performs module enrichment analysis, and summarizes the data.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param maxBlockSize Maximum number of genes processed simultaneously.
#' @param TOMType Type of Topological Overlap Matrix.
#' @param minModuleSize Minimum number of genes per module.
#' @param reassignThreshold p-value threshold for reassigning genes between modules.
#' @param mergeCutHeight Dendrogram cut height for merging modules.
#' @param pamRespectsDendro Logical. If TRUE, ensures PAM respects dendrogram structure.
#' @param corType Type of correlation to use: "pearson" or "bicor".
#' @param userThreshold p-value threshold for enrichment significance.
#' @param sources Functional term sources for enrichment analysis (e.g., "GO:BP").
#' @param organism Organism name (e.g., "hsapiens").
#' @param orgdb Bioconductor OrgDb package (e.g., "org.Hs.eg.db").
#' @param reducedTermsThreshold Similarity threshold for reducing enrichment terms.
#' @param filename Filename to save the enrichment graph as PNG.
#' @param ... Additional parameters for internal functions.
#'
#' @return A 'lacen' S3 object with updated `summdf`, `rrvgolist`, and `TOM`.
#' @export
#' @import WGCNA gprofiler2 rrvgo
summarizeAndEnrichModules <- function(lacenObject,
                                       maxBlockSize = 20000,
                                       TOMType = "unsigned",
                                       minModuleSize = 30,
                                       reassignThreshold = 0,
                                       mergeCutHeight = 0.3,
                                       pamRespectsDendro = FALSE,
                                       corType = "bicor",
                                       userThreshold = 0.05,
                                       sources = "BP",
                                       organism = "hsapiens",
                                       orgdb = "org.Hs.eg.db",
                                       reducedTermsThreshold = 0.7,
                                       filename = "5_enrichedgraph.png",
                                       ...) {
  UseMethod("summarizeAndEnrichModules")
}

#' @export
summarizeAndEnrichModules.lacen <- function(lacenObject,
                                             maxBlockSize = 20000,
                                             TOMType = "unsigned",
                                             minModuleSize = 30,
                                             reassignThreshold = 0,
                                             mergeCutHeight = 0.3,
                                             pamRespectsDendro = FALSE,
                                             corType = "bicor",
                                             userThreshold = 0.05,
                                             sources = "BP",
                                             organism = "hsapiens",
                                             orgdb = "org.Hs.eg.db",
                                             reducedTermsThreshold = 0.7,
                                             filename = "5_enrichedgraph.png",
                                             ...) {
  # Extract necessary data from lacenObject
  annotation_data <- lacenObject$annotationData
  expression_data <- lacenObject$datExpression
  nc_annotation <- lacenObject$ncAnnotation
  expr_data <- lacenObject$datExpr
  soft_power <- lacenObject$indicePower
  traits <- lacenObject$datTraits$Trait
  cut_bootstrap <- lacenObject$cutBootstrap
  bootstrap_stability <- lacenObject$bootstrapStability
  ontology <- sources
  
  # Capture additional arguments
  additional_args <- list(...)
  if (length(additional_args) > 0) {
    for (arg_name in names(additional_args)) {
      assign(arg_name, additional_args[[arg_name]])
    }
  }
  
  # Default values for additional parameters if not provided
  if (!exists("log")) log <- FALSE
  if (!exists("log_path")) log_path <- "log.txt"
  if (!exists("modPath")) modPath <- FALSE
  
  # Parameter validation
  if (isTRUE(log)) {
    cat("Starting summarizeAndEnrichModules function\n", file = log_path, append = TRUE)
  }
  
  # Internal function to validate parameters
  validate_parameter <- function(x, expected_class, param_name) {
    if (missing(x) || is.null(x) || length(x) == 0) {
      stop(paste0("Parameter '", param_name, "' is missing or NULL."))
    }
    if (expected_class == "list" && !is.list(x)) {
      stop(paste0("Parameter '", param_name, "' must be a list."))
    }
    if (expected_class == "character" && !is.character(x)) {
      stop(paste0("Parameter '", param_name, "' must be a character string."))
    }
    if (expected_class == "numeric" && !is.numeric(x)) {
      stop(paste0("Parameter '", param_name, "' must be numeric."))
    }
    if (expected_class == "logical" && !is.logical(x)) {
      stop(paste0("Parameter '", param_name, "' must be logical (TRUE/FALSE)."))
    }
  }
  
  # Validate essential parameters
  validate_parameter(annotation_data, "list", "annotationData")
  validate_parameter(nc_annotation, "list", "ncAnnotation")
  validate_parameter(expression_data, "matrix", "datExpression")
  validate_parameter(expr_data, "matrix", "datExpr")
  validate_parameter(soft_power, "numeric", "indicePower")
  
  if (!all(c("gene_id", "gene_name") %in% colnames(annotation_data))) {
    stop("annotationData must contain 'gene_id' and 'gene_name' columns.")
  }
  
  if (!all(c("gene_id", "gene_name") %in% colnames(nc_annotation))) {
    stop("ncAnnotation must contain 'gene_id' and 'gene_name' columns.")
  }
  
  if (!all(c("gene_id", "log2FC", "pval") %in% colnames(lacenObject$datExpression))) {
    stop("datExpression must contain 'gene_id', 'log2FC', and 'pval' columns.")
  }
  
  # Check organism consistency
  if (organism != "hsapiens") {
    warning("Organism is not 'hsapiens'. Ensure that 'orgdb' matches the specified organism.")
  }
  
  # Map ontology sources
  ontology_map <- list(
    BP = "GO:BP",
    MF = "GO:MF",
    CC = "GO:CC"
  )
  
  if (!ontology %in% names(ontology_map)) {
    stop("Ontology must be one of 'BP', 'MF', or 'CC'.")
  }
  
  sources <- ontology_map[[ontology]]
  
  if (isTRUE(log)) {
    cat("Loading functions\n", file = log_path, append = TRUE)
  }
  
  # Internal function to summarize data
  summarize_data <- function(all_degrees, modules, bootstrap_stability, cut_bootstrap, 
                             nc_annotation, coding_genes, differential_expression) {
    # Assign module information
    module_info <- modules
    
    # Apply bootstrap stability cutoff
    if(!isFALSE(bootstrap_stability)){
      stability_filter <- as.vector(bootstrap_stability > cut_bootstrap/100)
    } else{
      stability_filter <- TRUE
    }
    
    # Combine gene names with coding status
    nc_annotation$is_nc <- TRUE
    coding_genes$is_nc <- FALSE
    combined_gene_names <- rbind(nc_annotation, coding_genes)
    combined_gene_names <- combined_gene_names[combined_gene_names$gene_id %in% rownames(all_degrees), ]
    combined_gene_names <- combined_gene_names[match(rownames(all_degrees), combined_gene_names$gene_id), ]
    
    # Merge differential expression data
    deg_data <- differential_expression[match(rownames(all_degrees), differential_expression$gene_id), ]
    
    # Compile summarized dataframe
    summary_df <- cbind(
      combined_gene_names,
      module = module_info,
      cutBootstrap = stability_filter,
      all_degrees,
      deg_data
    )
    
    return(summary_df)
  }
  
  # Internal function to perform enrichment analysis using gprofiler
  perform_enrichment <- function(summary_df, user_thresh, sources, organism) {
    background_genes <- summary_df[summary_df$cutBootstrap, ]$gene_id
    
    enrichment_results <- list()
    
    for(mod in sort(unique(summary_df$module))) {
      if (mod != 0) {
        query_genes <- summary_df$gene_id[summary_df$module == mod & summary_df$cutBootstrap]
        if (length(query_genes) > 0) {
          result <- tryCatch({
            gprofiler2::gost(
              query = query_genes,
              evcodes = TRUE,
              multi_query = FALSE,
              ordered_query = TRUE,
              user_threshold = user_thresh,
              custom_bg = background_genes,
              sources = sources,
              organism = organism
            )
          }, error = function(e) {
            warning(paste("gprofiler2::gost failed for module", mod, ":", e$message))
            return(NULL)
          })
          enrichment_results[[as.character(mod)]] <- result
        }
      }
    }
    
    # Remove NULL results
    enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]
    return(enrichment_results)
  }
  
  # Internal function to reduce enrichment terms using rrvgo
  reduce_enrichment_terms <- function(enrichment_list, summary_df, expr_data, power, threshold, org_db, ontology) {
    
    adjacency_matrix <- abs(WGCNA::cor(expr_data))^soft_power
    
    ncData <- summary_df[summary_df$cutBootstrap & summary_df$is_nc, c("gene_id", "gene_name")]
    
    
    convert_id_nc <- function(ncData, nc_id){
      if(length(nc_id) == 1){
        if(nc_id %in% ncData$gene_id){
          return(ncData[ncData$gene_id == nc_id, "gene_name"])
        } else {
          return(nc_id)
        }
      } else {
        result <- c()
        for(i in nc_id){
          if(i %in% ncData$gene_id){
            result <- c(result, ncData[ncData$gene_id == i, "gene_name"])
          } else {
            result <- c(result, i)
          }
        }
        return(result)
      }
    }
    
    
    reduced_list <- list()
    
    for (mod in names(enrichment_list)) {
      enrichment_data <- enrichment_list[[mod]]$result
      if (nrow(enrichment_data) > 2) {
        # Calculate similarity matrix
        orgdb <- org_db
        ont <- ontology
        sim_matrix <- rrvgo::calculateSimMatrix(enrichment_data$term_id, 
                                                orgdb = org_db, 
                                                semdata = GOSemSim::godata(annoDb = orgdb, ont = ont, keytype = "ENTREZID"),
                                                ont = ontology, 
                                                method = "Rel")
        
        # Assign scores based on p-values
        scores <- setNames(-log10(enrichment_data$p_value), enrichment_data$term_id)
        
        # Reduce the similarity matrix
        reduced_terms <- rrvgo::reduceSimMatrix(sim_matrix, 
                                                scores, 
                                                threshold = threshold, 
                                                orgdb = org_db)
        
        clustered_genes <- clustered_genes[order(clustered_genes$cluster), ]
        
        # Aggregate genes per cluster
        submodules <- lapply(unique(clustered_genes$cluster), function(cluster_id) {
          intersect_genes <- unique(
            unlist(
              strsplit(
                clustered_genes$intersection[clustered_genes$cluster == cluster_id], ",")
            )
          )
          return(intersect_genes)
        })
        
        
        # Aggregate genes per cluster
        submodules <- lapply(unique(clustered_genes$cluster), function(cluster_id) {
          intersect_genes <- unique(unlist(strsplit(clustered_genes$intersection[clustered_genes$cluster == cluster_id], ",")))
          return(intersect_genes)
        })
        
        #Genes is noncoding
        genes_in_module <- summary_df[summary_df$module == mod, ]$gene_id
        nc_in_module <- ncData[ncData$gene_id %in% genes_in_module, "gene_id"]
        
        #loop in submodule?
        mod_adj <- adjacency_matrix[summary_df$gene_id %in% unique(unlist(submodules)), summary_df$gene_id %in% nc_in_module]
        
        
        #get module kwithin median
        median_kwithin_submod <- vector()
        
        
        if(!is.vector(mod_adj) & dim(as.data.frame(mod_adj))[1] > 1 & dim(as.data.frame(mod_adj))[2] > 1){
          
          for(submodule in seq(1, length(submodules))){
            submod_adj <- adjacency_matrix[summary_df$gene_id %in% submodules[[submodule]],summary_df$gene_id %in% submodules[[submodule]]]
            kwithin_submodule_genes <- apply(submod_adj, 2, sum) - 1
            median_kwithin <- stats::median(kwithin_submodule_genes)
            median_kwithin_submod[as.character(submodule)] <- median_kwithin
          }
          
          #for each submodule, subset the module adj dataframe and find the lncrnas submodule kwithin/submodule
          nc_kwithin_submod_list <- list()
          
          for(submodule in seq(1, length(submodules))){
            submod_adj <- mod_adj[rownames(mod_adj) %in% submodules[[submodule]],]
            nc_kwithin_submod <- apply(submod_adj, 2, sum)
            nc_kwithin_submod_list[[submodule]] <- nc_kwithin_submod
          }
          
          nc_kwithin_submod_df <- as.data.frame(nc_kwithin_submod_list)
          colnames(nc_kwithin_submod_df) <- seq(1, length(submodules))
          nc_kwithin_submod_df <- t(nc_kwithin_submod_df)
          
          #compare each kwithin(median_kwithin_submod) with
          #the median connectivity of each submodule (nc_kwithin_submod_df)
          lnc_over_median_kwithin <- list()
          
          for(submodule in seq(1, length(submodules))){
            lnc_over_median <- which(nc_kwithin_submod_df[submodule,] > median_kwithin_submod[as.character(submodule)])
            if(!length(names(lnc_over_median)) == 0){
              lnc_over_median_kwithin[[as.character(submodule)]] <- names(lnc_over_median)
            } else {
              lnc_over_median_kwithin[[as.character(submodule)]] <- NA
            }
            
          }
          
          pre_reduced_term <- lnc_over_median_kwithin[reduced_terms$cluster]
          DE_genes <- summary_df[abs(summary_df$log2FC) > 1 & summary_df$pval < 0.05, "gene_id"]
          DE_genes <- DE_genes[!is.na(DE_genes)]
          
          n_lnc <- lapply(pre_reduced_term, function(x) ifelse(is.na(x), 0, length(x)))
          n_lnc <- lapply(n_lnc, unique)
          id_lnc <- lapply(pre_reduced_term, paste, collapse = ", ")
          name_lnc <- lapply(pre_reduced_term, function(x) convert_id_nc(ncData, x))
          name_lnc <- lapply(name_lnc, paste, collapse = ", ")
          n_de <- lapply(pre_reduced_term, function(x) sum(x %in% DE_genes))
          n_de <- lapply(n_de, function(x) ifelse(x == 0 & length(x) == 1, 0, x))
          ID_de0 <- lapply(pre_reduced_term, function(x) x[x %in% DE_genes])
          ID_de <- lapply(ID_de0, function(x) ifelse(length(x) == 0, NA, paste(x, collapse = ", ")))
          name_de <- lapply(ID_de0, function(x) convert_id_nc(ncData, x))
          name_de <- lapply(name_de, function(x) ifelse(is.null(x), NA, paste(x, collapse = ", ")))
          
          reduced_terms$number_lnc_highly_connected <- unlist(n_lnc)
          reduced_terms$id_lnc_highly_connected <- id_lnc
          reduced_terms$name_lnc_highly_connected <- name_lnc
          reduced_terms$number_lnc_DE <- unlist(n_de)
          reduced_terms$id_lnc_DE <- ID_de
          reduced_terms$name_lnc_DE <- name_de
          
          
          reduced_terms$parentTerm2 <- reduced_terms$parentTerm
          reduced_terms$parentTerm <- paste(reduced_terms$cluster, " - ", reduced_terms$parentTerm, " (",n_de , "/", n_lnc, ")", sep = "" )
          
          # Store the reduced terms
          reduced_list[[mod]] <- reduced_terms
        }
      }
    }
    
    return(reduced_list)
  }
  
  # Internal function to save the enrichment graph as a treemap
  save_enrichment_graph <- function(reduced_enrichment, summary_df, expr_data, traits, filename, mod_path) {
    # Prepare genes to select
    genes_selected <- summary_df$gene_id[summary_df$cutBootstrap]
    
    # Extract expression data for selected genes
    expr_selected <- expr_data[, colnames(expr_data) %in% genes_selected]
    
    # Define layout for multiple treemaps
    num_graphs <- length(reduced_enrichment)
    grid_size <- ceiling(sqrt(num_graphs))
    
    grDevices::png(filename = filename, width = 800 * grid_size, height = 800 * grid_size)
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(grid_size, grid_size)))
    
    counter <- 1
    for (mod in names(reduced_enrichment)) {
      if (counter > num_graphs) break
      module <- mod
      genes_in_module <- summary_df$module == module
      
      # Calculate module eigengenes
      module_eigengenes <- WGCNA::moduleEigengenes(expr_selected, as.numeric(genes_in_module))$eigengenes
      ordered_MEs <- WGCNA::orderMEs(module_eigengenes)
      module_trait_cor <- WGCNA::cor(ordered_MEs, traits, use = "p")
      module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nrow(expr_selected))
      correlation <- module_trait_cor[[2]]
      p_value <- module_trait_pval[[2]]
      
      # Plot treemap using rrvgo
      rrvgo::treemapPlot(reduced_enrichment[[mod]], size = "score",
                         title = paste0("Module ", module, 
                                        " (corr: ", round(correlation, 2), 
                                        ", pValue: ", format(p_value, scientific = TRUE, digits = 2), ")"),
                         vp = grid::viewport(layout.pos.row = ceiling(counter / grid_size),
                                             layout.pos.col = counter %% grid_size),
                         fontsize.title = 22,
                         fontsize.labels = 14)
      counter <- counter + 1
    }
    
    grDevices::dev.off()
    
    # Optionally save individual module plots
    if (!isFALSE(mod_path)) {
      counter <- 1
      for (mod in names(reduced_enrichment)) {
        if (counter > num_graphs) break
        module <- mod
        genes_in_module <- summary_df$module == module
        
        # Calculate module eigengenes
        module_eigengenes <- WGCNA::moduleEigengenes(expr_selected, as.numeric(genes_in_module))$eigengenes
        ordered_MEs <- WGCNA::orderMEs(module_eigengenes)
        module_trait_cor <- WGCNA::cor(ordered_MEs, traits, use = "p")
        module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nrow(expr_selected))
        correlation <- module_trait_cor[[2]]
        p_value <- module_trait_pval[[2]]
        
        # Define file path for the module plot
        module_filename <- file.path(mod_path, paste0("module_", module, ".png"))
        
        # Save treemap for the module
        grDevices::png(filename = module_filename, width = 800, height = 800)
        rrvgo::treemapPlot(reduced_enrichment[[mod]], size = "score",
                           title = paste0("Module ", module, 
                                          " (corr: ", round(correlation, 2), 
                                          ", pValue: ", format(p_value, scientific = TRUE, digits = 2), ")"),
                           fontsize.title = 22,
                           fontsize.labels = 14)
        grDevices::dev.off()
        counter <- counter + 1
      }
    }
  }
  
  # Internal function to generate Topological Overlap Matrix (TOM)
  generate_TOM <- function(expr_data, power, summary_df) {
    tom_matrix <- WGCNA::TOMsimilarityFromExpr(expr_data, power = power, TOMType = "unsigned")
    rownames(tom_matrix) <- colnames(expr_data)
    colnames(tom_matrix) <- colnames(expr_data)
    # Exclude genes not assigned to any module
    tom_filtered <- tom_matrix[!summary_df$module == 0, !summary_df$module == 0]
    return(tom_filtered)
  }
  
  if (isTRUE(log)) {
    cat("Making the network\n", file = log_path, append = TRUE)
  }
  
  # Construct the co-expression network and identify modules
  network <- WGCNA::blockwiseModules(
    expr_data, 
    power = soft_power, 
    maxBlockSize = maxBlockSize,
    TOMType = TOMType, 
    minModuleSize = minModuleSize,
    reassignThreshold = reassignThreshold, 
    mergeCutHeight = mergeCutHeight,
    numericLabels = TRUE, 
    pamRespectsDendro = pamRespectsDendro,
    verbose = 3, 
    corType = corType
  )
  
  modules <- network$colors
  
  if (isTRUE(log)) {
    cat("Getting connectivities\n", file = log_path, append = TRUE)
  }
  
  # Calculate intramodular connectivity
  intramodular_conn <- WGCNA::intramodularConnectivity.fromExpr(
    datExpr = expr_data, 
    power = soft_power, 
    colors = modules
  )
  intramodular_conn$gene_id <- colnames(expr_data)
  rownames(intramodular_conn) <- colnames(expr_data)
  
  # Extract coding gene names
  coding_gene_names <- annotation_data[!annotation_data$gene_id %in% nc_annotation$gene_id, ]
  
  if (isTRUE(log)) {
    cat("Summarizing the data\n", file = log_path, append = TRUE)
  }
  
  # Summarize the data without enrichment
  summary_df <- summarize_data(
    all_degrees = intramodular_conn,
    modules = modules,
    bootstrap_stability = bootstrap_stability,
    cut_bootstrap = cut_bootstrap,
    nc_annotation = nc_annotation,
    coding_genes = coding_gene_names,
    differential_expression = expression_data
  )
  
  if (isTRUE(log)) {
    cat("Enriching the modules\n", file = log_path, append = TRUE)
  }
  
  # Perform enrichment analysis
  enrichment_results <- perform_enrichment(summary_df, userThreshold, sources, organism)
  
  if (isTRUE(log)) {
    cat("Reducing the enriched modules\n", file = log_path, append = TRUE)
  }
  
  # Reduce enrichment terms
  reduced_enrichment <- reduce_enrichment(
    enrichment_results, 
    summary_df, 
    expr_data, 
    soft_power, 
    reducedTermsThreshold, 
    orgdb, 
    ontology
  )
  
  if (length(reduced_enrichment) > 0) {
    if (isTRUE(log)) {
      cat("Generating plot\n", file = log_path, append = TRUE)
    }
    
    # Save the enrichment treemap plots
    save_enrichment_graph(reduced_enrichment, summary_df, expr_data, traits, filename, modPath)
    
    if (isTRUE(log)) {
      cat("Final Summarizing\n", file = log_path, append = TRUE)
    }
    
    # Optionally, summarize submodules (not fully implemented here)
    # summary_df <- summarizeSubmodules(summary_df, reduced_enrichment, enrichment_results)
    
    if (isTRUE(log)) {
      cat("Preparing Topological Overlap Matrix\n", file = log_path, append = TRUE)
    }
    
    # Generate the Topological Overlap Matrix
    tom_matrix <- generate_TOM(expr_data, soft_power, summary_df)
    
    # Update the lacenObject with enrichment results
    lacenObject$summdf <- summary_df
    lacenObject$rrvgolist <- reduced_enrichment
    lacenObject$TOM <- tom_matrix
    
    return(lacenObject)
  } else {
    stop("No module could be enriched. Consider increasing the number of transcripts or lowering filter parameters.")
  }
  
  if (isTRUE(log)) {
    cat("Done!\n", file = log_path, append = TRUE)
  }
}