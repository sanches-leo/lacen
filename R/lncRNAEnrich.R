#' lncRNA Enrichment Analysis
#'
#' Performs functional enrichment analysis for a specified lncRNA based on top correlated genes.
#'
#' @param lncName Gene name of the lncRNA for enrichment analysis.
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param nGenes Number of top correlated genes to include in the analysis.
#' @param nHighlight Number of top enriched terms to highlight in the enrichment graph.
#' @param sources Functional term sources for enrichment analysis.
#' @param organism Organism name for functional annotations.
#' @param nGenesNet Number of genes to visualize in the network plot.
#' @param nTerm Number of top enriched terms to display in the enrichment graph.
#' @param lncHighlight Logical. If TRUE, highlights all lncRNAs in the network visualization.
#' @param ... Additional parameters for internal functions.
#'
#' @return None (plots and saves enrichment results and network).
#' @export
#' @import gprofiler2 igraph ggraph scatterpie grid Polychrome
lncRNAEnrich <- function(lncName,
                         lacenObject,
                         nGenes = 100,
                         nHighlight = 10,
                         sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                         organism = "hsapiens",
                         nGenesNet = 10,
                         nTerm = 3,
                         lncHighlight = FALSE,
                         ...) {
  UseMethod("lncRNAEnrich2")
                         }


#' @export
lncRNAEnrich.lacen<- function(lncName,
                         lacenObject,
                         nGenes = 100,
                         nHighlight = 10,
                         sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                         organism = "hsapiens",
                         nGenesNet = 10,
                         nTerm = 3,
                         lncHighlight = FALSE,
                         ...) {
  # Validate input parameters
  if (!is.character(lncName) || length(lncName) != 1) {
    stop("'lncName' must be a single character string.")
  }
  
  if (!inherits(lacenObject, "lacen")) {
    stop("'lacenObject' must be an object of class 'lacen'.")
  }
  
  if (!lncName %in% lacenObject$annotationData$gene_name) {
    stop(paste("lncName '", lncName, "' not found in annotationData.", sep = ""))
  }
  
  if (!is.numeric(nGenes) || length(nGenes) != 1) {
    stop("'nGenes' must be a single numeric value.")
  }
  
  if (!is.numeric(nHighlight) || length(nHighlight) != 1) {
    stop("'nHighlight' must be a single numeric value.")
  }
  
  if (!is.character(sources)) {
    stop("'sources' must be a character vector.")
  }
  
  if (!is.character(organism) || length(organism) != 1) {
    stop("'organism' must be a single character string.")
  }
  
  if (!is.numeric(nGenesNet) || length(nGenesNet) != 1) {
    stop("'nGenesNet' must be a single numeric value.")
  }
  
  if (!is.numeric(nTerm) || length(nTerm) != 1) {
    stop("'nTerm' must be a single numeric value.")
  }
  
  if (!is.logical(lncHighlight) || length(lncHighlight) != 1) {
    stop("'lncHighlight' must be a single logical value (TRUE/FALSE).")
  }
  
  # Extract gene ID for the specified lncRNA
  lnc_id <- lacenObject$summdf$gene_id[lacenObject$summdf$gene_name == lncName & lacenObject$summdf$is_nc][1]
  
  if (is.na(lnc_id)) {
    stop(paste("lncRNA '", lncName, "' is not classified as non-coding.", sep = ""))
  }
  
  # Extract TOM values for the lncRNA
  tom_values <- lacenObject$TOM[rownames(lacenObject$TOM) == lnc_id, ]
  
  # Identify the module of the lncRNA
  module <- lacenObject$summdf$module[lacenObject$summdf$gene_id == lnc_id]
  
  # Extract TOM values within the module
  tom_module <- tom_values[ names(tom_values) %in% lacenObject$summdf$gene_id[lacenObject$summdf$module == module] ]
  
  # Calculate median connectivity within the module
  median_connectivity <- stats::median(tom_module, na.rm = TRUE)
  
  # Select genes with connectivity above the median
  tom_selected <- tom_module[tom_module > median_connectivity]
  
  # Adjust nGenes based on availability
  if (nGenes == "OM") {
    nGenes <- length(tom_selected)
  } else if (length(tom_selected) < nGenes) {
    message(paste("Only ", length(tom_selected), " genes available with connectivity above the median. Setting nGenes to this value.", sep = ""))
    nGenes <- length(tom_selected)
  }
  
  # Select top nGenes based on connectivity
  tom_selected <- sort(tom_selected, decreasing = TRUE)[1:nGenes]
  
  # Extract top correlated genes
  genes_to_enrich <- names(tom_selected)
  
  # Perform enrichment analysis using gprofiler2
  enrichment <- tryCatch({
    gprofiler2::gost(
      query = genes_to_enrich,
      evcodes = TRUE,
      multi_query = FALSE,
      ordered_query = TRUE,
      user_threshold = 0.05,
      custom_bg = lacenObject$summdf$gene_id,
      sources = sources,
      organism = organism
    )
  }, error = function(e) {
    stop("gprofiler2::gost failed: ", e$message)
  })
  
  if (is.null(enrichment)) {
    stop("No enrichment terms found.")
  }
  
  # Select top enriched terms to highlight
  if (nHighlight > nrow(enrichment$result)) {
    nHighlight <- nrow(enrichment$result)
    warning("nHighlight exceeds the number of enrichment terms. Adjusted to ", nHighlight)
  }
  highlight_terms <- enrichment$result$term_id[1:nHighlight]
  
  # Save the enrichment plot
  enr_path <- ifelse(is.null(list(...)[["enrPath"]]), 
                     paste0("./", lncName, "_enr.png"), 
                     list(...)[["enrPath"]])
  
  grDevices::png(filename = enr_path, width = 1200, height = 1000)
  gprofiler2::publish_gostplot(
    gprofiler2::gostplot(enrichment, interactive = FALSE),
    highlight_terms = highlight_terms,
    filename = enr_path,
    width = 12,
    height = 10
  )
  grDevices::dev.off()
  
  # Prepare connectivity dataframe for export
  connectivity_df <- data.frame(
    gene_id = names(tom_selected), 
    connectivity = tom_selected,
    gene_name = lacenObject$annotationData$gene_name[match(names(tom_selected), lacenObject$annotationData$gene_id)],
    is_nc = lacenObject$summdf$is_nc[match(names(tom_selected), lacenObject$summdf$gene_id)],
    module = lacenObject$summdf$module[match(names(tom_selected), lacenObject$summdf$gene_id)],
    kTotal = lacenObject$summdf$kTotal[match(names(tom_selected), lacenObject$summdf$gene_id)],
    kWithin = lacenObject$summdf$kWithin[match(names(tom_selected), lacenObject$summdf$gene_id)],
    log2FC = lacenObject$summdf$log2FC[match(names(tom_selected), lacenObject$summdf$gene_id)],
    pval = lacenObject$summdf$pval[match(names(tom_selected), lacenObject$summdf$gene_id)],
    is_deg = abs(lacenObject$summdf$log2FC[match(names(tom_selected), lacenObject$summdf$gene_id)]) >= 1 & 
              lacenObject$summdf$pval[match(names(tom_selected), lacenObject$summdf$gene_id)] <= 0.05,
    stringsAsFactors = FALSE
  )
  
  # Save connectivity data as CSV
  connec_path <- ifelse(is.null(list(...)[["connecPath"]]), "./connectivities.csv", list(...)[["connecPath"]])
  utils::write.csv(connectivity_df, file = connec_path, row.names = FALSE)
  
  # Save enrichment results as CSV
  enr_csv_path <- ifelse(is.null(list(...)[["enrCsvPath"]]), "./enrichment.csv", list(...)[["enrCsvPath"]])
  utils::write.csv(enrichment$result, file = enr_csv_path, row.names = FALSE)
  
  # Prepare the TOM for network visualization
  tom_subset <- tom_matrix[rownames(tom_matrix) %in% genes_to_enrich,
                           colnames(tom_matrix) %in% genes_to_enrich]
  tom_subset[tom_subset < stats::median(tom_matrix, na.rm = TRUE)] <- 0
  diag(tom_subset) <- 0  # Remove self-connections
  
  # Create igraph object from TOM
  graph_net <- igraph::graph_from_adjacency_matrix(tom_subset, mode = "undirected", weighted = TRUE)
  
  # Limit the number of terms if it exceeds the maximum allowed
  if (nTerm > 32) {
    warning("nTerm exceeds the maximum allowed (32). Setting nTerm to 32.")
    nTerm <- 32
  }
  
  if (!is.null(enrichment)) {
    # Select top nTerm enrichment terms
    top_terms <- head(enrichment$result, n = nTerm)
    
    # Assign module attributes to nodes based on enrichment
    for (i in 1:nrow(top_terms)) {
      term <- top_terms$term_id[i]
      genes_in_term <- unlist(strsplit(top_terms$intersection[i], ","))
      V(graph_net)$color[V(graph_net)$name %in% genes_in_term] <- "red"
    }
    
    # Set layout for the network
    layout_coords <- ggraph::create_layout(graph_net, layout = "fr")
    
    # Generate network plot using ggraph
    network_plot <- ggraph::ggraph(layout_coords) +
      ggraph::geom_edge_link(aes(width = weight), alpha = 0.5) +
      ggraph::geom_node_point(aes(color = color), size = 5) +
      ggraph::geom_node_text(aes(label = name, color = color), repel = TRUE) +
      ggplot2::scale_color_manual(values = c("red" = "red", "black" = "black")) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(paste("Network for lncRNA:", lncName)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Save network plot
    net_path <- ifelse(is.null(list(...)[["netPath"]]), 
                       paste0("./", lncName, "_net.png"), 
                       list(...)[["netPath"]])
    grDevices::png(filename = net_path, width = 1800, height = 1200, res = 150)
    print(network_plot)
    grDevices::dev.off()
  }
}