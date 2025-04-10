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
#' @import gprofiler2 igraph ggraph scatterpie grid
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

  if (!lncName %in% lacenObject$summdf$gene_name) {
    stop(paste("lncName '", lncName, "' not found.", sep = ""))
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
  subDF <- lacenObject$summdf[lacenObject$summdf$gene_name == lncName, ]

  lnc_id <- subDF$gene_id[subDF$is_nc]

  lnc_module <- subDF$module > 0

  lnc_id <- lnc_id[lnc_module]

  if(length(lnc_id) > 1){
    warning(paste("More than 1 ID is identified for this lncRNA. Using ", lnc_id[1]))
    lnc_id <- lnc_id[1]
  }

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

  # Arrange the pValue order of enrichment to crescent
  enrichment$result <- enrichment$result[order(enrichment$result$p_value, decreasing = FALSE), ]

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
  enr_df <- as.data.frame(apply(enrichment$result,2,as.character))

  utils::write.csv(enr_df, file = enr_csv_path, row.names = FALSE)

  # Prepare the TOM for network visualization
  tom_matrix <- lacenObject$TOM
  tom_subset <- tom_matrix[rownames(tom_matrix) %in% genes_to_enrich[1:nGenesNet],
                           colnames(tom_matrix) %in% genes_to_enrich[1:nGenesNet]]
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
    top_terms <- utils::head(enrichment$result, n = nTerm)

    # Assign module attributes to nodes based on enrichment

    # for (i in 1:nrow(top_terms)) {
    #   term <- top_terms$term_id[i]
    #   genes_in_term <- unlist(strsplit(top_terms$intersection[i], ","))
    #   igraph::V(graph_net)$color[igraph::V(graph_net)$name %in% genes_in_term] <- "red"
    # }

    enrichment_categories <- c()
    no_term <- matrix(nrow = 0, ncol = nGenesNet)
    for(i in 1:nrow(top_terms)){
      term <- top_terms[i, "term_name"]
      enr_int <- unlist(strsplit(top_terms[i, "intersection"], ","))
      intersection <- as.numeric(rownames(tom_subset) %in% enr_int)
      term <- paste(term, " (", length(enr_int),"/", length(genes_to_enrich),")", sep = "")
      graph_net <- igraph::set_vertex_attr(graph_net, term, index = igraph::V(graph_net), intersection)
      enrichment_categories <- c(enrichment_categories, term)
      no_term <- rbind(no_term, intersection)
    }

    # Set layout for the network
    layout_coords <- graphlayouts::layout_with_stress(graph_net)
    igraph::V(graph_net)$x <- layout_coords[, 1]
    igraph::V(graph_net)$y <- layout_coords[, 2]

    # Change the nodes names to gene name from geneID
    igraph::V(graph_net)$label <- lacenObject$annotationData$gene_name[match(igraph::V(graph_net)$name,
                                                                             lacenObject$annotationData$gene_id)]

    # Set the nodes words colors
    igraph::V(graph_net)$lcolor <- "black"
    if(isTRUE(lncHighlight)){
      igraph::V(graph_net)$lcolor[igraph::V(graph_net)$label %in% lacenObject$ncAnnotation$gene_id] <- "red"
    }
    # igraph::V(graph_net)$lcolor[igraph::V(graph_net)$label == lncName] <- "red4"

    # Assign lncRNA query to network (red)
    # is_query <- as.numeric(lncName == igraph::V(graph_net)$label)
    # graph_net <- igraph::set_vertex_attr(graph_net, "lncRNA query", index = igraph::V(graph_net), is_query)
    # enrichment_categories <- c("lncRNA query", enrichment_categories)

    # Assign 0 to the non enriched nodes
    is_enriched_node <- colSums(no_term)
    is_enriched_node <- ifelse(is_enriched_node == 0, 1, 0)
    is_enriched_node[which((lncName == igraph::V(graph_net)$label))] <- 0
    graph_net <- igraph::set_vertex_attr(graph_net, "no enriched term", index = igraph::V(graph_net), is_enriched_node)
    enrichment_categories <- c("no enriched term", enrichment_categories)

    # Create a dataframe with the coordinates of nc nodes
    igraph::V(graph_net)$is_nc <- as.numeric(connectivity_df$is_nc[match(igraph::V(graph_net)$label,
                                                                         connectivity_df$gene_name)])
    subgraph_nc <- data.frame(x = igraph::V(graph_net)$x,
                              y = igraph::V(graph_net)$y,
                              label = igraph::V(graph_net)$label,
                              is_nc = as.numeric(igraph::V(graph_net)$is_nc))
    subgraph_nc <- subgraph_nc[subgraph_nc$is_nc == 1, ]

    # Create and filter the color palette
    path_pal <- c('#38871C',
                  '#3D7CAC',
                  '#AF4BB4',
                  '#A1683B',
                  '#EE002E',
                  '#2A26FB',
                  '#D20D7F',
                  '#3B4738',
                  '#C74D00',
                  '#BE16F6',
                  '#620038',
                  '#0016AD',
                  '#886A90',
                  '#D80DB2',
                  '#9F3242',
                  '#7A7A2E',
                  '#1C1C62',
                  '#0D866C',
                  '#C60045',
                  '#B44288',
                  '#9856D7',
                  '#2675DE',
                  '#7960AD',
                  '#C90DD1',
                  '#94686A',
                  '#22606E',
                  '#661669',
                  '#5A2A00',
                  '#8422E7',
                  '#472E45',
                  '#AD5870',
                  '#3B6322',
                  '#B65C38',
                  '#CF2A0D',
                  '#A35697',
                  '#510080')
    path_pal <- path_pal[seq(1,(length(enrichment_categories)-1))]

    # Generate network plot using ggraph

    suppressWarnings({
      network_plot <- ggraph::ggraph(graph_net,
                                     "manual",
                                     x = igraph::V(graph_net)$x,
                                     y = igraph::V(graph_net)$y) +
        ggraph::geom_edge_link0(alpha = .25,
                                ggplot2::aes(width = weight),
                                show.legend = F) +
        ggraph::scale_edge_width(range = c(0.5, 2.5))+
        scatterpie::geom_scatterpie(
          cols = enrichment_categories,
          data = igraph::as_data_frame(graph_net, "vertices"),
          colour = NA,
          pie_scale = 2,
          legend_name = "GO"
        ) +
        ggplot2::scale_fill_manual(values = c("#808080",
                                              path_pal)) +
        ggraph::geom_node_point(data = subgraph_nc,
                                ggplot2::aes(x = x, y = y),
                                shape = 23,
                                size = 30,
                                show.legend = FALSE,
                                fill = "#808080",
                                colour = "#808080") +
        ggplot2::coord_fixed() +
        ggraph::geom_node_text(ggplot2::aes(label = label),
                               repel = FALSE,
                               colour = igraph::V(graph_net)$lcolor,
                               fontface = "bold",
                               size = 5
        ) +
        ggraph::theme_graph(base_family="sans")+
        ggplot2::theme(legend.position = "right", text = ggplot2::element_text(size = 24))
    })

    # Save network plot
    net_path <- ifelse(is.null(list(...)[["netPath"]]),
                       paste0("./", lncName, "_net.png"),
                       list(...)[["netPath"]])
    grDevices::png(filename = net_path, width = 1800, height = 1200, res = 100)
    print(network_plot)
    invisible(grDevices::dev.off())
  }
}

