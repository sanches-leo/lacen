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
    message(paste("Only ", length(tom_selected), " genes available with connectivity above the median. Setting nGenes and nGenesNet to this value.", sep = ""))
    nGenes <- length(tom_selected) }
  if (nGenesNet > nGenes){
    nGenesNet <- nGenes
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
  tom_subset <- tom_matrix[rownames(tom_matrix) %in% genes_to_enrich[1:nGenes],
                           colnames(tom_matrix) %in% genes_to_enrich[1:nGenes]]
  tom_subset[tom_subset < stats::median(tom_matrix, na.rm = TRUE)] <- 0
  diag(tom_subset) <- 0  # Remove self-connections

  # Create igraph object from TOM
  graph_net <- igraph::graph_from_adjacency_matrix(tom_subset, mode = "undirected", weighted = TRUE)

  # Limit the number of terms if it exceeds the maximum allowed
  if (nTerm > 20) {
    warning("nTerm exceeds the maximum allowed (20). Setting nTerm to 20.")
    nTerm <- 20
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

    # if(base::startsWith(x = sources, prefix = "GO") & organism == "hsapiens"){
    #   ont <- unlist(base::strsplit(sources, split = ":"))[2]
    #   org_db <- "org.Hs.eg.db"
    #   sim_matrix <- rrvgo::calculateSimMatrix(top_terms$term_id,
    #                                           orgdb = org_db,
    #                                           semdata = GOSemSim::godata(annoDb = org_db, ont = ont, keytype = "ENTREZID"),
    #                                           ont = ont,
    #                                           method = "Rel")
    #
    #   # Assign scores based on p-values
    #   scores <- stats::setNames(-log10(top_terms$p_value), top_terms$term_id)
    #
    #   # Reduce the similarity matrix
    #   reduced_terms <- rrvgo::reduceSimMatrix(sim_matrix,
    #                                           scores,
    #                                           threshold = 0.7,
    #                                           orgdb = org_db)
    #
    # }


    for(i in 1:nrow(top_terms)){
      term <- top_terms[i, "term_name"]
      term <- substr(term, 1, 40)
      enr_int <- unlist(strsplit(top_terms[i, "intersection"], ","))
      intersection <- as.numeric(rownames(tom_subset) %in% enr_int)
      term <- paste(term, " (", length(enr_int),"/", length(genes_to_enrich),")", sep = "")
      graph_net <- igraph::set_vertex_attr(graph_net, term, index = igraph::V(graph_net), intersection)
      enrichment_categories <- c(enrichment_categories, term)
      no_term <- rbind(no_term, intersection)
    }

    # # Set layout for the network
    # layout_coords <- graphlayouts::layout_with_stress(graph_net)
    # igraph::V(graph_net)$x <- layout_coords[, 1]
    # igraph::V(graph_net)$y <- layout_coords[, 2]

    # Identify index of the node of interest
    center_name <- lnc_id
    center_index <- which(igraph::V(graph_net)$name == center_name)

    # Reorder the graph so the center node is first
    graph_net <- igraph::permute(graph_net, c(center_index, setdiff(seq_len(igraph::vcount(graph_net)), center_index)))

    # Create layout: center node at origin, others in a circle
    theta <- seq(0, 2 * pi, length.out = igraph::vcount(graph_net))
    layout_coords <- matrix(0, nrow = igraph::vcount(graph_net), ncol = 2)
    layout_coords[1, ] <- c(0, 0)  # center node at (0, 0)
    layout_coords[-1, ] <- cbind(cos(theta[-1]), sin(theta[-1]))

    # Assign coordinates to vertex attributes
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
    # is_enriched_node[which((lncName == igraph::V(graph_net)$label))] <- 0
    is_nc <- connectivity_df$gene_name[connectivity_df$is_nc]
    is_enriched_node[igraph::V(graph_net)$label %in% is_nc] <- 0
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
    path_pal <- c("#22FD00",
                  "#F6DC9A",
                  "#FE0DD0",
                  "#00FCE8",
                  "#AE0DFF",
                  "#A44563",
                  "#006216",
                  "#E4DDF4",
                  "#2690FF",
                  "#424080",
                  "#F099FE",
                  "#C8F01C",
                  "#FD9F2A",
                  "#F41688",
                  "#7fd1ee",
                  "#6EEC90",
                  "#635C5A",
                  "#8E5600",
                  "#4d915b",
                  "#DB3B0D")
    path_pal <- path_pal[seq(1,(length(enrichment_categories)-1))]

    # Setting the edge colors of lnc to black
    # edge_list <- igraph::get.edgelist(graph_net)
    edge_list <- igraph::as_edgelist(graph_net)
    edge_color <- ifelse(edge_list[,1] == lnc_id | edge_list[,2] == lnc_id,
                         "grey75", "grey30")
    igraph::E(graph_net)$edge_color <- edge_color

    # Set the edge alpha based on the nodes number
    nGenesToScale <- nGenesNet

    if(nGenesToScale < 20){nGenesToScale <- 20}

    alpha <- 40 * 0.35 / nGenesToScale
    scale_range <- c(0.5, 2.5) * 30 / nGenesToScale

    # Setting the edge alpha of lnc to 0.15
    edge_alpha <- ifelse(edge_list[,1] == lnc_id | edge_list[,2] == lnc_id,
                         0.5, alpha)
    igraph::E(graph_net)$edge_alpha <- edge_alpha

    # Setting the edge weights to less connected genes to 0
    edge_weight <- igraph::E(graph_net)$weight
    sel <- (edge_weight < quantile(edge_weight, 0.9)) & (edge_list[,1] != lnc_id) & (edge_list[,2] != lnc_id)
    edge_weight[sel] <- NA
    igraph::E(graph_net)$weight <- edge_weight

    # Setting the node size
    scalaling_factor <- 32/nGenesToScale

    node_size <- 24 * scalaling_factor

    lnc_node_scale <- 2 * scalaling_factor

    text_size <- 5 * scalaling_factor

    # Generate network plot using ggraph
    suppressWarnings({
      network_plot <- ggraph::ggraph(graph_net,
                                     "manual",
                                     x = igraph::V(graph_net)$x,
                                     y = igraph::V(graph_net)$y) +
        ggraph::geom_edge_link0(ggplot2::aes(width = weight),
                                color = igraph::E(graph_net)$edge_color,
                                alpha = igraph::E(graph_net)$edge_alpha,
                                show.legend = FALSE) +
        ggraph::scale_edge_width(range = scale_range)+
        scatterpie::geom_scatterpie(
          cols = enrichment_categories,
          data = igraph::as_data_frame(graph_net, "vertices"),
          colour = NA,
          pie_scale = lnc_node_scale,
          legend_name = "GO"
        ) +
        ggplot2::scale_fill_manual(values = c("#808080",
                                              path_pal)) +
        ggraph::geom_node_point(data = subgraph_nc,
                                ggplot2::aes(x = x, y = y),
                                shape = 23,
                                size = node_size,
                                show.legend = FALSE,
                                fill = "#808080",
                                colour = "#808080") +
        ggplot2::coord_fixed() +
        ggraph::geom_node_text(ggplot2::aes(label = label),
                               repel = FALSE,
                               colour = igraph::V(graph_net)$lcolor,
                               fontface = "bold",
                               size = text_size
        ) +
        ggraph::theme_graph(base_family="sans")+
        ggplot2::theme(legend.position = "right", text = ggplot2::element_text(size = 16)) +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

    })

    #Save network plot
    lncNameClean <- gsub("[^A-Za-z0-9]", "", lncName)
    net_path <- ifelse(is.null(list(...)[["netPath"]]),
                       paste0("./", lncName, "_net.png"),
                       list(...)[["netPath"]])
    grDevices::png(filename = net_path, width = 1800, height = 1200, res = 100)
    print(network_plot)
    invisible(grDevices::dev.off())
  }
}

