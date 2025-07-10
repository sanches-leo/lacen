#' Heatmap of Top Connectivity
#'
#' Generates a heatmap of a specified module or submodule, showing interconnectedness based on the TOM.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param module Numeric value indicating the module to analyze.
#' @param submodule Optional. Numeric value indicating the submodule to analyze. Set to FALSE to ignore.
#' @param hmDimensions Numeric value indicating the number of genes to display in the heatmap.
#' @param filename Filename to save the heatmap as PNG.
#' @param removeNonDEG Logical. If TRUE, displays only differentially expressed genes.
#' @param outTSV Logical or character string. If FALSE, does not save. If a path is provided, saves the TSV with heatmap values.
#' @param ... Additional parameters for internal functions.
#'
#' @return None (plots and saves the heatmap).
#' @export
#' @import WGCNA igraph ggraph grid ggplot2
heatmapTopConnectivity <- function(lacenObject,
                                   module,
                                   submodule = FALSE,
                                   hmDimensions = FALSE,
                                   filename = FALSE,
                                   removeNonDEG = FALSE,
                                   outTSV = FALSE,
                                   ...) {
  UseMethod("heatmapTopConnectivity")
}

#' @export
heatmapTopConnectivity.lacen <- function(
    lacenObject,
    module,
    submodule = FALSE,
    hmDimensions = FALSE,
    filename = FALSE,
    removeNonDEG = FALSE,
    outTSV = FALSE,
    ...) {
  # Extract required data from the lacenObject
  datExpr <- lacenObject$datExpr           # Gene expression data
  traits <- lacenObject$datTraits$Trait      # Trait information (numeric vector)
  TOM <- lacenObject$TOM                     # Topological Overlap Matrix
  summdf <- lacenObject$summdf               # Summary data frame with module and connectivity info
  rrvgolist <- lacenObject$rrvgolist         # List of reduced enrichment results per module
  
  # Helper function to check if a parameter is invalid.
  # Returns TRUE if x is missing, NULL, or empty, or if it does not meet the expected class.
  invalidPar <- function(x, class) {
    if (missing(x) || is.null(x) || length(x) == 0) {
      return(TRUE)
    }
    else if (is.list(x) & class == "list") {
      return(FALSE)
    }
    else if (is.vector(x) & class == "cvector" & is.character(x)) {
      return(FALSE)
    }
    else if (is.character(x) & class == "character") {
      return(FALSE)
    }
    else if (is.numeric(as.numeric(x)) & class == "numeric") {
      return(FALSE)
    }
    else if (is.logical(x) & class == "logical") {
      return(FALSE)
    }
    else if (is.vector(x) & class == "nvector" & is.numeric(x)) {
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
  
  # Validate module parameter: must be a positive numeric value and exist in summdf
  if (invalidPar(module, "numeric") | module <= 0) {
    stop("invalid format. Module should be a positive integer")
  } else if (!module %in% names(table(summdf$module))) {
    possible <- base::sort(base::unique(summdf$module))
    possible <- possible[possible > 0]
    possible <- base::paste(possible, collapse = ", ")
    stop("Invalid module. Possible modules are: ", possible)
  }
  
  # Validate submodule if provided (must be a positive integer and exist in the reduced enrichment list)
  if (!isFALSE(submodule)) {
    if (invalidPar(submodule, "numeric") | submodule <= 0) {
      stop("invalid format. Module should be a positive integer")
    }
    else if (!as.character(submodule) %in% as.character(unique(rrvgolist[[as.character(module)]]$cluster))) {
      possible <- base::sort(unique(rrvgolist[[as.character(module)]]$cluster))
      if(length(possible) == 0){
        stop("This module couldn't be enriched, thus there are no submodules.")
      } else {
        possible <- base::paste(possible, collapse = ", ")
        stop("Invalid submodule. Possible submodules are: ", possible)
      }
    }
  }
  
  # Validate filename: must be FALSE or a character string
  if (!isFALSE(filename)) {
    if (invalidPar(filename, "character")) {
      stop("filename should be FALSE or a string")
    }
  }
  
  # Validate removeNonDEG: must be a logical value
  if (invalidPar(removeNonDEG, "logical")) {
    stop("removeNonDEG should be a logical value")
  }
  
  # Validate outTSV: must be logical or a character string (file path)
  if (invalidPar(outTSV, "logical") & !is.character(outTSV)) {
    stop("outTSV should be a logical value or a path and filename")
  }
  
  # Validate traits: must be a numeric vector
  if (invalidPar(traits, "nvector")) {
    stop("traits should be a numeric vector")
  }
  
  # Helper function: dataToColor
  # Returns a vector of colors for each gene based on either differential expression (if logDEG==TRUE)
  # or connectivity values (if logDEG==FALSE). Colors differ based on thresholds.
  dataToColor <- function(x, summdf, module, logDEG = FALSE, mod) {
    if (logDEG) {
      # Extract DEG data for genes in x and ensure the order matches
      degdf <- summdf[summdf$gene_id %in% x, c("gene_id", "log2FC", "pval")]
      degdf <- degdf[match(x, degdf$gene_id), ]
      colors <- rep("black", length(x))
      # Assign red shades for upregulated genes and blue shades for downregulated ones
      colors[which(degdf$log2FC > 1 & degdf$pval <= 0.05)] <- "#FF0000"
      colors[which((degdf$log2FC > 1 & degdf$pval > 0.05) | (degdf$log2FC > 0 & degdf$log2FC < 1))] <- "#FF9696"
      colors[which(degdf$log2FC < -1 & degdf$pval <= 0.05)] <- "#0000FF"
      colors[which((degdf$log2FC < -1 & degdf$pval > 0.05) | (degdf$log2FC < 0 & degdf$log2FC > -1))] <- "#9696FF"
    }
    else {
      # For connectivity values, limit x to the 10th and 90th percentiles
      stats <- summdf[summdf$module == mod, ]$kWithin
      maxvalue <- stats::quantile(stats, probs = 0.9)
      minvalue <- stats::quantile(stats, probs = 0.1)
      x[x < minvalue] <- minvalue
      x[x > maxvalue] <- maxvalue
      # Create a palette from light green to dark green and assign colors based on intervals
      pal <- grDevices::colorRampPalette(c("#c8fac8", "#003200"))(128)
      colors <- pal[findInterval(x, seq(minvalue, maxvalue, length.out = length(pal)), all.inside = FALSE)]
    }
    return(colors)
  }
  
  # Helper function: termsToOut
  # Merges reduced enrichment terms into the output table for each gene.
  termsToOut <- function(out, rrvgolist, module) {
    module <- as.character(module)
    rrvgoterms <- rrvgolist[[module]]
    rrvgoterms <- rrvgoterms[, c("cluster", "parentTerm2")]
    rrvgoterms <- rrvgoterms[!base::duplicated(temp), ]
    rrvgoterms <- rrvgoterms[order(rrvgoterms$cluster), ]
    rownames(rrvgoterms) <- NULL
    termsCol <- NULL
    # For each row in the output table, concatenate the term(s) from the reduced enrichment results
    for (l in 1:nrow(out)) {
      sel <- unlist(out[l, which(colnames(out) == "1"):dim(out)[2]])
      if (sum(sel) > 0) {
        pastedTerm <- paste(rrvgoterms$parentTerm2[sel], collapse = ", ")
      }
      else {
        pastedTerm <- " "
      }
      termsCol <- c(termsCol, pastedTerm)
    }
    out <- out[, c(1:4, 7, 11, 15)]
    out$submodules <- termsCol
    return(out)
  }
  
  # Define a custom heatmap function (heatmap.3)
  # This function extends standard heatmap plotting to include additional annotations.
  # (For brevity, internal details follow standard procedures to generate heatmaps with dendrograms,
  # color keys, and side color bars.)
  heatmap.3 <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = stats::dist, hclustfun = stats::hclust,
                        dendrogram = c("both", "row", "column", "none"), symm = FALSE,
                        scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv, "Rowv"),
                        add.expr, breaks, symbreaks = max(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors",
                        colsep, rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
                        notecol = "cyan", na.color = graphics::par("bg"), trace = c("none", "column", "row", "both"),
                        tracecol = "cyan", hline = stats::median(breaks), vline = stats::median(breaks),
                        linecol = tracecol, margins = c(0, 0), ColSideColors, RowSideColors,
                        side.height.fraction = 0.3, cexRow = mcr, cexCol = mcc, labRow = NULL, labCol = NULL,
                        key = FALSE, keysize = 1.5, density.info = c("none", "histogram", "density"),
                        denscol = tracecol, symkey = max(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
                        main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
                        ColSideColorsSize = 1, RowSideColorsSize = 1, KeyValueName = "Value", keyvalues = TRUE, ...) {
    # (Extensive code to compute dendrograms, layout, scaling, and plotting is implemented here.)
    # For clarity, the internal details are not modified but follow standard procedures.
    # The function returns a list with information about the heatmap.
    invalid <- function(x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      (x - low) / (high - low)
    }
    retval <- list()
    scale <- if (symm && missing(scale)) "none" else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale when breaks are specified can produce unpredictable results.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    cr <- 0.2 + 1/log10(nr)
    mcr <- ifelse(cr > 1.3, 1.1, cr)
    cc <- 0.2 + 1/log10(nc)
    mcc <- ifelse(cc > 1.3, 1.3, cc)
    if (nr <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dendrogram <- "none"
        warning("Rowv is FALSE; omitting row dendrogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Colv is FALSE; omitting column dendrogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- stats::order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- stats::as.dendrogram(hcr)
      ddr <- stats::reorder(ddr, Rowv)
      rowInd <- stats::order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("Row dendrogram ordering gave wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- stats::as.dendrogram(hcr)
      ddr <- stats::reorder(ddr, Rowv)
      rowInd <- stats::order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("Row dendrogram ordering gave wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- stats::order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- stats::order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm) x else t(x)))
      ddc <- stats::as.dendrogram(hcc)
      ddc <- stats::reorder(ddc, Colv)
      colInd <- stats::order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("Column dendrogram ordering gave wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm) x else t(x)))
      ddc <- stats::as.dendrogram(hcc)
      ddc <- stats::reorder(ddc, Colv)
      colInd <- stats::order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("Column dendrogram ordering gave wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, stats::sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, stats::sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (inherits(col, "function"))
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)) {
        graphics::par(mar = c(margins[1], 0, 0, 0))
        graphics::image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      }
      else {
        graphics::par(mar = c(margins[1], 0, 0, 0))
        rsc <- t(RowSideColors[, rowInd, drop = F])
        rsc.colors <- matrix()
        rsc.names <- names(table(rsc))
        rsc.i <- 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] <- rsc.name
          rsc[rsc == rsc.name] <- rsc.i
          rsc.i <- rsc.i + 1
        }
        rsc <- matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        graphics::image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          graphics::axis(1, 0:(dim(rsc)[2] - 1)/max(1, (dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    if (!missing(ColSideColors)) {
      if (!is.matrix(ColSideColors)) {
        graphics::par(mar = c(0.5, 0, 0, margins[2]))
        graphics::image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      }
      else {
        graphics::par(mar = c(0.5, 0, 0, margins[2]))
        csc <- ColSideColors[colInd, , drop = F]
        csc.colors <- matrix()
        csc.names <- names(table(csc))
        csc.i <- 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] <- csc.name
          csc[csc == csc.name] <- csc.i
          csc.i <- csc.i + 1
        }
        csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
        graphics::image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          graphics::axis(4, 0:(dim(csc)[2] - 1)/max(1, (dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    graphics::par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    graphics::image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr),
                    axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
      mmat <- ifelse(is.na(x), 1, NA)
      graphics::image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color, add = TRUE)
    }
    graphics::axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    if (!is.null(xlab))
      graphics::mtext(xlab, side = 1, line = margins[1] - 1.25)
    graphics::axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab))
      graphics::mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) graphics::rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)),
                                          xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep),
                                          lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) graphics::rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5,
                                          xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
                                          lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          graphics::abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          graphics::abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      graphics::text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)
    graphics::par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
    }
    if (keyvalues) {
      {
        pal1 <- grDevices::colorRampPalette(c("#c8fac8", "#003200"))(128)
        pal2 <- rev(grDevices::grey.colors(128))
        plot(0, 0, xlim = c(0, 2), ylim = c(0, 128), col = "white", axes = FALSE, ylab = "", xlab = "")
        graphics::axis(2, at = 1, labels = "Low", tick = 0, las = 2)
        graphics::axis(2, at = 127, labels = "High", tick = 0, las = 2)
        graphics::axis(3, at = 0.5, labels = "Connectivity", tick = 0, las = 2)
        graphics::axis(3, at = 1.5, labels = "Topological Overlap", tick = 0, las = 2)
        for (i in 0:127) {
          graphics::polygon(c(0:1, 1:0, 0) + 0, c(0, 0:1, 1:0) + i, col = pal1[i], border = FALSE)
          graphics::polygon(c(0:1, 1:0, 0) + 1, c(0, 0:1, 1:0) + i, col = pal2[i], border = FALSE)
        }
      }
      if (isFALSE(removeNonDEG)) {
        plot(0, 0, xlim = c(0, 2), ylim = c(0, 3), col = "white", axes = FALSE, ylab = "", xlab = "")
        graphics::axis(2, at = 2.5, labels = "upregulated", tick = 0, las = 2)
        graphics::axis(2, at = 1.5, labels = "downregulated", tick = 0, las = 2)
        graphics::axis(3, at = 0.5, labels = "significant", tick = 0, las = 2)
        graphics::axis(3, at = 1.5, labels = "non-significant", tick = 0, las = 2)
        graphics::axis(2, at = 0.5, labels = "no DEG data", tick = 0, las = 2)
        graphics::polygon(c(0, 1, 1, 0) + 0, c(0, 0, 1, 1) + 1, col = "#0000FF", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 1, c(0, 0, 1, 1) + 1, col = "#9696FF", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 0, c(0, 0, 1, 1) + 2, col = "#FF0000", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 1, c(0, 0, 1, 1) + 2, col = "#FF9696", border = FALSE)
        graphics::polygon(c(0.5, 1.5, 1.5, 0.5) + 0, c(0, 0, 1, 1) + 0, col = "black", border = FALSE)
      }
      else {
        plot(0, 0, xlim = c(0, 2), ylim = c(0, 2), col = "white", axes = FALSE, ylab = "", xlab = "")
        graphics::axis(2, at = 1.5, labels = "upregulated", tick = 0, las = 2)
        graphics::axis(2, at = 0.5, labels = "downregulated", tick = 0, las = 2)
        graphics::axis(3, at = 0.5, labels = "significant", tick = 0, las = 2)
        graphics::axis(3, at = 1.5, labels = "non-significant", tick = 0, las = 2)
        graphics::polygon(c(0, 1, 1, 0) + 0, c(0, 0, 1, 1) + 0, col = "#0000FF", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 1, c(0, 0, 1, 1) + 0, col = "#9696FF", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 0, c(0, 0, 1, 1) + 1, col = "#FF0000", border = FALSE)
        graphics::polygon(c(0, 1, 1, 0) + 1, c(0, 0, 1, 1) + 1, col = "#FF9696", border = FALSE)
      }
    }
    else graphics::plot.new()
    graphics::par(mar = c(margins[1], 0, 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else graphics::plot.new()
    if (!is.null(main))
      graphics::mtext(main, cex = 3, adj = 0.6, outer = TRUE)
    if (key) {
      graphics::par(mar = c(0, 0, 0, 0), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      z <- seq(min.raw, max.raw, length = length(col))
      graphics::image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, xaxt = "n", yaxt = "n")
      graphics::par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      graphics::axis(1, at = xv, labels = lv)
      if (scale == "row")
        graphics::mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        graphics::mtext(side = 1, "Column Z-Score", line = 2)
      else graphics::mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- stats::density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        graphics::lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, lwd = 1)
        graphics::axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        graphics::title("Color Key\nand Density Plot")
        graphics::mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- graphics::hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        graphics::lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", col = denscol)
        graphics::axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        graphics::title("Color Key\nand Histogram")
        graphics::mtext(side = 2, "Count", line = 2)
      }
      else graphics::title("Color Key")
    }
    else graphics::plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  
  # If filename is FALSE, set a default filename based on module number (only if submodule is FALSE)
  if (isFALSE(filename)) {
    if (isFALSE(submodule)) {
      filename <- paste("7_heatmap_", module, ".png", sep = "")
    }
  }
  
  # Determine the set of genes to plot:
  # If no submodule is specified, take protein-coding genes in the module.
  # If submodule is specified, filter based on the submodule flag in summdf.
  if (isFALSE(submodule)) {
    genes_in_mod_and_subm <- summdf$module == module & !summdf$is_nc
  } else {
    genes_in_mod_and_subm <- summdf$module == module & summdf[, as.character(submodule)] & !summdf$is_nc
  }
  genes_in_mod_and_submod <- summdf$gene_id[genes_in_mod_and_subm]
  
  # For the given module, calculate the median connectivity (kWithin)
  degrees_mod <- summdf[summdf$module == module, ]
  median_mod <- stats::median(degrees_mod$kWithin)
  # Identify lncRNAs with connectivity above the module median
  lnc_over_median <- degrees_mod$gene_id[degrees_mod$kWithin > median_mod & degrees_mod$is_nc]
  
  # Subset the Topological Overlap Matrix (TOM) for lncRNAs (rows) and protein-coding genes (columns)
  adj_mod_cluster <- TOM[rownames(TOM) %in% lnc_over_median,
                         colnames(TOM) %in% genes_in_mod_and_submod]
  
  # If only one lncRNA is above the median, print info and halt execution.
  if(is.null(dim(lnc_over_median))){
    if(length(lnc_over_median) == 1){
      print(paste("Only the lncRNA ", lnc_over_median, " has connectivity over the module median", sep = ""))
      print("Top 10 connected pc genes: ")
      print(sort(adj_mod_cluster, decreasing = TRUE)[1:10])
      stop("Execution halted")
    }
  } else if (!(nrow(adj_mod_cluster) > 0 & ncol(adj_mod_cluster) > 0)) {
    stop("No lncrnas with connectivity over the module median")
  }
  
  # Compute module eigengenes for the selected protein-coding genes and order them
  MEs0 <- WGCNA::moduleEigengenes(datExpr, as.numeric(genes_in_mod_and_subm))$eigengenes
  MEs <- WGCNA::orderMEs(MEs0)
  
  # Calculate correlation between module eigengenes and traits (only the first value is used)
  submodTraitCor <- unlist(WGCNA::cor(MEs, traits, use = "p"))[1]
  submoduleTraitPvalue <- unlist(WGCNA::corPvalueStudent(submodTraitCor, nrow(datExpr)))[1]
  
  # Retrieve connectivity and DEG (differential expression) data for protein-coding genes
  pc_connect <- summdf[summdf$gene_id %in% colnames(adj_mod_cluster),
                       c("kWithin", "gene_id")]
  pc_connect[is.na(pc_connect)] <- 0
  pc_connect <- stats::setNames(pc_connect$kWithin, pc_connect$gene_id)
  pc_connect <- sort(pc_connect, decreasing = TRUE)
  
  pc_DEG <- summdf[summdf$gene_id %in% colnames(adj_mod_cluster),
                   c("log2FC", "gene_id")]
  pc_DEG[is.na(pc_DEG)] <- 0
  pc_DEG <- pc_DEG$gene_id
  pc_DEG <- pc_DEG[match(names(pc_connect), pc_DEG)]
  
  # Retrieve connectivity and DEG data for lncRNAs
  lnc_connect <- summdf[summdf$gene_id %in% rownames(adj_mod_cluster),
                        c("kWithin", "gene_id")]
  lnc_connect[is.na(lnc_connect)] <- 0
  lnc_connect <- stats::setNames(lnc_connect$kWithin, lnc_connect$gene_id)
  lnc_connect <- sort(lnc_connect, decreasing = TRUE)
  
  lnc_DEG <- summdf[summdf$gene_id %in% rownames(adj_mod_cluster),
                    c("log2FC", "gene_id")]
  lnc_DEG[is.na(lnc_DEG)] <- 0
  lnc_DEG <- lnc_DEG$gene_id
  lnc_DEG <- lnc_DEG[match(names(lnc_connect), lnc_DEG)]
  
  # Reorder the TOM matrix to match the sorted order of lncRNA and protein-coding connectivity
  adj_mod_cluster <- adj_mod_cluster[match(names(lnc_connect), rownames(adj_mod_cluster)),
                                     match(names(pc_connect), colnames(adj_mod_cluster))]
  
  # Assign colors for protein-coding connectivity (and DEG) using the dataToColor helper function
  pc_connect_color <- dataToColor(unlist(pc_connect), summdf, logDEG = FALSE, mod = module)
  pc_DEG_color <- dataToColor(unlist(pc_DEG), summdf, logDEG = TRUE, mod = module)
  lnc_connect_color <- dataToColor(unlist(lnc_connect), summdf, logDEG = FALSE, mod = module)
  lnc_DEG_color <- dataToColor(unlist(lnc_DEG), summdf, logDEG = TRUE, mod = module)
  
  # Create matrices for side colors (used for annotation in the heatmap)
  colside <- as.matrix(cbind(lnc_DEG_color, lnc_connect_color))
  rowside <- as.matrix(cbind(pc_connect_color, pc_DEG_color))
  
  # If hmDimensions is provided (not FALSE), adjust the heatmap dimensions by subsetting the matrix
  if (!isFALSE(hmDimensions)) {
    hmDimensions <- suppressWarnings(as.integer(hmDimensions))
    if (hmDimensions <= 0 | is.na(hmDimensions)) {
      stop("hmDimensions should be an integer bigger than 0")
    }
    else if (invalidPar(hmDimensions, "numeric")) {
      stop("hmDimensions should be FALSE or an integer bigger than 0")
    }
    else {
      hmDimensions1 <- hmDimensions2 <- hmDimensions
      if (hmDimensions > dim(adj_mod_cluster)[1]) {
        warning("hmDimensions exceeds the lnc gene number")
        hmDimensions1 <- dim(adj_mod_cluster)[1]
      }
      else if (hmDimensions > dim(adj_mod_cluster)[2]) {
        warning("hmDimensions exceeds the protein coding gene number")
        hmDimensions2 <- dim(adj_mod_cluster)[2]
      }
      adj_mod_cluster <- adj_mod_cluster[1:hmDimensions1, 1:hmDimensions2]
      colside <- colside[1:hmDimensions1, ]
      rowside <- rowside[1:hmDimensions2, ]
    }
  }
  
  # If removeNonDEG is TRUE, remove genes with no DEG data (indicated by black color)
  if (isTRUE(removeNonDEG)) {
    adj_mod_cluster <- adj_mod_cluster[!(colside[, 1] == "black"), !(rowside[, 2] == "black")]
    colside <- colside[!(colside[, 1] == "black"), ]
    rowside <- rowside[!(rowside[, 2] == "black"), ]
  }
  
  # Replace the row and column names of the TOM matrix with gene names from summdf
  rnames <- rownames(adj_mod_cluster)
  temp <- summdf[summdf$gene_id %in% rnames, c("gene_id", "gene_name")]
  temp <- temp[match(rnames, temp$gene_id), "gene_name"]
  temp[is.na(temp)] <- rnames[is.na(temp)]
  rownames(adj_mod_cluster) <- temp
  
  cnames <- colnames(adj_mod_cluster)
  temp <- summdf[summdf$gene_id %in% cnames, c("gene_id", "gene_name")]
  temp <- temp[match(cnames, temp$gene_id), "gene_name"]
  temp[is.na(temp)] <- cnames[is.na(temp)]
  colnames(adj_mod_cluster) <- temp
  
  # Label the side color matrices
  colnames(colside) <- c("Non-coding exp.", "Non-coding conn.")
  colnames(rowside) <- c("Protein-coding conn.", "Protein-coding exp.")
  
  # Prepare title lines for the heatmap
  hmline1 <- paste("Module:", as.character(module))
  hmline2 <- paste("\nSubmodule:", rrvgolist[[as.character(module)]][rrvgolist[[as.character(module)]]$cluster == submodule, "parentTerm2"][1])
  hmline3 <- ""  # Additional info (e.g., correlation, p-value) can be added here
  
  # Capture any additional arguments passed via ... and assign them locally
  args <- as.list(match.call(expand.dots = FALSE)$`...`)
  if(length(args) > 0){
    for(i in 1:length(args)){
      assign(names(args)[i], args[[i]])
    }
  }
  
  # Set default flag for plotting the heatmap if not provided
  if(!exists("plothm")){
    plothm <- TRUE
  }
  
  # Plot the heatmap if plothm is TRUE
  if(isTRUE(plothm)){
    if (isFALSE(submodule)) {
      titlehm <- paste(hmline1, hmline3)
    }
    else {
      titlehm <- paste(hmline1, hmline2, hmline3)
    }
    # Open a PNG device to save the heatmap image
    grDevices::png(filename = filename, width = 1800, height = 1800, units = "px")
    # Define a layout matrix for the heatmap and side color bars
    m <- cbind(rep(0, 22), c(c(rep(0, 6), rep(4, 6), rep(0, 5), 5, seq(7, 10))), rep(0, 22),
               c(0, 0, rep(1, 20)), c(6, 2, rep(3, 20)))
    graphics::layout(m, widths = c(0.2, 0.5, 0.2, 0.5, 10),
                     heights = c(0.5, 0.5, rep(0.5, 20)), respect = FALSE)
    graphics::par(oma = c(10, 8, 7, 10), cex = 1.75)
    # Generate the heatmap using the custom heatmap.3 function
    heatmap.3(t(adj_mod_cluster), col = rev(grDevices::grey.colors(100)),
              main = paste(titlehm), ColSideColors = colside, RowSideColors = t(rowside),
              dendrogram = "none", key = FALSE, Rowv = FALSE, Colv = FALSE,
              keyvalues = TRUE)
    grDevices::dev.off()
  }
  
  # If outTSV is specified, generate an output table with DEG data and enrichment terms.
  if (!isFALSE(outTSV)) {
    if(module %in% names(rrvgolist)){
      out <- summdf[summdf$gene_id %in% c(lnc_DEG, pc_DEG), ]
      out <- out[order(out$is_nc, out$kWithin, decreasing = TRUE), ]
      if (isTRUE(removeNonDEG)) {
        out <- out[!is.na(out$pval), ]
      }
      out <- termsToOut(out, rrvgolist, module)
      tsvname <- outTSV
      utils::write.table(out, file = tsvname, sep = "\t")
    }
  }
}
