#' Check Data
#'
#' Validates the format of data within a lacen object.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param orgdb A Bioconductor `OrgDb` object corresponding to the organism of interest (e.g., `org.Hs.eg.db` for human).
#'
#' @return Prints a confirmation message if data format is correct; otherwise, issues warnings.
#' @export
checkData <- function(lacenObject,
                      orgdb = "org.Hs.eg.db") {
  UseMethod("checkData")
}

#' @export
checkData.lacen <- function(lacenObject,
                            orgdb = "org.Hs.eg.db") {
  is_valid <- TRUE

  # Helper function to check column names
  check_columns <- function(df, expected_cols, df_name) {
    if (!all(colnames(df) == expected_cols)) {
      warning(paste0(df_name, ' should have columns: "', paste(expected_cols, collapse = '", "'), '".'))
      return(FALSE)
    }
    return(TRUE)
  }

  # Validate annotationData
  if (!is.null(lacenObject$annotationData)) {
    is_valid <- check_columns(lacenObject$annotationData, c("gene_id", "gene_name"), "annotationData") && is_valid
  } else {
    warning("annotationData is NULL.")
    is_valid <- FALSE
  }

  # Validate ncAnnotation
  if (!is.null(lacenObject$ncAnnotation)) {
    is_valid <- check_columns(lacenObject$ncAnnotation, c("gene_id", "gene_name"), "ncAnnotation") && is_valid
  } else {
    warning("ncAnnotation is NULL.")
    is_valid <- FALSE
  }

  # Validate datTraits
  if (!is.null(lacenObject$datTraits)) {
    if (!check_columns(lacenObject$datTraits, c("Sample", "Trait"), "datTraits")) {
      is_valid <- FALSE
    }
    if (!is.numeric(lacenObject$datTraits$Trait)) {
      warning('datTraits "Trait" column should be numeric.')
      is_valid <- FALSE
    }
  } else {
    warning("datTraits is NULL.")
    is_valid <- FALSE
  }

  # Validate datCounts
  if (!is.null(lacenObject$datCounts)) {
    if (!all(rownames(lacenObject$datCounts) %in% lacenObject$annotationData$gene_id)) {
      warning("All gene IDs in datCounts should be present in annotationData.")
      is_valid <- FALSE
    }
    if (!all(colnames(lacenObject$datCounts) %in% lacenObject$datTraits$Sample)) {
      warning("All sample names in datCounts should be present in datTraits.")
      is_valid <- FALSE
    }
  } else {
    warning("datCounts is NULL.")
    is_valid <- FALSE
  }

  # Validate datExpression
  if (!is.null(lacenObject$datExpression)) {
    if (!all(c("gene_id", "log2FC", "pval") %in% colnames(lacenObject$datExpression))) {
      warning('datExpression should have columns: "gene_id", "log2FC", "pval".')
      is_valid <- FALSE
    }
    if (!all(lacenObject$datExpression$gene_id %in% lacenObject$annotationData$gene_id)) {
      warning("All gene IDs in datExpression should be present in annotationData.")
      is_valid <- FALSE
    }
  } else {
    warning("datExpression is NULL.")
    is_valid <- FALSE
  }

  if (isTRUE(is_valid)) {
    genes <- lacenObject$annotationData$gene_name

    if (!requireNamespace(orgdb, quietly = TRUE)) {
      warning(sprintf("Package '%s' is required but not installed. Please install it via Bioconductor.", orgdb))
      is_valid <- FALSE
    } else {
      tryCatch({
        orgdb_obj <- get(orgdb, envir = asNamespace(orgdb))

        mapped <- AnnotationDbi::select(orgdb_obj,
                                        keys = genes,
                                        columns = c("ENTREZID"),
                                        keytype = "SYMBOL")

        positiveGenes <- sum(!is.na(mapped$ENTREZID))
        allGenes <- length(genes)
        geneRatio <- round(100 * positiveGenes / allGenes)

        if (geneRatio < 50) {
          warning("Fewer than 50% of the gene symbols were successfully mapped. Please verify that the gene symbols and the annotation database are correct; otherwise, the enrichment analysis may yield unreliable results.")
          is_valid <- FALSE
        }
      }, error = function(e) {
        warning("An error occurred during gene symbol mapping. Please verify the annotation database and gene symbols.")
        is_valid <- FALSE
      })
    }
  }


  if (is_valid) {
    message("The data in the lacen object is in the correct format.")
  }
  return(is_valid)
}
