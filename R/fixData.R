#' Process and Format Input Data for a lacen Object
#'
#' This function takes raw count data, differential expression results, sample
#' traits, and annotation data, then processes and formats them into a
#' standardized `lacen` object for downstream analysis.
#'
#' @details
#' During processing, genes may be filtered out if they are not present in the
#' annotation data. Additionally, the column names of the input data frames
#' may be standardized to ensure consistency.
#'
#' @param datCounts Data counts dataframe with genes as rows and samples as columns.
#' The gene/transcript IDs should be in the row names and sample names/IDs in the column names.
#' @param datExpression Differential expression dataframe. Can be from Limma output or a dataframe with columns "gene_id", "log2FC", and "pval".
#' @param datTraits Conditions/Traits data. A two-column dataframe with sample IDs in "Sample" and condition codes in "Trait".
#' @param annotationData Annotation dataframe with "gene_id" and "gene_name" columns.
#' @param ncAnnotation Subset of "annotationData" containing only long non-coding RNAs.
#' @param verbose A logical value (`TRUE` or `FALSE`) indicating whether to print
#'   progress messages to the console. Defaults to `TRUE`.
#'
#' @return A 'lacen' S3 object containing the cleaned,
#'   formatted, and validated data, ready for further analysis.
#'
#' @export

fixData <- function(annotationData = annotation_data,
                    datCounts      = rawExpressionData,
                    datExpression  = expressionDGEData,
                    datTraits      = metadata,
                    ncAnnotation   = ncAnnotation,
                    verbose        = TRUE) {

  # --------------------------- tiny logger ---------------------------
  .emit <- function(level = "info", msg = "") {
    if (!isTRUE(verbose)) return(invisible())
    prefix <- switch(level,
                     step   = "▶",
                     change = "•",
                     note   = "•",
                     warn   = "!",
                     done   = "✓",
                     info   = "•")
    cat(sprintf("%s %s\n", prefix, msg))
  }
  .log <- list()
  .add <- function(key, value) .log[[length(.log)+1L]] <<- list(key = key, value = value)

  # --------------------------- helpers ---------------------------
  .lower <- function(x) tolower(gsub("\\s+", "", x))
  .is_char <- function(x) is.character(x) || is.factor(x)
  .as_char <- function(x) { if (is.factor(x)) as.character(x) else as.character(x) }
  .notna   <- function(x) x[!is.na(x)]
  .strip_versions <- function(ids) sub("\\.\\d+$", "", ids)

  common_gene_symbols <- c(
    "TP53","EGFR","GAPDH","ACTB","BRCA1","BRCA2","MYC","FOS","JUN","PTEN",
    "VEGFA","IL6","TNF","ALB","INS","MKI67","CDKN1A","APOE","PPIA","RPLP0","RPS18",
    "Trp53","Actb","Gapdh","Myc","Pten","Vegfa","Il6","Tnf","Alb","Ins1","Mki67","Cdkn1a"
  )
  id_patterns <- c("ENSG\\d+", "ENSMUSG\\d+", "ENS\\w+G\\d+", "FBgn\\d+", "GeneID:\\d+")
  name_syns   <- c("genename","gene_name","genesymbol","symbol","external_gene_name","hgnc_symbol","gene")
  id_syns     <- c("ensembl","geneid","gene_id","ensemblgene","ensembl_id","entrezid","entrez_id")

  .looks_like_idvec <- function(v) {
    v <- .as_char(v)
    any(grepl(paste(id_patterns, collapse="|"), v, perl=TRUE, ignore.case=FALSE)) ||
      (mean(grepl("^[0-9]+$", v)) > 0.6)
  }
  .looks_like_namevec <- function(v) {
    v <- .as_char(v)
    any(v %in% common_gene_symbols) ||
      mean(grepl("^[A-Za-z][A-Za-z0-9\\.-]*$", v)) > 0.7
  }

  .guess_anno_cols <- function(df, label="annotationData") {
    nms  <- names(df); nmsl <- .lower(nms)
    name_hit <- which(nmsl %in% name_syns | grepl("name|symbol", nmsl))
    id_hit   <- which(nmsl %in% id_syns   | grepl("ensembl|geneid|entrez", nmsl))
    gene_name_col <- if (length(name_hit)) name_hit[1] else NA_integer_
    gene_id_col   <- if (length(id_hit))   id_hit[1]   else NA_integer_

    if (is.na(gene_name_col) || is.na(gene_id_col)) {
      for (i in seq_along(nms)) {
        v <- df[[i]]
        if (is.na(gene_name_col) && .is_char(v) && .looks_like_namevec(v)) gene_name_col <- i
        if (is.na(gene_id_col)   && .is_char(v) && .looks_like_idvec(v))   gene_id_col   <- i
      }
    }
    if ((is.na(gene_name_col) || is.na(gene_id_col)) && ncol(df) == 2) {
      s1 <- .looks_like_namevec(df[[1]]); s2 <- .looks_like_namevec(df[[2]])
      if (is.na(gene_name_col) && (s1 || !s2)) gene_name_col <- 1L
      if (is.na(gene_name_col) && (!s1 && s2)) gene_name_col <- 2L
      if (is.na(gene_id_col))   gene_id_col <- setdiff(1:2, gene_name_col)
    }
    if (is.na(gene_name_col) || is.na(gene_id_col)) {
      stop(sprintf("%s: could not infer 'gene_id' and 'gene_name' columns.", label))
    }

    out <- data.frame(
      gene_id   = .as_char(df[[gene_id_col]]),
      gene_name = .as_char(df[[gene_name_col]]),
      stringsAsFactors = FALSE
    )
    if (anyDuplicated(out$gene_id)) {
      dups <- sum(duplicated(out$gene_id))
      .emit("note", sprintf("%s: %d duplicated gene_id -> keeping first occurrence.", label, dups))
      .add(paste0(label, ".duplicate_gene_id_removed"), dups)
      out <- out[!duplicated(out$gene_id), , drop=FALSE]
    }
    out
  }

  .guess_traits <- function(dt) {
    nms  <- names(dt); nmsl <- .lower(nms)
    sample_hit <- which(nmsl %in% c("sample","sampleid","id","run","barcode","cell","library","sample_name","samplename"))
    trait_hit  <- which(nmsl %in% c("trait","condition","group","status","phenotype","treatment","class","cluster"))
    Sample <- if (length(sample_hit)) .as_char(dt[[ sample_hit[1] ]]) else NULL
    Trait  <- if (length(trait_hit))  dt[[ trait_hit[1]  ]]          else NULL

    if (is.null(Sample) || is.null(Trait)) {
      uniq_counts <- vapply(dt, function(x) length(unique(.notna(.as_char(x)))), integer(1))
      cand_trait <- which(uniq_counts <= 10 & uniq_counts >= 2)
      if (length(cand_trait)) {
        trait_idx <- cand_trait[ which.min(uniq_counts[cand_trait]) ]
        Trait  <- dt[[trait_idx]]
        cand_sample <- which.max(uniq_counts)
        Sample <- .as_char(dt[[cand_sample]])
        .emit("note", sprintf("datTraits: inferred Trait='%s' (<=10 classes), Sample='%s'.",
                              names(dt)[trait_idx], names(dt)[cand_sample]))
      } else {
        Sample <- .as_char(dt[[1]]); Trait <- dt[[2]]
        .emit("note", "datTraits: using first two columns as Sample/Trait (no clear candidates).")
      }
    }
    before <- Sample
    Sample <- gsub("[^A-Za-z0-9]", "", Sample)
    changed <- sum(before != Sample)
    if (changed) {
      .emit("change", sprintf("datTraits$Sample: sanitized %d sample name(s) (removed non-alphanumerics).", changed))
      .add("traits.sample_names_sanitized", changed)
    }
    Trait <- as.numeric(as.factor(.as_char(Trait)))
    data.frame(Sample = Sample, Trait = Trait, stringsAsFactors = FALSE)
  }

  .counts_all_numeric <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    pre_non_numeric <- names(df)[!vapply(df, is.numeric, logical(1))]
    for (j in seq_len(ncol(df))) {
      if (is.factor(df[[j]])) df[[j]] <- as.character(df[[j]])
      suppressWarnings(df[[j]] <- as.numeric(df[[j]]))
    }
    n_na <- sum(is.na(as.matrix(df)))
    if (n_na > 0) df[is.na(df)] <- 0
    list(data = df, na_introduced = n_na, coerced_cols = pre_non_numeric)
  }

  # --------------------------- header ---------------------------
  .emit("step", "Checking/repairing inputs …")
  .emit("info", sprintf("annotationData: %d rows", NROW(annotationData)))
  .emit("info", sprintf("datCounts: %d genes x %d samples", NROW(datCounts), NCOL(datCounts)))
  .emit("info", sprintf("datExpression: %d rows", NROW(datExpression)))
  .emit("info", sprintf("datTraits: %d rows", NROW(datTraits)))

  # -------------------- 0) Ensure/repair annotationData --------------------
  if (!all(c("gene_id","gene_name") %in% names(annotationData))) {
    .emit("step", "Inferring gene_id/gene_name columns in annotationData …")
    annotationData <- .guess_anno_cols(annotationData, "annotationData")
    .add("annotationData.inferred_columns", TRUE)
  } else {
    annotationData$gene_id   <- .as_char(annotationData$gene_id)
    annotationData$gene_name <- .as_char(annotationData$gene_name)
  }

  # -------------------- 0b) Ensure/repair ncAnnotation (if any) -----------
  if (!is.null(ncAnnotation)) {
    if (!all(c("gene_id","gene_name") %in% names(ncAnnotation))) {
      .emit("step", "Inferring gene_id/gene_name columns in ncAnnotation …")
      ncAnnotation <- .guess_anno_cols(ncAnnotation, "ncAnnotation")
      .add("ncAnnotation.inferred_columns", TRUE)
    } else {
      ncAnnotation$gene_id   <- .as_char(ncAnnotation$gene_id)
      ncAnnotation$gene_name <- .as_char(ncAnnotation$gene_name)
    }
  }

  # ------------ 0c) Special rule about gene_id if names look like symbols --
  if (!is.null(rownames(datCounts))) {
    rn <- rownames(datCounts)
    if (sum(rn %in% annotationData$gene_name) > length(rn)/10) {
      .emit("note", "Many count rownames match gene_name; using gene_name as gene_id in annotation tables.")
      annotationData$gene_id <- annotationData$gene_name
      if (!is.null(ncAnnotation)) ncAnnotation$gene_id <- ncAnnotation$gene_name
      .add("annotationData.use_gene_name_as_id", TRUE)
    }
  }

  # -------------------- 1) Traits: enforce Sample + Trait ------------------
  .emit("step", "Validating datTraits (Sample/Trait) …")
  if (!all(c("Sample","Trait") %in% names(datTraits))) {
    datTraits <- .guess_traits(datTraits)
    .add("traits.inferred", TRUE)
  } else {
    before <- datTraits$Sample
    datTraits$Sample <- gsub("[^A-Za-z0-9]", "", .as_char(datTraits$Sample))
    changed <- sum(before != datTraits$Sample)
    if (changed) {
      .emit("change", sprintf("datTraits$Sample: sanitized %d sample name(s).", changed))
      .add("traits.sample_names_sanitized", changed)
    }
    datTraits$Trait  <- as.numeric(as.factor(.as_char(datTraits$Trait)))
  }

  # -------------------- 2) Counts: numeric + clean sample names -----------
  .emit("step", "Processing datCounts …")
  dc <- as.data.frame(datCounts, check.names = FALSE)
  if ("gene_id" %in% names(dc)) {
    if (anyDuplicated(dc$gene_id)) stop("datCounts$gene_id contains duplicates.")
    rownames(dc) <- .as_char(dc$gene_id); dc$gene_id <- NULL
  } else if (is.null(rownames(dc))) {
    stop("datCounts must have rownames as gene IDs or contain a 'gene_id' column.")
  }

  before_cols <- colnames(dc)
  colnames(dc) <- gsub("[^A-Za-z0-9]", "", colnames(dc))
  changed_cols <- sum(before_cols != colnames(dc))
  if (changed_cols) {
    .emit("change", sprintf("datCounts colnames: sanitized %d sample name(s).", changed_cols))
    .add("counts.sample_names_sanitized", changed_cols)
  }

  cn <- .counts_all_numeric(dc)
  dc <- cn$data
  if (length(cn$coerced_cols)) {
    .emit("change", sprintf("datCounts: coerced %d non-numeric column(s) to numeric.", length(cn$coerced_cols)))
    .add("counts.columns_coerced_to_numeric", length(cn$coerced_cols))
  }
  if (cn$na_introduced > 0) {
    .emit("note", sprintf("datCounts: replaced %d NA value(s) created during coercion with 0.", cn$na_introduced))
    .add("counts.na_replaced_with_zero", cn$na_introduced)
  }

  # -------------------- 3) Expression: flexible mapping -------------------
  .emit("step", "Processing datExpression …")
  de <- as.data.frame(datExpression, check.names = FALSE)
  nms_de  <- names(de); nms_del <- .lower(nms_de)

  if (!any(nms_del %in% c("gene_id","geneid","gene_name","genename"))) {
    sym_hit <- which(nms_del %in% c("genesymbol","symbol","hgnc_symbol","external_gene_name"))
    if (length(sym_hit)) {
      names(de)[sym_hit[1]] <- "gene_name"
      .emit("change", sprintf("datExpression: renamed '%s' -> 'gene_name'.", nms_de[sym_hit[1]]))
      .add("expression.rename_to_gene_name", nms_de[sym_hit[1]])
    }
  }
  if (!("gene_id" %in% names(de)) && !("gene_name" %in% names(de)) && !is.null(rownames(de))) {
    rn <- rownames(de)
    if (any(grepl(paste(id_patterns, collapse="|"), rn, perl=TRUE))) {
      de$gene_id <- rn
      .emit("note", "datExpression: used rownames as gene_id.")
    } else {
      de$gene_name <- rn
      .emit("note", "datExpression: used rownames as gene_name.")
    }
  }

  if (!("log2FC" %in% names(de))) {
    if ("logFC" %in% names(de)) {
      de$log2FC <- de$logFC
      .emit("change", "datExpression: mapped 'logFC' -> 'log2FC'.")
    } else if ("log2fc" %in% nms_del) {
      de$log2FC <- de[[ which(nms_del == "log2fc")[1] ]]
      .emit("change", "datExpression: normalized 'log2fc' -> 'log2FC'.")
    } else if ("logfc" %in% nms_del) {
      de$log2FC <- de[[ which(nms_del == "logfc")[1] ]]
      .emit("change", "datExpression: normalized 'logfc' -> 'log2FC'.")
    } else {
      stop("datExpression must provide 'log2FC' or 'logFC'.")
    }
  }

  p_hits <- c("fdr","adjpval","adj.p.val","padj","pvalue","pval")
  have_p <- intersect(p_hits, nms_del)
  if (length(have_p)) {
    de$pval <- de[[ which(nms_del == have_p[1])[1] ]]
    if (have_p[1] != "pval") .emit("change", sprintf("datExpression: mapped '%s' -> 'pval'.", have_p[1]))
  } else {
    stop("datExpression must provide FDR/adj.P.Val/padj or PValue/pval (mapped to 'pval').")
  }

  if (!("gene_id" %in% names(de))) {
    map_df <- unique(annotationData[, c("gene_id","gene_name")])
    de$gene_name <- .as_char(de$gene_name)
    pre_n <- nrow(de)
    # Suppress any messages the merge or downstream might trigger
    de <- suppressMessages(merge(de, map_df, by = "gene_name", all.x = TRUE, sort = FALSE))
    n_na_gid <- sum(is.na(de$gene_id))
    if (n_na_gid > 0) {
      .emit("change", sprintf("datExpression: removed %d row(s) with gene_name not in annotation.", n_na_gid))
      .add("expression.rows_removed_no_gene_id", n_na_gid)
      de <- de[!is.na(de$gene_id), , drop = FALSE]
    }
    tab <- table(de$gene_name)
    if (any(tab > 1)) {
      g_multi <- sum(tab > 1)
      .emit("note", sprintf("datExpression: %d gene_name matched multiple gene_id (kept all).", g_multi))
      .add("expression.multiple_ids_per_name", g_multi)
    }
    post_n <- nrow(de)
    .emit("info", sprintf("datExpression: %d -> %d row(s) after mapping to gene_id.", pre_n, post_n))
  }
  de <- de[, c("gene_id","log2FC","pval"), drop = FALSE]
  de$gene_id <- .as_char(de$gene_id)

  # -------------------- 4) Harmonize IDs & overlaps -----------------------
  .emit("step", "Harmonizing IDs and filtering by annotation …")
  n_ver_counts <- sum(grepl("\\.\\d+$", rownames(dc)))
  n_ver_anno   <- sum(grepl("\\.\\d+$", annotationData$gene_id))
  n_ver_de     <- sum(grepl("\\.\\d+$", de$gene_id))
  if (n_ver_counts + n_ver_anno + n_ver_de > 0) {
    .emit("change", sprintf("Removed version suffixes from IDs (counts:%d, annotation:%d, DE:%d).",
                            n_ver_counts, n_ver_anno, n_ver_de))
    .add("ids.version_suffix_removed", c(counts=n_ver_counts, annotation=n_ver_anno, de=n_ver_de))
  }
  base_ids <- .strip_versions(rownames(dc))
  annotationData$gene_id <- .strip_versions(annotationData$gene_id)
  if (!is.null(ncAnnotation)) ncAnnotation$gene_id <- .strip_versions(ncAnnotation$gene_id)
  de$gene_id <- .strip_versions(de$gene_id)

  dup_n <- sum(duplicated(base_ids))
  if (dup_n > 0) {
    pre_g <- nrow(dc)
    dc <- as.data.frame(rowsum(as.matrix(dc), group = base_ids, reorder = FALSE))
    .emit("change", sprintf("datCounts: collapsed %d duplicated gene_id(s) after version strip (%d -> %d genes).",
                            dup_n, pre_g, nrow(dc)))
    .add("counts.collapsed_duplicates", dup_n)
  } else {
    rownames(dc) <- base_ids
  }

  pre_counts_genes <- nrow(dc)
  dc <- dc[rownames(dc) %in% annotationData$gene_id, , drop = FALSE]
  drop_counts <- pre_counts_genes - nrow(dc)
  if (drop_counts > 0) {
    .emit("change", sprintf("datCounts: filtered %d gene(s) not present in annotation.", drop_counts))
    .add("counts.filtered_not_in_annotation", drop_counts)
  }

  pre_de_genes <- nrow(de)
  de <- de[de$gene_id %in% annotationData$gene_id, , drop = FALSE]
  drop_de <- pre_de_genes - nrow(de)
  if (drop_de > 0) {
    .emit("change", sprintf("datExpression: filtered %d row(s) not present in annotation.", drop_de))
    .add("expression.filtered_not_in_annotation", drop_de)
  }

  if (nrow(dc) == 0L) stop("After filtering, datCounts has 0 genes matching annotationData.")
  if (nrow(de) == 0L) stop("After filtering, datExpression has 0 rows matching annotationData.")

  # -------------------- 5) Align samples with traits ----------------------
  .emit("step", "Aligning samples with traits …")
  before_dt_n <- nrow(datTraits)
  sample_intersect <- intersect(colnames(dc), datTraits$Sample)
  rm_counts_samples <- ncol(dc) - length(sample_intersect)
  rm_traits_rows    <- before_dt_n - length(sample_intersect)
  if (!length(sample_intersect)) stop("No overlap between datCounts columns and datTraits$Sample after sanitization.")

  if (rm_counts_samples > 0)
    .emit("change", sprintf("Dropped %d sample(s) from datCounts not in datTraits.", rm_counts_samples))
  if (rm_traits_rows > 0)
    .emit("change", sprintf("Ignored %d sample(s) in datTraits not in datCounts.", rm_traits_rows))
  .add("samples.dropped_from_counts", rm_counts_samples)
  .add("samples.ignored_in_traits", rm_traits_rows)

  dc <- dc[, sample_intersect, drop = FALSE]
  datTraits <- datTraits[match(sample_intersect, datTraits$Sample), , drop = FALSE]

  # -------------------- 6) Make datExpression and counts share genes ------
  .emit("step", "Synchronizing gene sets between counts and DE …")
  common_genes <- intersect(rownames(dc), de$gene_id)
  if (!length(common_genes)) {
    .emit("warn", "No overlap between datCounts genes and datExpression genes.")
  } else {
    drop_de2 <- nrow(de) - length(common_genes)
    drop_dc2 <- nrow(dc) - length(common_genes)
    if (drop_de2 > 0)
      .emit("change", sprintf("datExpression: restricted to %d shared genes (dropped %d).",
                              length(common_genes), drop_de2))
    if (drop_dc2 > 0)
      .emit("change", sprintf("datCounts: restricted to %d shared genes (dropped %d).",
                              length(common_genes), drop_dc2))
    .add("genes.shared_between_counts_and_de", length(common_genes))
    de <- de[match(common_genes, de$gene_id), , drop = FALSE]
    dc <- dc[common_genes, , drop = FALSE]
  }

  # -------------------- 7) Build object & check (suppress noisy messages) ----
  .emit("step", "Building lacen object and running checks …")
  lacenObject <- suppressMessages(initLacen(annotationData = annotationData,
                                            datCounts      = dc,
                                            datExpression  = de,
                                            datTraits      = datTraits,
                                            ncAnnotation   = ncAnnotation))
  ok <- suppressMessages(tryCatch(checkData(lacenObject), error = function(e) e))
  if (isTRUE(ok)) {
    .emit("done", "checkData passed.")
  } else {
    stop("checkData failed; please inspect your inputs. Details: ", as.character(ok))
  }

  # -------------------- summary & audit trail -----------------------------
  .emit("step", "Summary:")
  .emit("info", sprintf("Genes (counts): %d", nrow(dc)))
  .emit("info", sprintf("Genes (DE): %d", nrow(de)))
  .emit("info", sprintf("Samples: %d", ncol(dc)))

  attr(lacenObject, "fixData_log") <- .log
  invisible(lacenObject)
}
