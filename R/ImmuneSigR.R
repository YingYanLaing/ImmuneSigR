#' @title Internal function to read the built-in GMT file
#' @description Safely reads the GMT file stored in inst/extdata.
#' @return A list of gene signatures.
#' @noRd
immune_sigr_extdata <- function(file_name) {
  file_path <- system.file("extdata", file_name, package = "ImmuneSigR")
  if (nzchar(file_path)) return(file_path)

  dev_path <- file.path(getwd(), "inst", "extdata", file_name)
  if (file.exists(dev_path)) return(dev_path)

  stop("ImmuneSigR data file not found: ", file_name, call. = FALSE)
}

#' @noRd
immune_sigr_grepl <- function(pattern, x, ignore_case = TRUE, fixed = FALSE) {
  x[is.na(x)] <- ""
  if (isTRUE(fixed) && isTRUE(ignore_case)) {
    return(grepl(tolower(pattern), tolower(x), fixed = TRUE))
  }
  grepl(pattern, x, ignore.case = ignore_case, fixed = fixed)
}

#' @noRd
read_internal_gmt <- function(file_path = NULL) {
  if (is.null(file_path)) {
    file_path <- immune_sigr_extdata("cellmarker_gmt_0503.gmt")
  }
  if (!file.exists(file_path)) stop("GMT file not found: ", file_path, call. = FALSE)

  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  res <- list()
  sig_names <- character(length(lines))
  for(i in seq_along(lines)) {
    parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    sig_names[i] <- parts[1]
    genes <- unique(trimws(parts[seq.int(3, length(parts))]))
    res[[i]] <- genes[nzchar(genes) & !is.na(genes)]
  }
  keep <- nzchar(sig_names) & lengths(res) > 0
  res <- res[keep]
  sig_names <- sig_names[keep]
  names(res) <- make.unique(sig_names)
  res
}

#' @title Internal function to read the built-in metadata CSV
#' @description Safely reads the metadata CSV stored in inst/extdata.
#' @return A data frame containing metadata.
#' @noRd
read_internal_meta <- function(include_markers = FALSE) {
  file_path <- immune_sigr_extdata("cellmarker_meta.csv")
  meta <- utils::read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )

  if (!include_markers && ncol(meta) > 4) {
    marker_counts <- rowSums(!is.na(meta[, -(1:4), drop = FALSE]) & meta[, -(1:4), drop = FALSE] != "")
    meta <- meta[, 1:4, drop = FALSE]
    meta$Marker_Count <- marker_counts
  }

  meta
}

#' @title Export the built-in GMT file
#' @description Copies the internal ImmuneSigR GMT database to a specified local directory.
#' @param out_dir A character string specifying the output directory. Defaults to current working directory.
#' @param create_dir Logical. If TRUE, creates `out_dir` when it does not exist.
#' @return Invisibly returns the exported GMT file path.
#' @export
#' @examples
#' \dontrun{
#' Export_ImmuneSigR_GMT(out_dir = tempdir())
#' }
Export_ImmuneSigR_GMT <- function(out_dir = ".", create_dir = TRUE) {
  if (!dir.exists(out_dir)) {
    if (!create_dir) stop("Output directory does not exist: ", out_dir, call. = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  file_path <- immune_sigr_extdata("cellmarker_gmt_0503.gmt")
  dest_path <- file.path(out_dir, "ImmuneSigR_signatures.gmt")
  ok <- file.copy(file_path, dest_path, overwrite = TRUE)
  if (!ok) stop("Failed to export GMT file to: ", dest_path, call. = FALSE)
  message("ImmuneSigR GMT file successfully exported to: ", dest_path)
  invisible(dest_path)
}

#' @title Search the ImmuneSigR Database
#' @description Search for specific immune cell signatures based on metadata such as cell type, literature title, or PMID.
#' @param keyword A character string specifying the search term (e.g., "Macrophage"). If NULL, returns all records.
#' @param search_by A character string specifying the column to search. Options include "Cell_Type", "Title", "cell_name", "PMID".
#' @param ignore_case Logical. If TRUE, matching ignores case.
#' @param fixed Logical. If TRUE, treats `keyword` as plain text instead of a regular expression.
#' @param max_results Maximum number of rows to return.
#' @param include_markers Logical. If TRUE, includes all marker columns from the metadata CSV.
#' @return A data frame containing the search results.
#' @export
#' @examples
#' \dontrun{
#' # Search for B cell signatures
#' b_cell_info <- Search_ImmuneSigR(keyword = "B cell", search_by = "Cell_Type")
#' }
Search_ImmuneSigR <- function(keyword = NULL, search_by = "Cell_Type", ignore_case = TRUE,
                              fixed = FALSE, max_results = Inf, include_markers = FALSE) {
  meta <- read_internal_meta(include_markers = include_markers)
  if (!search_by %in% colnames(meta)) {
    stop(
      "Invalid 'search_by' column. Please choose from: ",
      paste(colnames(meta), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(keyword) || !nzchar(keyword)) {
    res <- meta
  } else {
    hits <- immune_sigr_grepl(keyword, meta[[search_by]], ignore_case = ignore_case, fixed = fixed)
    res <- meta[hits, , drop = FALSE]
  }

  if (is.finite(max_results) && nrow(res) > max_results) {
    res <- res[seq_len(max_results), , drop = FALSE]
  }

  if(nrow(res) == 0) {
    message("No matching records found.")
  } else {
    message(sprintf("Found %d matching cell signatures.", nrow(res)))
  }
  return(res)
}

#' @title Get Marker Genes for Specific Cell Types
#' @description Retrieve the marker gene lists for specified immune cell subpopulations.
#' @param cell_type A character vector specifying cell type keywords. If NULL, returns the entire database.
#' @param ignore_case Logical. If TRUE, matching ignores case.
#' @param fixed Logical. If TRUE, treats `cell_type` as plain text instead of a regular expression.
#' @param min_genes Minimum number of genes required for a returned signature.
#' @param gmt_file Optional path to a custom GMT file.
#' @return A list containing the marker genes.
#' @export
#' @examples
#' \dontrun{
#' # Get all markers
#' all_markers <- Get_Markers()
#' # Get specific markers
#' t_cell_markers <- Get_Markers("T cell")
#' }
Get_Markers <- function(cell_type = NULL, ignore_case = TRUE, fixed = FALSE,
                        min_genes = 1, gmt_file = NULL) {
  sigs <- read_internal_gmt(gmt_file)
  sigs <- sigs[lengths(sigs) >= min_genes]
  if (is.null(cell_type) || all(!nzchar(cell_type))) return(sigs)

  matched <- unique(unlist(lapply(cell_type[nzchar(cell_type)], function(query) {
    names(sigs)[immune_sigr_grepl(query, names(sigs), ignore_case = ignore_case, fixed = fixed)]
  }), use.names = FALSE))
  if(length(matched) == 0) stop(paste("No signatures found matching:", cell_type))
  sigs[matched]
}

#' @noRd
immune_sigr_prepare_signatures <- function(target_cells = NULL, min_genes = 5,
                                           gmt_file = NULL, object_genes = NULL) {
  sigs <- Get_Markers(target_cells, min_genes = min_genes, gmt_file = gmt_file)
  if (!is.null(object_genes)) {
    sigs <- lapply(sigs, intersect, y = object_genes)
    sigs <- sigs[lengths(sigs) >= min_genes]
  }

  if (length(sigs) == 0) {
    stop(
      "No matching cell signatures with at least ",
      min_genes,
      " usable genes were found.",
      call. = FALSE
    )
  }

  names(sigs) <- paste0("ImmuneSigR_", names(sigs))
  sigs
}

#' @noRd
immune_sigr_score_matrix <- function(expr, signatures, score_name = "_score") {
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) {
    stop("Expression matrix must have gene names as row names.", call. = FALSE)
  }

  scores <- vapply(signatures, function(genes) {
    colMeans(expr[genes, , drop = FALSE])
  }, numeric(ncol(expr)))

  if (is.null(dim(scores))) {
    scores <- matrix(scores, ncol = 1)
  }

  colnames(scores) <- paste0(names(signatures), score_name)
  rownames(scores) <- colnames(expr)
  as.data.frame(scores, check.names = FALSE)
}

#' @noRd
immune_sigr_rank_score_one_cell <- function(values, signature_indices, max_rank) {
  ranks <- rank(-values, ties.method = "average", na.last = "keep")
  ranks[is.na(ranks)] <- max_rank
  ranks <- pmin(ranks, max_rank)

  vapply(signature_indices, function(idx) {
    sig_ranks <- ranks[idx]
    n_genes <- length(sig_ranks)
    best_sum <- n_genes * (n_genes + 1) / 2
    worst_sum <- n_genes * max_rank
    if (worst_sum <= best_sum) return(NA_real_)

    score <- 1 - ((sum(sig_ranks) - best_sum) / (worst_sum - best_sum))
    max(0, min(1, score))
  }, numeric(1))
}

#' @noRd
immune_sigr_score_matrix_rank <- function(expr, signatures, score_name = "_score",
                                          max_rank = 1500) {
  expr <- as.matrix(expr)
  if (is.null(rownames(expr))) {
    stop("Expression matrix must have gene names as row names.", call. = FALSE)
  }

  max_rank <- min(as.integer(max_rank), nrow(expr))
  if (is.na(max_rank) || max_rank < 2) {
    stop("max_rank must be at least 2 after considering the number of genes.", call. = FALSE)
  }

  signature_indices <- lapply(signatures, match, table = rownames(expr))
  signature_indices <- lapply(signature_indices, function(idx) idx[!is.na(idx)])
  signature_indices <- signature_indices[lengths(signature_indices) > 0]

  scores <- apply(
    expr,
    2,
    immune_sigr_rank_score_one_cell,
    signature_indices = signature_indices,
    max_rank = max_rank
  )

  if (is.null(dim(scores))) {
    scores <- matrix(scores, ncol = 1)
  }

  scores <- t(scores)
  colnames(scores) <- paste0(names(signature_indices), score_name)
  rownames(scores) <- colnames(expr)
  as.data.frame(scores, check.names = FALSE)
}

#' @title Create a Custom GMT File
#' @description Converts a user-defined R list of markers into a standard GMT file format.
#' @param marker_list A list where names are cell types and elements are character vectors of genes.
#' @param file_name A character string specifying the output file name.
#' @return Invisibly returns the created GMT file path.
#' @export
#' @examples
#' \dontrun{
#' my_markers <- list(Custom_T = c("CD3D", "CD8A"), Custom_B = c("CD19", "MS4A1"))
#' Create_Custom_GMT(my_markers, file_name = file.path(tempdir(), "custom.gmt"))
#' }
Create_Custom_GMT <- function(marker_list, file_name = "My_Custom_Signatures.gmt") {
  if(!is.list(marker_list) || is.null(names(marker_list))) {
    stop("marker_list must be a named list of character vectors.", call. = FALSE)
  }
  if (any(!nzchar(names(marker_list)))) {
    stop("All marker_list entries must have non-empty names.", call. = FALSE)
  }

  con <- file(file_name, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  for (name in names(marker_list)) {
    genes <- unique(trimws(as.character(marker_list[[name]])))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    if (length(genes) == 0) {
      warning("Skipping empty marker set: ", name, call. = FALSE)
      next
    }
    writeLines(paste(c(name, "Custom_Signature", genes), collapse = "\t"), con = con)
  }
  message("Custom GMT file created successfully: ", file_name)
  invisible(file_name)
}

#' @title Score Expression Data using ImmuneSigR
#' @description Scores expression matrices directly with dependency-free rank or mean methods.
#' @param expr An expression matrix or data frame with genes in rows and cells in columns.
#' @param target_cells A character vector of target cell subpopulations to score. If NULL, scores all signatures.
#' @param min_genes Minimum number of genes required after filtering to genes present in the object.
#' @param gmt_file Optional path to a custom GMT file.
#' @param score_name Suffix used for generated score columns.
#' @param method Scoring method. "rank" is a dependency-free UCell-like rank score; "mean" uses average expression; "auto" uses "rank".
#' @param max_rank Maximum rank considered by the rank scoring method.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param ... Reserved for future extensions.
#' @return A data frame of scores with cells in rows and signatures in columns.
#' @export
#' @examples
#' \dontrun{
#' scores <- Score_ImmuneSigR(expr_matrix, target_cells = c("B cell", "Plasma cell"), method = "rank")
#' }
Score_ImmuneSigR <- function(expr, target_cells = NULL, min_genes = 5,
                             gmt_file = NULL, score_name = "_score",
                             method = c("auto", "rank", "mean"),
                             max_rank = 1500,
                             verbose = TRUE, ...) {
  method <- match.arg(method)
  is_matrix_like <- is.matrix(expr) || is.data.frame(expr)

  if (method == "auto") {
    method <- "rank"
  }

  if (!is_matrix_like) {
    stop(
      "Score_ImmuneSigR() requires a matrix or data frame. For Seurat objects, pass an expression matrix extracted from the object.",
      call. = FALSE
    )
  }

  sigs <- immune_sigr_prepare_signatures(
    target_cells = target_cells,
    min_genes = min_genes,
    gmt_file = gmt_file,
    object_genes = rownames(expr)
  )

  if (isTRUE(verbose)) {
    message(sprintf("Scoring %d signatures using %s scoring...", length(sigs), method))
  }
  if (method == "rank") {
    return(immune_sigr_score_matrix_rank(expr, sigs, score_name = score_name, max_rank = max_rank))
  }

  immune_sigr_score_matrix(expr, sigs, score_name = score_name)
}

# Backward-compatible aliases for scripts written against the original package.

#' @rdname Export_ImmuneSigR_GMT
#' @export
Export_GMT <- Export_ImmuneSigR_GMT

#' @rdname Search_ImmuneSigR
#' @export
Search_CellSigR <- Search_ImmuneSigR

#' @rdname Score_ImmuneSigR
#' @export
Score_CellSigR <- Score_ImmuneSigR
