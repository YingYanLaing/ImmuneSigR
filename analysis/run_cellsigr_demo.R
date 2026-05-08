# =====================================================================
# ImmuneSigR Comprehensive Local Demonstration and Validation Script
# =====================================================================
# This script demonstrates the core workflow of ImmuneSigR. 
# It includes database retrieval, GMT management, base algorithm tests,
# and real-world scRNA-seq validation suitable for scientific publication.
# =====================================================================

source("R/ImmuneSigR.R")

# 1. Setup Output Directory
out_dir <- file.path(getwd(), "analysis", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Starting ImmuneSigR demo. Outputs will be saved to:", out_dir, "\n\n")

# ---------------------------------------------------------------------
# PART 1: Database Exploration & Search
# ---------------------------------------------------------------------
cat("--- PART 1: Database Exploration ---\n")
all_markers <- Get_Markers()
markers_min5 <- Get_Markers(min_genes = 5)

# Retrieve basic records for testing
b_cell_records <- Search_ImmuneSigR("B cell", search_by = "Cell_Type", fixed = TRUE)
cd8_records <- Search_ImmuneSigR("CD8+", search_by = "cell_name", fixed = TRUE)
t_nk_markers <- Get_Markers(c("T cell", "NK cell"), min_genes = 5)

cat("Successfully retrieved markers and records.\n\n")

# ---------------------------------------------------------------------
# PART 2: GMT File Management
# ---------------------------------------------------------------------
cat("--- PART 2: GMT Exports ---\n")
exported_gmt <- Export_ImmuneSigR_GMT(out_dir)
custom_gmt <- Create_Custom_GMT(
  list(
    Custom_T = c("CD3D", "CD3E", "CD8A"),
    Custom_B = c("CD19", "MS4A1", "CD79A")
  ),
  file_name = file.path(out_dir, "custom_demo_signatures.gmt")
)
cat("Exported built-in GMT to:", exported_gmt, "\n")
cat("Created custom GMT at:", custom_gmt, "\n\n")

# ---------------------------------------------------------------------
# PART 3: Base R Scoring (Dummy Matrix Validation)
# ---------------------------------------------------------------------
cat("--- PART 3: Base R Scoring (Dummy Data for CRAN Tests) ---\n")
# Using original dummy data to test base algorithms independently without heavy dependencies
demo_genes <- unique(unlist(Get_Markers(c("B cell", "T cell"), min_genes = 5)[1:8]))
demo_genes <- demo_genes[seq_len(min(120, length(demo_genes)))]
set.seed(1)
expr_matrix_dummy <- matrix(
  stats::rpois(length(demo_genes) * 12, lambda = 2),
  nrow = length(demo_genes),
  dimnames = list(demo_genes, paste0("cell_", seq_len(12)))
)

matrix_rank_scores <- Score_ImmuneSigR(
  expr = expr_matrix_dummy, target_cells = c("B cell", "T cell"), min_genes = 5, method = "rank"
)
matrix_mean_scores <- Score_ImmuneSigR(
  expr = expr_matrix_dummy, target_cells = c("B cell", "T cell"), min_genes = 5, method = "mean"
)
cat("Successfully scored dummy matrix using rank and mean methods.\n\n")

# ---------------------------------------------------------------------
# PART 4: Real Single-Cell Validation Workflow (Seurat Integration)
# ---------------------------------------------------------------------
cat("--- PART 4: Real Single-Cell Validation ---\n")
# CRAN Standard: Check for heavy external packages before running
if (requireNamespace("Seurat", quietly = TRUE) && 
    requireNamespace("SeuratData", quietly = TRUE) && 
    requireNamespace("ggplot2", quietly = TRUE)) {
  
  suppressWarnings({
    library(Seurat)
    library(SeuratData)
    library(ggplot2)
    
    cat("Loading pbmc3k dataset for validation...\n")
    data("pbmc3k")
    pbmc <- UpdateSeuratObject(pbmc3k)
    pbmc <- NormalizeData(pbmc, verbose = FALSE)
    pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
    pbmc <- ScaleData(pbmc, verbose = FALSE)
    pbmc <- RunPCA(pbmc, verbose = FALSE)
    pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
    
    # 4.1 Funnel Search Strategy (Scientific Logic)
    cat("\nExecuting funnel search strategy for targeting...\n")
    all_db_records <- Search_ImmuneSigR()
    cat("Step 1: Available major cell lineages (Cell_Type):\n")
    print(unique(all_db_records$Cell_Type))
    
    cat("\nStep 2 & 3: Filtering for Plasma cells within the B cell pool...\n")
    plasma_records <- b_cell_records[grepl("Plasma", b_cell_records$cell_name, ignore.case = TRUE), ]
    cat("Top candidates based on literature sources:\n")
    print(head(plasma_records[, c("cell_name", "Title", "PMID")]))
    
    # 4.2 Extract matrix and score using precise targets from literature
    real_targets <- c("Plasma cell", "Conventional Plasma cells") 
    expr_matrix_real <- as.matrix(pbmc[["RNA"]]$data)
    cat("\nScoring pbmc3k using targeted Plasma signatures...\n")
    real_scores <- Score_ImmuneSigR(expr_matrix_real, target_cells = real_targets, min_genes = 5, method = "rank")
    
    # 4.3 Add metadata and plot results
    pbmc <- AddMetaData(pbmc, metadata = real_scores)
    score_cols <- colnames(real_scores)
    
    if(length(score_cols) >= 2) {
      p_umap <- FeaturePlot(pbmc, features = score_cols[1:2], ncol = 2, pt.size = 0.8) & 
        scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))
      
      # Save Vector PDF for Scientific Publication
      pdf_path <- file.path(out_dir, "ImmuneSigR_Validation_UMAP.pdf")
      ggsave(filename = pdf_path, plot = p_umap, width = 12, height = 5, device = "pdf")
      
      # Save High-DPI PNG for GitHub Display
      png_path <- file.path(out_dir, "ImmuneSigR_Validation_UMAP.png")
      ggsave(filename = png_path, plot = p_umap, width = 12, height = 5, dpi = 300)
      
      cat("Validation outputs successfully saved to:", out_dir, "(PDF and PNG formats)\n\n")
    }
  })
} else {
  cat("Seurat, SeuratData, or ggplot2 are not installed. Skipping real data UMAP generation.\n\n")
}

# ---------------------------------------------------------------------
# PART 5: Summary Generation and CSV Exports
# ---------------------------------------------------------------------
cat("--- PART 5: Exporting Summary Statistics ---\n")
meta <- ImmuneSigR:::read_internal_meta() 
signature_lengths <- lengths(all_markers)

summary_df <- data.frame(
  metric = c("total_signatures", "signatures_min_5_genes", "b_cell_records", 
             "cd8_cell_name_records", "t_or_nk_signatures_min_5_genes", 
             "matrix_demo_cells", "matrix_rank_demo_score_columns", 
             "matrix_mean_demo_score_columns", "min_signature_genes", 
             "median_signature_genes", "max_signature_genes"),
  value = c(length(all_markers), length(markers_min5), nrow(b_cell_records), 
            nrow(cd8_records), length(t_nk_markers), nrow(matrix_rank_scores), 
            ncol(matrix_rank_scores), ncol(matrix_mean_scores), 
            min(signature_lengths), stats::median(signature_lengths), max(signature_lengths))
)

cell_type_counts <- as.data.frame(sort(table(meta$Cell_Type), decreasing = TRUE))
names(cell_type_counts) <- c("Cell_Type", "Signature_Count")

# Export all CSVs (CRAN testing artifacts)
utils::write.csv(summary_df, file.path(out_dir, "immunesigr_demo_summary.csv"), row.names = FALSE)
utils::write.csv(cell_type_counts, file.path(out_dir, "cell_type_counts.csv"), row.names = FALSE)
utils::write.csv(b_cell_records, file.path(out_dir, "b_cell_search_results.csv"), row.names = FALSE)
utils::write.csv(cd8_records, file.path(out_dir, "cd8_search_results.csv"), row.names = FALSE)
utils::write.csv(matrix_rank_scores, file.path(out_dir, "matrix_rank_scores.csv"), row.names = TRUE)
utils::write.csv(matrix_mean_scores, file.path(out_dir, "matrix_mean_scores.csv"), row.names = TRUE)

cat("ImmuneSigR local demo completed successfully.\n")
cat("Output directory:", out_dir, "\n")
