# Supporting functions for run_tutorial.R.
# Scientific choices stay visible in the run script. No function downloads data,
# installs packages, reads an old Seurat object, or changes the CIQ project.

read_counts <- function(file) {
  counts <- Seurat::Read10X_h5(file, use.names = FALSE, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  stopifnot(inherits(counts, "sparseMatrix"), !anyDuplicated(rownames(counts)),
    !anyDuplicated(colnames(counts)))
  counts
}

raw_barcode_summary <- function(file) {
  # Read the H5 dimensions and sparse-column pointers, not a dense raw matrix.
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  shape <- h5[["matrix/shape"]][]
  pointers <- h5[["matrix/indptr"]][]
  data.frame(raw_gene_rows = shape[1], raw_barcode_columns = shape[2],
    raw_barcodes_with_counts = sum(diff(pointers) > 0))
}

fit_embedding <- function(object, seed = 100, variable_genes = 3000,
                          pcs = 30, resolution = 0.4) {
  # These are the same seven commands shown explicitly for WT0 in the tutorial.
  set.seed(seed)
  object <- Seurat::NormalizeData(object, scale.factor = 10000, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst",
    nfeatures = min(variable_genes, nrow(object) - 1L), verbose = FALSE)
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  pcs <- min(pcs, length(Seurat::VariableFeatures(object)) - 1L, ncol(object) - 1L)
  object <- Seurat::RunPCA(object, npcs = pcs, seed.use = seed, verbose = FALSE)
  object <- Seurat::FindNeighbors(object, dims = seq_len(pcs), verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = resolution,
    random.seed = seed, verbose = FALSE)
  Seurat::RunUMAP(object, dims = seq_len(pcs), seed.use = seed, verbose = FALSE)
}

call_doublets <- function(object, pK, expected_rate = 0.05, seed = 100) {
  # Five percent is an assumed loading rate, not an observed doublet frequency.
  # Homotypic adjustment discounts cells resembling the same expression cluster.
  if (any(grepl("^DF.classifications_", names(object[[]])))) {
    stop("Use the RNA-gated object, not an object already carrying doublet calls.")
  }
  homotypic <- DoubletFinder::modelHomotypic(as.character(Seurat::Idents(object)))
  expected_raw <- round(expected_rate * ncol(object))
  expected_adjusted <- max(1L, round(expected_raw * (1 - homotypic)))
  set.seed(seed)
  # DoubletFinder repeats preprocessing using this object's recorded commands.
  object <- DoubletFinder::doubletFinder(object,
    PCs = seq_len(ncol(Seurat::Embeddings(object, "pca"))), pN = 0.25, pK = pK,
    nExp = expected_adjusted, reuse.pANN = NULL, sct = FALSE)
  scores <- grep("^pANN_", names(object[[]]), value = TRUE)
  calls <- grep("^DF.classifications_", names(object[[]]), value = TRUE)
  stopifnot(length(scores) == 1L, length(calls) == 1L)
  object$doublet_score <- object[[]][[scores]]
  object$doublet_call <- object[[]][[calls]]
  object@misc$doublet_settings <- data.frame(pK, pN = 0.25, expected_rate,
    homotypic_fraction = homotypic, expected_raw, expected_adjusted)
  print(object@misc$doublet_settings)
  object
}

rebuild_singlets <- function(counts, called_object, sample_id, seed = 100) {
  # Rebuild from counts. Do not merely delete points from the old UMAP.
  keep <- colnames(called_object)[called_object$doublet_call == "Singlet"]
  object <- Seurat::CreateSeuratObject(counts[, keep, drop = FALSE],
    project = sample_id, min.cells = 3, min.features = 0)
  object$sample_id <- sample_id
  object$input_UMIs <- Matrix::colSums(counts[, keep, drop = FALSE])
  object$input_genes <- Matrix::colSums(counts[, keep, drop = FALSE] > 0)
  object$doublet_score <- called_object$doublet_score[match(keep, colnames(called_object))]
  object$barcode <- sub(paste0("^", sample_id, "_"), "", keep)
  fit_embedding(object, seed = seed)
}

cell_summary <- function(counts, pass_rna, called_object, singlets, raw_summary) {
  retained <- counts[, colnames(singlets), drop = FALSE]
  data.frame(sample = unique(singlets$sample_id), raw_summary,
    CellRanger_cells = ncol(counts), RNA_pass = sum(pass_rna),
    predicted_doublets = sum(called_object$doublet_call == "Doublet"),
    final_singlets = ncol(singlets), input_gene_rows = nrow(counts),
    genes_detected_in_final_cells = sum(Matrix::rowSums(retained > 0) > 0),
    embedding_genes_detected_in_3_cells = nrow(singlets),
    total_UMIs = sum(retained), median_UMIs = median(Matrix::colSums(retained)),
    median_genes = median(Matrix::colSums(retained > 0)))
}

rna_deg <- function(object, case_cells, control_cells, annotation, minimum_cells = 20) {
  # RNA only. Integration changes coordinates, never the values tested here.
  stopifnot(Seurat::DefaultAssay(object) == "RNA")
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Run setup.R to restore the pinned limma Wilcoxon backend.")
  }
  if (min(length(case_cells), length(control_cells)) < minimum_cells) {
    stop("Need at least ", minimum_cells, " cells in each comparison group.")
  }
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  cells <- c(case_cells, control_cells)
  normalized <- counts[, cells, drop = FALSE] %*%
    Matrix::Diagonal(x = 10000 / Matrix::colSums(counts[, cells, drop = FALSE]))
  dimnames(normalized) <- list(rownames(counts), cells)
  a <- Matrix::rowMeans(normalized[, case_cells, drop = FALSE])
  b <- Matrix::rowMeans(normalized[, control_cells, drop = FALSE])
  effects <- data.frame(feature = rownames(counts),
    gene_id = sub("^Bdiv-", "Bdiv_", rownames(counts)),
    log2FC = log2((a + 1) / (b + 1)),
    pct_case = Matrix::rowMeans(counts[, case_cells, drop = FALSE] > 0),
    pct_control = Matrix::rowMeans(counts[, control_cells, drop = FALSE] > 0))

  # Force the same limma Wilcoxon backend used by CIQ when presto was absent.
  # There is no FC pre-screen. min.pct decides which genes receive a test.
  tests <- Seurat::FindMarkers(object[["RNA"]], cells.1 = case_cells,
    cells.2 = control_cells, test.use = "wilcox_limma", slot = "data",
    min.pct = 0.10, logfc.threshold = 0, max.cells.per.ident = Inf, verbose = FALSE)
  tests <- data.frame(feature = rownames(tests), P_value = tests$p_val,
    Bonferroni_P = tests$p_val_adj)
  result <- dplyr::left_join(effects, tests, by = "feature")
  result <- dplyr::left_join(result, annotation, by = "gene_id")
  result$tested <- !is.na(result$P_value)
  result$called <- result$tested & result$Bonferroni_P < 0.05 & abs(result$log2FC) > log2(1.5)
  result$direction <- ifelse(!result$tested, "not_tested", ifelse(!result$called, "ns",
    ifelse(result$log2FC > 0, "up", "down")))
  result$case_cells <- length(case_cells)
  result$control_cells <- length(control_cells)
  result$prespecified_candidate <- result$gene_id %in% c("Bdiv_023810", "Bdiv_031680c", "Bdiv_038820")
  dplyr::arrange(result, factor(direction, c("up", "down", "ns", "not_tested")),
    dplyr::desc(abs(log2FC)), Bonferroni_P, gene_id)
}

volcano_plot <- function(result) {
  # Named top genes are labels for readability, not a new selection criterion.
  d <- result[result$tested, ]
  d$evidence <- -log10(pmax(d$Bonferroni_P, 1e-300))
  named <- d[d$called & !grepl("hypothetical|uncharacterized", d$description, ignore.case = TRUE), ]
  named <- named[order(named$Bonferroni_P, -abs(named$log2FC)), ]
  label_ids <- unique(c(head(named$gene_id, 6), d$gene_id[d$prespecified_candidate]))
  labels <- d[d$gene_id %in% label_ids, ]
  candidate_names <- c(Bdiv_023810 = "Rb KD target", Bdiv_031680c = "SAMS", Bdiv_038820 = "EAPP-like")
  labels$label <- ifelse(labels$prespecified_candidate,
    paste0(candidate_names[labels$gene_id], "\n", labels$gene_id), labels$gene_id)
  ggplot2::ggplot(d, ggplot2::aes(log2FC, evidence, color = direction)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * log2(1.5), linetype = "dashed", color = "gray45") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray45") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label),
      size = 4.3, seed = 100, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(up = "#C64D5E", down = "#247FA8", ns = "#C7CDD1")) +
    ggplot2::labs(x = "WT 7 d cold - WT 0 h log2FC", y = "-log10(Bonferroni-adjusted P)",
      title = "Whole-sample cold versus warm expression", color = NULL,
      subtitle = "Cell-level Wilcoxon tests, one biological sample per condition") +
    ggplot2::theme_classic(base_size = 17)
}

show_figure <- function(plot, file, save = TRUE, width = 11, height = 7) {
  # Always show the figure in RStudio; saving is a separate visible choice.
  print(plot)
  if (save) {
    ggplot2::ggsave(file, plot, width = width, height = height, dpi = 220, bg = "white")
    ggplot2::ggsave(sub("\\.png$", ".pdf", file), plot,
      width = width, height = height, bg = "white")
  }
  invisible(plot)
}

# Optional for NEW samples. The tutorial uses reviewed pK values, not this sweep.
choose_doublet_pk <- function(object, seed = 100) {
  # Same pinned DoubletFinder worker as CIQ, using ALL real cells and pN=.1/.2/.3.
  # Official paramSweep can cap real cells at 10,000; this explicit worker avoids it.
  if (as.character(utils::packageVersion("DoubletFinder")) != "2.0.6") {
    stop("This full-cell sweep is written for the pinned DoubletFinder 2.0.6.")
  }
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  real_cells <- colnames(object)
  pN <- c(0.10, 0.20, 0.30)
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.30, by = 0.01))
  minimum_artificial <- round(length(real_cells) / (1 - min(pN)) - length(real_cells))
  pK <- pK[round(pK * minimum_artificial) >= 1]
  worker <- getFromNamespace("parallel_paramSweep", "DoubletFinder")
  sweep <- list()
  for (i in seq_along(pN)) {
    set.seed(seed + i)
    z <- worker(n = i, n.real.cells = length(real_cells), real.cells = real_cells,
      pK = pK, pN = pN, data = counts, orig.commands = object@commands,
      PCs = seq_len(ncol(Seurat::Embeddings(object, "pca"))), sct = FALSE)
    names(z) <- paste("pN", pN[i], "pK", pK, sep = "_")
    sweep <- c(sweep, z)
  }
  stopifnot(all(vapply(sweep, nrow, integer(1)) == length(real_cells)))
  stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  metrics <- DoubletFinder::find.pK(stats)
  metrics$pK <- as.numeric(as.character(metrics$pK))
  metrics <- metrics[is.finite(metrics$BCmetric), ]
  metrics <- metrics[order(-metrics$BCmetric, metrics$pK), ]
  list(selected_pK = metrics$pK[1], metrics = metrics)
}
