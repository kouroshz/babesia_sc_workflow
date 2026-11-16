# Babesia scRNA-seq: original WT 0 h versus seven days at 4 C.
# Open Babesia_scRNA_Tutorial.Rproj. Run one numbered section at a time.
# Everything starts from count matrices; no saved Seurat RDS is an input.
# The supplied WT0 phase-label CSV is an explicit annotated reference input.

# 0. Obtain the data in the RStudio Terminal, NOT the R Console -----------------
# Kourosh's approved login route:
# ssh -J kouroshz@markov.math.umb.edu kourosh.zarringhalam@chimera.umb.edu
# Molly and Aryan should substitute their own authorized institutional accounts.
#
# Run these two rsync commands from this repository's root. They copy only the
# filtered/raw H5 files and 10x reports. No BAM or FASTQ is needed for this lesson.
#
# rsync -av --partial -e "ssh -J kouroshz@markov.math.umb.edu" \
#   --include='filtered_feature_bc_matrix.h5' --include='raw_feature_bc_matrix.h5' \
#   --include='web_summary.html' --include='metrics_summary.csv' --exclude='*' \
#   kourosh.zarringhalam@chimera.umb.edu:/mathbiospace/data01/mathbio_lab/parasites/acdc_data/bdivergens/scRNAseq/2022-07-02_bd_cold_shock_10x/bd0hr/result_cr7/ref_PiroplasmaDB-68_Bdivergens1802A/count/outs/ \
#   input/WT0/
#
# rsync -av --partial -e "ssh -J kouroshz@markov.math.umb.edu" \
#   --include='filtered_feature_bc_matrix.h5' --include='raw_feature_bc_matrix.h5' \
#   --include='web_summary.html' --include='metrics_summary.csv' --exclude='*' \
#   kourosh.zarringhalam@chimera.umb.edu:/mathbiospace/data01/mathbio_lab/parasites/acdc_data/bdivergens/scRNAseq/2022-07-02_bd_cold_shock_10x/bd7d/result_cr7/ref_PiroplasmaDB-68_Bdivergens1802A/count/outs/ \
#   input/WT7/
#
# These lines are comments: sourcing this R script never connects to Chimera.
# Read web_summary.html for each sample before interpreting our extra filtering.
# Check the downloaded H5 hashes against input/data_files.csv in the Terminal:
# shasum -a 256 input/WT0/*.h5 input/WT7/*.h5

# 1. Start the session and inspect the choices ---------------------------------
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(limma)                    # Baseline Seurat Wilcoxon backend; pinned in renv.
source("tutorial_functions.R")
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3) # Allow the two-library integration objects.
set.seed(100)

SAVE_OUTPUTS <- TRUE              # FALSE: view figures in RStudio without saving.
WARM_MIN_GENES <- 200             # Strict: genes > 200.
WARM_MIN_UMIS <- 250              # Inclusive: UMIs >= 250.
COLD_MIN_GENES <- 75              # Lower gate preserves low-RNA cold cells.
COLD_MIN_UMIS <- 100
EXPECTED_DOUBLET_RATE <- 0.05      # Assumption, not a measured doublet rate.
WT0_PK <- 0.005                   # Frozen choices from each sample's BCmvn sweep.
WT7_PK <- 0.01                    # Do not copy these values to unrelated samples.

for (folder in c("figures", "tables", "objects")) {
  dir.create(file.path("output", folder), recursive = TRUE, showWarnings = FALSE)
}
annotation <- read.csv("input/gene_annotation.csv", stringsAsFactors = FALSE)
reference_labels <- read.csv("input/WT0_phase_labels.csv", stringsAsFactors = FALSE)
markers <- read.csv("input/phase_markers.csv", stringsAsFactors = FALSE)

# 2. Read integer UMI counts and inspect what 10x called ------------------------
wt0_counts <- read_counts("input/WT0/filtered_feature_bc_matrix.h5")
wt7_counts <- read_counts("input/WT7/filtered_feature_bc_matrix.h5")
dim(wt0_counts)                   # Gene rows by Cell Ranger-called barcodes.
dim(wt7_counts)
# Both matrices should use the same reference gene list. Never silently mix
# matrices aligned to different references or with reordered gene rows.
stopifnot(identical(rownames(wt0_counts), rownames(wt7_counts)))
wt0_raw <- raw_barcode_summary("input/WT0/raw_feature_bc_matrix.h5")
wt7_raw <- raw_barcode_summary("input/WT7/raw_feature_bc_matrix.h5")
wt0_raw                          # Raw barcode columns are NOT known cells.
wt7_raw

# Prefix cell IDs so the same 10x barcode can occur in both libraries safely.
colnames(wt0_counts) <- paste0("WT0_", colnames(wt0_counts))
colnames(wt7_counts) <- paste0("WT7_", colnames(wt7_counts))
wt0_all <- CreateSeuratObject(wt0_counts, project = "WT0", min.cells = 0, min.features = 0)
wt7_all <- CreateSeuratObject(wt7_counts, project = "WT7", min.cells = 0, min.features = 0)
wt0_all$input_genes <- Matrix::colSums(wt0_counts > 0)
wt0_all$input_UMIs <- Matrix::colSums(wt0_counts)
wt7_all$input_genes <- Matrix::colSums(wt7_counts > 0)
wt7_all$input_UMIs <- Matrix::colSums(wt7_counts)

# Seurat replaces underscores in feature IDs with dashes. Gene tables keep the
# original Bdiv_ IDs; rna_deg() translates the prefix back when joining names.
head(rownames(wt0_all))
head(wt0_all[[]])

# 3. Review RNA complexity BEFORE applying our gates ----------------------------
qc0 <- ggplot(wt0_all[[]], aes(input_UMIs, input_genes)) +
  geom_point(size = 0.5, alpha = 0.45, color = "#237E83") +
  geom_vline(xintercept = WARM_MIN_UMIS, linetype = "dashed") +
  geom_hline(yintercept = WARM_MIN_GENES, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() + theme_classic(base_size = 17) +
  labs(title = "WT 0 h: all 10x-called cells", x = "Observed UMIs (10x-called cells)", y = "Detected genes")
qc7 <- ggplot(wt7_all[[]], aes(input_UMIs, input_genes)) +
  geom_point(size = 0.5, alpha = 0.45, color = "#416DA5") +
  geom_vline(xintercept = COLD_MIN_UMIS, linetype = "dashed") +
  geom_hline(yintercept = COLD_MIN_GENES, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() + theme_classic(base_size = 17) +
  labs(title = "WT 7 d cold: all 10x-called cells", x = "Observed UMIs (10x-called cells)", y = "Detected genes")
show_figure(qc0 + qc7, "output/figures/01_RNA_cutoffs.png", SAVE_OUTPUTS, width = 14)

# These are the fixed CIQ baseline gates, not universal Seurat recommendations.
# No additional EmptyDrops, clipping-score, upper-UMI or mitochondrial proxy gate.
wt0_pass <- wt0_all$input_genes > WARM_MIN_GENES & wt0_all$input_UMIs >= WARM_MIN_UMIS
wt7_pass <- wt7_all$input_genes > COLD_MIN_GENES & wt7_all$input_UMIs >= COLD_MIN_UMIS
table(wt0_pass)
table(wt7_pass)

# Cell gates use ALL genes. The >=3-cell GENE rule below affects embedding only.
wt0_gated <- CreateSeuratObject(wt0_counts[, wt0_pass], project = "WT0", min.cells = 3)
wt7_gated <- CreateSeuratObject(wt7_counts[, wt7_pass], project = "WT7", min.cells = 3)

# 4. WT0 normalization, variable genes, scaling and PCA: visible Seurat calls ---
# NormalizeData divides by each cell's total counts, multiplies by 10,000 and
# takes log1p. It changes relative expression, not the original UMI counts.
wt0_gated <- NormalizeData(wt0_gated, scale.factor = 10000, verbose = FALSE)
wt0_gated <- FindVariableFeatures(wt0_gated, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
wt0_gated <- ScaleData(wt0_gated, features = VariableFeatures(wt0_gated), verbose = FALSE)
wt0_gated <- RunPCA(wt0_gated, npcs = 30, seed.use = 100, verbose = FALSE)
show_figure(ElbowPlot(wt0_gated, ndims = 30), "output/figures/02_WT0_PCA_elbow.png", SAVE_OUTPUTS)
show_figure(DimPlot(wt0_gated, reduction = "pca"), "output/figures/03_WT0_PCA.png", SAVE_OUTPUTS)

# 5. WT0 neighbors, clusters and UMAP ------------------------------------------
# Use all 30 fitted PCs for the independent embedding, matching the baseline.
wt0_gated <- FindNeighbors(wt0_gated, dims = 1:30, verbose = FALSE)
wt0_gated <- FindClusters(wt0_gated, resolution = 0.4, random.seed = 100, verbose = FALSE)
wt0_gated <- RunUMAP(wt0_gated, dims = 1:30, seed.use = 100, verbose = FALSE)
show_figure(DimPlot(wt0_gated, label = TRUE) + ggtitle("WT0 before doublet removal"),
  "output/figures/04_WT0_initial_UMAP.png", SAVE_OUTPUTS)

# For WT7, fit_embedding() repeats the same seven preprocessing commands above.
# Open that short function in tutorial_functions.R to inspect every call.
wt7_gated <- fit_embedding(wt7_gated)

# 6. DoubletFinder separately in each library ----------------------------------
# pK was chosen previously by BCmvn, using all real cells and pN=.1/.2/.3.
# The actual final DoubletFinder call uses pN=0.25 to generate artificial doublets.
# nExp is 5% of cells, adjusted by the sample's cluster-derived homotypic rate.
# For a new sample, optionally run choose_doublet_pk(object) and review metrics.
wt0_called <- call_doublets(wt0_gated, pK = WT0_PK, expected_rate = EXPECTED_DOUBLET_RATE)
wt7_called <- call_doublets(wt7_gated, pK = WT7_PK, expected_rate = EXPECTED_DOUBLET_RATE)
table(wt0_called$doublet_call)
table(wt7_called$doublet_call)
if (SAVE_OUTPUTS) {
  write.csv(data.frame(cell = colnames(wt0_called), wt0_called[[]]),
    "output/tables/WT0_DoubletFinder_all_cells.csv", row.names = FALSE)
  write.csv(data.frame(cell = colnames(wt7_called), wt7_called[[]]),
    "output/tables/WT7_DoubletFinder_all_cells.csv", row.names = FALSE)
}
p0 <- DimPlot(wt0_called, group.by = "doublet_call",
  cols = c(Doublet = "#C64D5E", Singlet = "#237E83")) + ggtitle("WT0 doublet calls")
p7 <- DimPlot(wt7_called, group.by = "doublet_call",
  cols = c(Doublet = "#C64D5E", Singlet = "#237E83")) + ggtitle("WT7 cold doublet calls")
show_figure(p0 + p7,
  "output/figures/05_DoubletFinder_calls.png", SAVE_OUTPUTS, width = 14)

# Rebuild from original counts and refit normalization, PCA and UMAP on singlets.
wt0 <- rebuild_singlets(wt0_counts, wt0_called, "WT0")
wt7 <- rebuild_singlets(wt7_counts, wt7_called, "WT7")
qc_table <- rbind(cell_summary(wt0_counts, wt0_pass, wt0_called, wt0, wt0_raw),
  cell_summary(wt7_counts, wt7_pass, wt7_called, wt7, wt7_raw))
qc_table                         # Expected baseline singlets: 9,532 and 5,571.
funnel <- qc_table |> select(sample, CellRanger_cells, RNA_pass, final_singlets) |>
  pivot_longer(-sample, names_to = "stage", values_to = "cells") |>
  mutate(stage = factor(stage, c("CellRanger_cells", "RNA_pass", "final_singlets")))
p <- ggplot(funnel, aes(stage, cells, fill = stage)) + geom_col(width = 0.65) +
  geom_text(aes(label = format(cells, big.mark = ",")), vjust = -0.35, size = 5.5) +
  facet_wrap(~sample) + scale_fill_manual(values = c("#8796AC", "#D1AA45", "#237E83")) +
  scale_x_discrete(labels = c("10x called", "RNA gate", "Final singlets")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 17) + theme(legend.position = "none") +
  labs(title = "Cell counts at each filtering step", x = NULL, y = "Cells")
show_figure(p, "output/figures/06_filtering_counts.png", SAVE_OUTPUTS, width = 13)
p0 <- DimPlot(wt0, label = TRUE) + ggtitle("WT0 singlets")
p7 <- DimPlot(wt7, label = TRUE) + ggtitle("WT7 cold singlets")
show_figure(p0 + p7,
  "output/figures/07_independent_singlet_UMAPs.png", SAVE_OUTPUTS, width = 14)
if (SAVE_OUTPUTS) {
  write.csv(qc_table, "output/tables/cell_and_gene_QC.csv", row.names = FALSE)
  saveRDS(wt0, "output/objects/WT0_singlets.rds")
  saveRDS(wt7, "output/objects/WT7_singlets.rds")
}

# 7. Inspect the annotated WT0 reference and TRANSFER phase labels --------------
# The provided CSV openly supplies the fixed project WT0 RNA-phase annotations.
# It is not raw experimental truth and does not rebuild the research clock.
# G/SM/MC/C are Babesia RNA states, not mammalian G1/S/G2/M measurements.
wt0$phase <- reference_labels$phase[match(wt0$barcode, reference_labels$barcode)]
table(wt0$phase, useNA = "ifany")
reference <- subset(wt0, cells = colnames(wt0)[!is.na(wt0$phase)])
# Changed QC may admit unannotated WT0 cells. Keep them visible as Unassigned;
# do not pretend they have known labels or silently remove them from the sample.
wt0$phase_display <- ifelse(is.na(wt0$phase), "Unassigned", wt0$phase)
anchors <- FindTransferAnchors(reference = reference, query = wt7,
  normalization.method = "LogNormalize", reference.reduction = "pca",
  reduction = "pcaproject", dims = 1:20,
  features = intersect(VariableFeatures(reference), rownames(wt7)), verbose = FALSE)
predictions <- TransferData(anchorset = anchors, refdata = reference$phase,
  dims = 1:20, verbose = FALSE)
wt7$phase <- as.character(predictions[colnames(wt7), "predicted.id"])
wt7$phase_score <- predictions[colnames(wt7), "prediction.score.max"]
wt7$phase_display <- ifelse(wt7$phase_score >= 0.50, wt7$phase, "Uncertain")
summary(wt7$phase_score)

# Native coordinates do not change when metadata are added.
phase_colors <- c(G = "#416DA5", SM = "#239481", MC = "#D1AA45", C = "#C64D5E",
  Uncertain = "#ADB5BA", Unassigned = "#ADB5BA")
wt0$phase_display <- factor(wt0$phase_display, levels = names(phase_colors))
wt7$phase_display <- factor(wt7$phase_display, levels = names(phase_colors))
p0 <- DimPlot(wt0, group.by = "phase_display", cols = phase_colors) +
  ggtitle("WT0: supplied reference labels")
p7 <- DimPlot(wt7, group.by = "phase_display", cols = phase_colors) +
  ggtitle("WT7: newly transferred labels")
show_figure(p0 + p7,
  "output/figures/08_phase_transfer_native_UMAPs.png", SAVE_OUTPUTS, width = 14)
show_figure(FeaturePlot(wt7, features = "phase_score") + ggtitle("WT7 transfer confidence"),
  "output/figures/09_phase_transfer_confidence.png", SAVE_OUTPUTS)

# 8. Joint integration is a SEPARATE view, not an expression matrix for DEG ----
# It can remove real condition differences. Always retain the independent UMAPs.
objects <- list(WT0 = wt0, WT7 = wt7)
features <- SelectIntegrationFeatures(objects, nfeatures = 3000)
features <- intersect(features, intersect(rownames(wt0), rownames(wt7)))
objects <- lapply(objects, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  RunPCA(x, features = features, npcs = 30, seed.use = 100, verbose = FALSE)
})
integration_anchors <- FindIntegrationAnchors(object.list = objects,
  anchor.features = features, reduction = "rpca", reference = 1,
  dims = 1:20, k.filter = NA, verbose = TRUE)
integrated <- IntegrateData(integration_anchors, dims = 1:20, verbose = TRUE)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, seed.use = 100, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:20, seed.use = 100, verbose = FALSE)
show_figure(DimPlot(integrated, group.by = "sample_id"),
  "output/figures/10_integrated_by_sample.png", SAVE_OUTPUTS)
show_figure(DimPlot(integrated, group.by = "phase_display", split.by = "sample_id",
    cols = phase_colors, ncol = 2),
  "output/figures/11_integrated_split_by_sample.png", SAVE_OUTPUTS, width = 14)

# 9. Merge full-gene RNA objects for DEG, without integration -------------------
# Keep only final singlets, but restore ALL input gene rows for each sample.
# The >=3-cell gene rule was for embedding, not the DEG gene universe.
wt0_metadata <- wt0[[]]
wt7_metadata <- wt7[[]]
# Preserve the embedding-assay QC summaries with explicit names. The new RNA
# objects recompute nCount_RNA and nFeature_RNA using the complete gene set.
names(wt0_metadata)[names(wt0_metadata) == "nCount_RNA"] <- "embedding_nCount_RNA"
names(wt0_metadata)[names(wt0_metadata) == "nFeature_RNA"] <- "embedding_nFeature_RNA"
names(wt7_metadata)[names(wt7_metadata) == "nCount_RNA"] <- "embedding_nCount_RNA"
names(wt7_metadata)[names(wt7_metadata) == "nFeature_RNA"] <- "embedding_nFeature_RNA"
rna0 <- CreateSeuratObject(wt0_counts[, colnames(wt0), drop = FALSE],
  project = "WT0", min.cells = 0, min.features = 0, meta.data = wt0_metadata)
rna7 <- CreateSeuratObject(wt7_counts[, colnames(wt7), drop = FALSE],
  project = "WT7", min.cells = 0, min.features = 0, meta.data = wt7_metadata)

# merge() combines cells and their metadata, without batch correction or summing
# samples. Existing WT0_/WT7_ prefixes already give every cell a unique ID.
rna <- merge(x = rna0, y = rna7, merge.data = FALSE, project = "WT0_WT7")
Layers(rna[["RNA"]])             # Separate sample count layers before joining.
rna <- JoinLayers(rna, assay = "RNA")
Layers(rna[["RNA"]])             # One counts layer; sample_id remains metadata.
# Keep the same G, SM, MC, C display order after merging sample metadata.
rna$phase_display <- factor(rna$phase_display, levels = names(phase_colors))
rna <- NormalizeData(rna, scale.factor = 10000, verbose = FALSE)
DefaultAssay(rna) <- "RNA"
dim(rna)                         # Full gene rows, not just embedding genes.

# Marker inspection is useful but not independent validation of the reference.
marker_groups <- split(markers$gene_id, factor(markers$phase, c("G", "SM", "MC", "C")))
marker_groups <- lapply(marker_groups, function(ids) sub("^Bdiv_", "Bdiv-", ids))
rna$sample_phase <- factor(paste(rna$sample_id, rna$phase, sep = ": "),
  levels = c("WT0: G", "WT0: SM", "WT0: MC", "WT0: C",
    "WT7: G", "WT7: SM", "WT7: MC", "WT7: C"))
show_figure(DotPlot(rna, features = marker_groups, group.by = "sample_phase") +
    RotatedAxis() + ggtitle("Published Babesia phase markers"),
  "output/figures/12_phase_marker_check.png", SAVE_OUTPUTS, width = 15)

# Count both forced phase assignments and confidence-aware display labels.
# The denominator is always all retained cells in that sample.
phase_counts <- rna[[]] |>
  count(sample_id, phase, name = "cells") |>
  group_by(sample_id) |> mutate(total_cells = sum(cells), fraction = cells / total_cells) |> ungroup()
phase_display_counts <- rna[[]] |>
  count(sample_id, phase_display, name = "cells") |>
  group_by(sample_id) |> mutate(total_cells = sum(cells), fraction = cells / total_cells) |> ungroup()
p <- ggplot(phase_display_counts, aes(sample_id, fraction, fill = phase_display)) +
  geom_col(width = 0.65) + geom_text(aes(label = ifelse(fraction >= 0.03,
    paste0(round(100 * fraction, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 5, color = "white") +
  scale_fill_manual(values = phase_colors) + theme_classic(base_size = 17) +
  labs(title = "Phase proportions in retained cells", y = "Fraction of all retained cells",
    x = NULL, fill = "RNA phase", subtitle = "Gray marks uncertain or unassigned labels, not excluded cells")
show_figure(p, "output/figures/13_phase_proportions.png", SAVE_OUTPUTS)

# 10. Whole-sample DEG: WT7 cold minus WT0 --------------------------------------
# This means all-cell RNA comparison, NOT a replicate-aware bulk-RNA experiment.
# There is one biological sample per condition; cell-level P-values do not
# provide biological replication. This is an exploratory teaching comparison.
# The result mixes condition, phase composition, RNA complexity and selection
# by different warm/cold gates; it does not isolate absolute cold repression.
whole_deg <- rna_deg(rna,
  case_cells = colnames(rna)[rna$sample_id == "WT7"],
  control_cells = colnames(rna)[rna$sample_id == "WT0"], annotation = annotation)
table(whole_deg$direction)
head(whole_deg[whole_deg$called, c("gene_id", "description", "log2FC", "Bonferroni_P")], 20)
whole_deg[whole_deg$prespecified_candidate,
  c("gene_id", "description", "log2FC", "Bonferroni_P", "direction")]
show_figure(volcano_plot(whole_deg), "output/figures/14_whole_sample_DEG.png", SAVE_OUTPUTS)

# 11. Phase-matched DEG: compare cold and warm WITHIN each transferred phase -----
# The exact same RNA test and cutoffs are used in G, SM, MC and C.
# Tests use phase assignments, not integrated expression or UMAP position.
# Low transfer scores are displayed above, not silently filtered out here.
# WT0 reference labels and WT7 predictions can have different boundary errors.
# That uncertainty can contribute to apparent phase-matched expression effects.
phase_results <- list()
coverage <- list()
for (phase_name in c("G", "SM", "MC", "C")) {
  case <- colnames(rna)[which(rna$sample_id == "WT7" & rna$phase == phase_name)]
  control <- colnames(rna)[which(rna$sample_id == "WT0" & rna$phase == phase_name)]
  eligible <- min(length(case), length(control)) >= 20
  coverage[[phase_name]] <- data.frame(phase = phase_name,
    WT7_cells = length(case), WT0_cells = length(control), eligible)
  if (!eligible) next             # Not tested is NOT zero significant genes.
  result <- rna_deg(rna, case, control, annotation)
  result$phase <- phase_name
  phase_results[[phase_name]] <- result
}
phase_deg <- bind_rows(phase_results)
coverage <- bind_rows(coverage)
coverage
deg_counts <- phase_deg |> filter(called) |> count(phase, direction, name = "genes") |>
  complete(phase = c("G", "SM", "MC", "C"), direction = c("up", "down"), fill = list(genes = 0L)) |>
  left_join(coverage, by = "phase") |> mutate(genes = ifelse(eligible, genes, NA_integer_),
    phase = factor(phase, levels = c("G", "SM", "MC", "C")),
    direction = factor(direction, levels = c("up", "down")))
p <- ggplot(deg_counts, aes(phase, genes, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7, na.rm = TRUE) +
  geom_text(aes(label = genes), position = position_dodge(0.8), vjust = -0.35,
    size = 5.5, na.rm = TRUE) +
  geom_text(data = subset(coverage, !eligible), aes(phase, 0, label = "NT"),
    inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c(up = "#C64D5E", down = "#247FA8")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + theme_classic(base_size = 17) +
  labs(title = "WT7 versus WT0: phase-matched DEGs", x = "Babesia RNA phase (WT0 reference / WT7 predicted)",
    y = "Called genes", fill = NULL,
    subtitle = "Bonferroni P < 0.05, >1.5-fold; at least 20 cells per group")
show_figure(p, "output/figures/15_phase_matched_DEG_counts.png", SAVE_OUTPUTS)

# 12. Save the foundations for your own next plot or analysis -------------------
# These RDS files are OUTPUTS of the code above, never undocumented inputs.
if (SAVE_OUTPUTS) {
  write.csv(whole_deg, "output/tables/WT7_vs_WT0_all_genes.csv", row.names = FALSE)
  write.csv(subset(whole_deg, called), "output/tables/WT7_vs_WT0_DEGs.csv", row.names = FALSE)
  write.csv(phase_deg, "output/tables/WT7_vs_WT0_phase_all_genes.csv", row.names = FALSE)
  write.csv(subset(phase_deg, called), "output/tables/WT7_vs_WT0_phase_DEGs.csv", row.names = FALSE)
  write.csv(coverage, "output/tables/phase_test_cell_counts.csv", row.names = FALSE)
  write.csv(phase_counts, "output/tables/phase_counts.csv", row.names = FALSE)
  write.csv(phase_display_counts, "output/tables/phase_counts_with_uncertainty.csv", row.names = FALSE)
  saveRDS(wt0, "output/objects/WT0_phase_reference.rds")
  saveRDS(wt7, "output/objects/WT7_phase_transferred.rds")
  saveRDS(integrated, "output/objects/WT0_WT7_integrated.rds")
  saveRDS(list(whole_sample = whole_deg, by_phase = phase_deg,
    phase_counts = phase_counts), "output/objects/DEG_results.rds")
  capture.output(sessionInfo(), file = "output/sessionInfo.txt")
}

# Next session: load ONLY what your new question needs. Examples:
# wt7 <- readRDS("output/objects/WT7_phase_transferred.rds")
# FeaturePlot(wt7, features = "Bdiv-031680c")
# deg <- readRDS("output/objects/DEG_results.rds")$whole_sample
# View(subset(deg, called))
#
# New samples: keep WT0 as a reference, QC and fit each new sample independently,
# select its own DoubletFinder pK, then reuse the label-transfer block. Do not
# assume that a cyclic UMAP, a transferred label or an abundance shift proves
# viability, measured cycle arrest, total RNA production, or cold survival.
