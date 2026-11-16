# Validation

Checked on **16 November 2026** on Kourosh's macOS computer. The main CIQ
analysis and its presentations were not modified.

## What Was Run

- Complete fresh processing from both original Cell Ranger filtered matrices.
  Raw H5s were read only for barcode accounting. No processed Seurat object was
  a starting input. The supplied WT0 annotation CSV is an explicit reference
  input, as explained in the README.
- Independent RNA QC, DoubletFinder, singlet rebuilding, normalization, PCA,
  clustering, native UMAPs, phase transfer, RPCA integration, marker checks,
  phase counts, whole-sample DEG and four phase-matched DEG comparisons.
- A separate project library restored from `renv.lock`, followed by another
  complete count-based run using that isolated library. `renv::status()` reports
  no inconsistencies. Most dependencies were linked from the existing package
  cache. This is not a test on a second operating system.
- The explicit `wilcox_limma` backend was checked against CIQ's original
  `wilcox` implementation with limma available and presto absent. Every
  whole-sample and phase test is numerically identical.
- All 15 numbered figure families were generated as PNG and PDF and inspected.

## Baseline Agreement

| Sample | 10x-called cells | RNA-gate pass | Predicted doublets | Final singlets |
|---|---:|---:|---:|---:|
| WT0 | 9,974 | 9,949 | 417 | 9,532 |
| WT7 | 5,812 | 5,810 | 239 | 5,571 |

The retained **barcode sets**, not only their counts, exactly match the CIQ
baseline. Final detected-gene and UMI summaries also match.

The whole-sample test restores all **4,129** input gene rows and reproduces
**79 up / 224 down** calls. Compared with the existing CIQ whole-sample table,
the maximum log2FC difference is approximately `2.22e-16` and the maximum
Bonferroni-P difference approximately `5.55e-17`, both numerical roundoff.
Repeated runs on this computer produced identical whole and phase DEG tables.

## Standard Seurat Merge Update

On **16 November 2026**, section 9 was changed from manual count-matrix binding
to two full-gene Seurat objects followed by `merge(merge.data=FALSE)`,
`JoinLayers()` and `NormalizeData()`. No filtering, phase assignment or DEG
settings changed.

The changed section and all five downstream DEG comparisons were rerun from
the original filtered count matrices, using the same final cells and phase
metadata from the unchanged upstream tutorial outputs. The upstream QC/PCA/
doublet/label-transfer steps were not repeated for this focused comparison.

- The complete **4,129 genes x 15,103 cells** count matrix is identical.
- Normalized RNA values and both phase-count tables are identical.
- All upstream per-cell metadata values are retained. Embedding-assay count
  and feature summaries are preserved as `embedding_nCount_RNA` and
  `embedding_nFeature_RNA`. The full-gene RNA assay computes its own
  `nCount_RNA` and `nFeature_RNA`.
- Whole-sample and all four phase DEG tables are identical to both the old
  construction and the previously saved tutorial results. Maximum differences
  in log2FC, raw P-value and Bonferroni-adjusted P-value are **zero**.
- The local test script, old/new code snapshots, log and machine-readable
  equivalence report are in `output/validation/merge_workflow/` and excluded
  from GitHub.
- A bounded read-only Claude Opus 5.5 review of the merge/metadata change found
  no blocking code defects. Numerical equivalence was checked separately in R.

## Additional Checks

- Toy examples checked the CP10k fold-change formula, the minimum-cell guard
  and refusal to use an integrated assay for DEG.
- Doublet-stage tests confirmed that calling doublets does not mutate the
  RNA-gated object, rerunning from that gated object is repeatable, and passing
  an already called object is rejected.
- The optional full-cell pK sweep ran successfully on a small synthetic object.
  A fresh full pK sweep of both real libraries was not repeated; the teaching
  run uses the previously reviewed sample-specific pK values.
- Every sample's phase proportions sum to one. Uncertain/unassigned display
  labels remain in the denominator. Both assignment-based and confidence-aware
  phase-count tables are saved.
- The newly transferred WT7 phase calls agree with the full research model for
  about **96.84%** of cells. They are deliberately not presented as identical
  clock-derived labels or as independent biological validation.
- Code, input provenance, package versions and detailed local review/run
  records are retained under `output/validation/` in the prepared local copy.
  Generated outputs and raw inputs are excluded from GitHub.

## Review and Limits

Claude Opus 5.5 independently reviewed the scripts and methods in read-only
mode. The review prompted distinct object names for each processing stage,
explicit backend control, input-gene-order checks and clearer phase/replication
caveats. The review is assistance, not a substitute for Kourosh's teaching QC.

The tested versions include R 4.4.2, Seurat 5.4.0, SeuratObject 5.3.0,
DoubletFinder 2.0.6, Matrix 1.7.4 and limma 3.62.2. R can report that some
cached binaries were built under R 4.4.3. Exact bitwise embeddings or borderline
doublet calls are not guaranteed on another operating system. Inspect changed
results rather than forcing the baseline counts.

Expected Seurat messages include replacing feature-name underscores with
dashes and integration warnings about corrected assays lacking count layers.
The DEG object is rebuilt separately from integer RNA counts. Unexpected RNA
input errors should not be dismissed as integration warnings.

This two-sample lesson does not reconstruct the advanced research clock,
include biological replication, reproduce the full slide deck or establish
arrest, absolute transcription suppression, survival or reactivation capacity.
