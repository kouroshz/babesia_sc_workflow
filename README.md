# Babesia scRNA-seq Workflow

A line-by-line RStudio lesson for Molly and Aryan using original, untagged
*Babesia divergens* WT at 37 C (0 h) and after seven days at 4 C.

Start with **[run_tutorial.R](run_tutorial.R)**. It contains the scientific
choices and the main Seurat commands. **[tutorial_functions.R](tutorial_functions.R)**
contains the short repeated operations. There is no hidden pipeline runner.

This is a **two-sample teaching workflow**, not the complete CIQ paper analysis
or all 61 slides. There is one biological sample per condition. Cell-level DEG
P-values do not supply biological replication.

## Start in RStudio

1. Clone this private repository and open **Babesia_scRNA_Tutorial.Rproj**.
2. Run **setup.R** once to restore the pinned packages, then restart R.
3. Obtain the eight count/report files using the Terminal commands at the top
   of **run_tutorial.R**. Use your own authorized Markov and Chimera accounts.
4. Open **run_tutorial.R** and run sections 1-3 to examine input and RNA QC.
5. Continue one section at a time. Pause at each figure and table before moving
   on. The first lesson should not start by sourcing the entire script.

Select a few lines and use **Cmd+Enter** on macOS or **Ctrl+Enter** on
Windows/Linux. Keep the Environment, Console and Plots panes visible. The
Terminal download commands assume macOS/Linux; Windows users need an rsync
installation such as WSL or an approved alternative transfer client.

The project does not load an old workspace. Opening it activates `renv` but
does not run an analysis. Package restoration is deliberate. The tested R
version is **4.4.2**, with Seurat **5.4.0** and DoubletFinder **2.0.6**.
Avoid updating packages during the lesson.

The Wilcoxon backend is also fixed: `limma` 3.62.2 is included. The code uses
Seurat's `test.use="wilcox_limma"`, explicitly selecting the same implementation
used by CIQ's Wilcoxon tests when `presto` was absent. Installing an optional
backend cannot silently change this choice.

A separate project library was restored from the pinned lockfile and tested
locally. Most packages came from the existing renv cache. A second operating
system has not been tested. Some Linux installations need compilers and HDF5
development libraries before `hdf5r` can install.

For a complete fresh rerun after reviewing the code:

```r
source("run_tutorial.R")
```

Every fresh run starts from the count matrices and overwrites this tutorial's
own output filenames. It never reads the main CIQ project's processed objects.
Copy the `output/` folder first if you want to preserve an earlier trial.

## Lesson Sequence

| Section | What you do | What to inspect |
|---|---|---|
| 0 | Download the two 10x outputs | Original HTML reports and file paths |
| 1-2 | Load packages, settings and integer UMI matrices | Gene rows, called cells and raw barcode counts |
| 3 | Plot detected genes versus UMIs and apply explicit RNA gates | Dashed cutoffs and cells lost |
| 4 | Normalize WT0, select variable genes, scale and fit PCA | PCA and elbow plots |
| 5 | Fit neighbors, clusters and independent UMAPs | Each sample's native structure |
| 6 | Call doublets independently, rebuild singlets and refit | Doublet locations and cell/gene counts |
| 7 | Attach the supplied WT0 annotations and transfer to WT7 | Native phase UMAPs and transfer confidence |
| 8 | Reference-based RPCA integration | Joint and sample-split UMAPs, alongside native views |
| 9 | Merge full-gene RNA objects, join layers and inspect markers | Published phase markers and state proportions |
| 10 | Compare all retained cold and warm cells | Whole-sample DEG table and volcano |
| 11 | Compare cold and warm within G, SM, MC and C | Cell denominators and up/down DEG counts |
| 12 | Save objects and tables for the next question | Readable CSVs and reusable RDS outputs |

## Folder Map

```text
babesia_sc_workflow/
  Babesia_scRNA_Tutorial.Rproj
  run_tutorial.R                  the lesson and visible controls
  tutorial_functions.R           repeated processing and plotting operations
  setup.R                        one-time package restoration
  README.md
  VALIDATION.md                  completed checks and known limits
  renv.lock                      pinned package versions
  renv/activate.R
  input/
    WT0/                         downloaded H5s and 10x reports, ignored by Git
    WT7/                         downloaded H5s and 10x reports, ignored by Git
    data_files.csv               exact source paths, sizes and SHA-256 checksums
    gene_annotation.csv          release-68 gene descriptions
    WT0_phase_labels.csv         explicit small annotated-reference input
    phase_markers.csv            published genes used for marker inspection
    reference_sources.csv        provenance of the small reference inputs
  output/
    figures/                     numbered PNG and PDF figures
    tables/                      QC counts and DEG tables
    objects/                     generated Seurat objects and DEG results
    sessionInfo.txt
```

Only code, documentation and small reference tables are committed. Count
matrices, vendor reports, generated figures/tables and RDS objects remain local.

## Data and Filtering

The updated alignment uses **PiroplasmaDB-68 Bdivergens1802A**. The count-matrix
input is Cell Ranger's **filtered** H5, which contains its called barcodes.
The raw H5 is used only to report pre-cell-calling barcode counts. A raw barcode
column is not a confirmed cell. No BAM/FASTQ is needed for this workflow.

Both remote locations and full `rsync` commands appear in section 0. These are
comments and never execute from R. Transfers use the approved Markov jump host.
Run the analysis locally or on an allocated compute node, not on an HPC login
node. The remote files were not rechecked during this teaching handoff. Their
exact locally staged versions and checksums are in `input/data_files.csv`.

| Sample | Detected-gene gate | Raw-UMI gate | Reviewed pK | Baseline final cells |
|---|---|---|---|---|
| WT0, 37 C | >200 | >=250 | 0.005 | 9,532 |
| WT7, 4 C for 7 d | >75 | >=100 | 0.010 | 5,571 |

These are fixed project choices, not universal Seurat thresholds. No extra
EmptyDrops gate, upper-count cutoff, fabricated mitochondrial percentage or
real-cell downsampling is applied. All cell gates use the complete input matrix.

DoubletFinder assumes a **5%** doublet loading rate, adjusted for homotypic
similarity from the sample's clusters. This is not a measured rate. The supplied
pK values came from reviewed BCmvn sweeps using all real cells and pN values
0.10, 0.20 and 0.30. The optional `choose_doublet_pk()` function provides that
full-cell sweep for a new sample. It is slower and not part of the routine
two-sample lesson. Do not assume the WT pK values fit a new library.

After doublet removal, the code rebuilds each sample from integer counts and
refits normalization, variable features, scaling, PCA, neighbors and UMAP.
Genes detected in at least three final singlets define the embedding matrix.
The final QC table also reports genes detected anywhere and raw-UMI summaries.
If your counts differ from the baseline, inspect them rather than forcing the
expected answer.

Object names distinguish stages: `wt0_all` is the complete 10x-called input,
`wt0_gated` passes the RNA gates, `wt0_called` carries doublet calls, and `wt0`
contains final singlets. WT7 follows the same naming. Re-running a gate or
doublet section therefore starts from the appropriate earlier-stage object.
After changing a setting, rerun that stage and every downstream stage.

## Where the Phase Labels Come From

`WT0_phase_labels.csv` is an **explicit processed annotation input**, not a
hidden Seurat object. It contains the project reference RNA-phase assignment
for the 9,532 baseline WT0 barcodes. Its exact source and checksum are recorded
in `reference_sources.csv`.

Those project labels came from the approved WT0 cyclic reference, bulk-adjusted
clock and Babesia atlas phase definitions. This tutorial **does not rebuild
that advanced clock or infer its boundaries from scratch**. It teaches how to
use an annotated reference and transfer its labels to a new sample. The full
research workflow remains in the separate CIQ project.

WT7 labels are recomputed with `FindTransferAnchors()` and `TransferData()`
using the newly processed WT0 RNA reference. They are not looked up by WT7
barcode. Because the teaching reference uses the newly processed WT0 sample
rather than the full research clock construction, the WT7 phase assignments
need not be identical to the presentation's clock-derived assignments.

G, SM, MC and C are **Babesia transcriptomic phase labels**, not mammalian
G1/S/G2/M measurements. Low-confidence WT7 transfers appear gray in display
figures. The phase DEG uses the assigned phase for every labeled cell, without
silently applying an extra confidence filter. Any newly admitted WT0 cells
without supplied reference labels remain visibly unassigned.

WT0 labels are supplied reference assignments, whereas WT7 labels are
predictions. Different boundary errors between these methods can contribute to
apparent phase-matched differences. Phase matching reduces one source of
composition differences but does not remove this uncertainty.

The marker table contains the top three available genes per phase, ordered by
published `avg_log2FC`, from the selected *B. divergens*, human-host markers in
[the Babesia atlas, S3 Table](https://doi.org/10.1371/journal.pbio.3001816).
The dot plot is a useful consistency check, not an independent validation of
labels that were developed using the same marker resource.

## Integration and DEG

Independent UMAPs come first. Joint RPCA integration then provides a separate
comparison view, using WT0 as the reference. Integration can remove genuine
cold effects because condition and library coincide in this two-sample example.
**Integrated expression is never used for DEG.**

For DEG, a separate full-gene Seurat object is rebuilt for each sample using
only its final singlets and all its cell metadata. Standard Seurat `merge()`
combines those objects without integration. `JoinLayers()` then joins the
sample count layers before normalization and `FindMarkers()`. Sample IDs and
phase labels remain in metadata. No cells or counts are summed across samples.

Restoring all input gene rows prevents the embedding-only gene filter from
changing the test denominator or RNA library totals. The earlier embedding
object's `nCount_RNA` and `nFeature_RNA` are preserved as
`embedding_nCount_RNA` and `embedding_nFeature_RNA`; the new RNA assay computes
its own summaries over the complete gene set. Other metadata, including
doublet scores and available transfer-confidence scores, are retained.
Saved cluster IDs remain sample-specific: WT0 cluster 0 is not automatically
the same population as WT7 cluster 0.

The direction is always **WT7 cold minus WT0**. Seurat's limma-backed Wilcoxon
test supplies P-values, with `min.pct=0.10`, `logfc.threshold=0` and no cell cap. The reported
project fold change uses linear per-cell counts normalized to 10,000:

```text
log2FC = log2[(mean cold CP10k + 1) / (mean warm CP10k + 1)]
```

This effect is calculated for every gene, including those failing Seurat's
detection-frequency test. It is kept separate from Seurat's returned
`avg_log2FC` definition. A called DEG requires **Bonferroni-adjusted P <0.05**
and **absolute log2FC >log2(1.5)**. Bonferroni uses the full input-gene count
within each comparison. No additional correction across phases is applied.

The whole-sample comparison uses all retained cells. The phase analysis repeats
the same comparison within each transferred phase, requiring at least 20 cells
per condition. NT means not tested, not zero DEGs. The tables distinguish
tested genes, nonsignificant genes and calls in each direction.

These are exploratory **cell-level** screens from one sample per condition,
not biological-replicate-aware bulk RNA-seq. Normalization also does not test
absolute transcription rate, viability or cold survival. The two samples alone
cannot reproduce the original cold-induction plus recovery-reversal candidate
ranking, which requires the recovery sample.

Whole-sample differences also combine condition, phase composition, RNA
complexity and selection by the different warm/cold gates. More down calls
alone do not establish absolute cold-induced repression or global RNA shutdown.

## Continue Your Own Exploration

```r
wt7 <- readRDS("output/objects/WT7_phase_transferred.rds")
FeaturePlot(wt7, features = "Bdiv-031680c")

deg <- readRDS("output/objects/DEG_results.rds")$whole_sample
View(subset(deg, called))
```

Keep the original lesson unchanged while learning. Add a new, clearly named
script for your next question. For a new library, use its own RNA QC, pK choice,
doublet calls and independent embedding, then apply the reference-transfer
section. Experimental metadata and biological replication must be reconsidered
before interpreting a new contrast.

## Method References

- [Seurat guided workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Seurat reference mapping](https://satijalab.org/seurat/articles/integration_mapping.html)
- [Pinned DoubletFinder source](https://github.com/chris-mcginnis-ucsf/DoubletFinder/tree/1b244d8f0d54b4b1cb4365639931bbb16f01e1cd)
- [Babesia atlas and supplementary phase markers](https://doi.org/10.1371/journal.pbio.3001816)

Current online Seurat documentation can describe a newer version. This lesson
uses the versions in `renv.lock`. See [VALIDATION.md](VALIDATION.md) for what was
actually rerun and checked before sharing.
