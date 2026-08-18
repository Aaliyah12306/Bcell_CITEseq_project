###############################################################################
# Reproduce CITE-seq panels Figure 2B, 2C, and 2D
#
# Paper: Atlas of human gut-associated lymphoid tissue reveals
#        immunomodulatory interactions of B cells
# DOI:   10.1126/sciimmunol.ady8948
#
# Input:  output/GALT_CITEseq_Bcells.rds
# Output: output/figures/
#
# This script reproduces only the CITE-seq panels. It does not create CosMx,
# spatial, ligand-receptor, BCR, or IMC panels.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

if (packageVersion("Seurat") < package_version("5.0.0")) {
  stop("Seurat v5 or newer is required.")
}

set.seed(20260818)


# =============================================================================
# 1. Paths
# =============================================================================

project_dir <- normalizePath(".", mustWork = TRUE)
reference_file <- file.path(project_dir, "output", "GALT_CITEseq_Bcells.rds")
figures_dir <- file.path(project_dir, "output", "figures")
globals_file <- file.path(project_dir, "config", "globals_Fig2BCD_inferred.R")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(reference_file)) {
  stop("Reference object not found: ", reference_file)
}

if (!file.exists(globals_file)) {
  stop("Inferred globals file not found: ", globals_file)
}

source(globals_file)


# =============================================================================
# 2. Reviewable label mapping - do not silently guess broad B-cell labels
# =============================================================================

# The paper used classification.2 with exactly 11 labels. The supplied reference
# uses B_cell_subset. This dictionary maps only exact or biologically unambiguous
# aliases. Broad labels such as "GC B", "Memory B", "Regulatory B", and
# "Age-associated B" are deliberately not forced into paper classes.
USER_LABEL_TO_PAPER_SUBSET <- c(
  "ACTB3" = "ActB3",
  "ACTIVATEDB3" = "ActB3",
  "ACTIVATEDBCELL3" = "ActB3",
  "ANAV" = "aNAV",
  "ACTIVATEDNAIVEB" = "aNAV",
  "ACTIVATEDNAIVEBCELL" = "aNAV",
  "CENTROBLAST" = "Centroblast",
  "CENTROBLASTB" = "Centroblast",
  "CENTROCYTE" = "Centrocyte",
  "CENTROCYTEB" = "Centrocyte",
  "CSM" = "CSM",
  "CLASSSWITCHEDMEMORYB" = "CSM",
  "CLASSSWITCHEDMEMORYBCELL" = "CSM",
  "SWITCHEDMEMORYB" = "CSM",
  "DN2" = "DN2",
  "DN2B" = "DN2",
  "DOUBLENEGATIVE2B" = "DN2",
  "DOUBLENEGATIVE2BCELL" = "DN2",
  "IGMONLY" = "IgM-only",
  "IGMONLYB" = "IgM-only",
  "IGMONLYMEMORYB" = "IgM-only",
  "MZ" = "MZ",
  "MZB" = "MZ",
  "MARGINALZONEB" = "MZ",
  "MARGINALZONEBCELL" = "MZ",
  "NAIVE" = "Naive",
  "NAIVEB" = "Naive",
  "NAIVEBCELL" = "Naive",
  "PB" = "PB",
  "PLASMABLAST" = "PB",
  "PLASMABLASTS" = "PB",
  "TS" = "TS",
  "TRANSITIONAL" = "TS",
  "TRANSITIONALB" = "TS",
  "TRANSITIONALBCELL" = "TS"
)

# Optional cluster-level overrides. Use these only after reviewing RNA/ADT
# markers for the specified cluster. Example: c("3" = "Centrocyte").
USER_CLUSTER_TO_PAPER_SUBSET <- c()

# TRUE preserves biological fidelity: the script stops if any cell cannot be
# assigned to one of the 11 paper subsets. Set FALSE only for exploratory output;
# unmatched cells will then be omitted from panels 2B-D and logged in the audit.
STRICT_PAPER_SUBSETS <- TRUE


# =============================================================================
# 3. Helpers
# =============================================================================

normalize_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

match_features_exact_or_normalized <- function(targets, available, assay_name) {
  target_key <- normalize_key(targets)
  available_key <- normalize_key(available)
  index <- match(target_key, available_key)

  if (anyNA(index)) {
    stop(
      assay_name,
      " is missing required Figure 2 features: ",
      paste(targets[is.na(index)], collapse = ", ")
    )
  }

  matched <- available[index]
  names(matched) <- targets
  matched
}

save_panel <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(figures_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(figures_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}


# =============================================================================
# 4. Load and validate the reference
# =============================================================================

ref_b <- readRDS(reference_file)

required_assays <- c("SCT", "ADT")
missing_assays <- setdiff(required_assays, Assays(ref_b))

if (length(missing_assays) > 0L) {
  stop("Missing required assays: ", paste(missing_assays, collapse = ", "))
}

if (!"B_cell_subset" %in% colnames(ref_b[[]])) {
  stop("The reference lacks the required B_cell_subset metadata column.")
}

if ("wnn.umap" %in% Reductions(ref_b)) {
  umap_reduction <- "wnn.umap"
} else if ("umap" %in% Reductions(ref_b)) {
  warning("wnn.umap is absent; using umap. This is not an exact panel 2B match.")
  umap_reduction <- "umap"
} else {
  stop("Neither wnn.umap nor umap is present in the reference.")
}

adt_features <- match_features_exact_or_normalized(
  targets = GALT_FIG2$adt_features,
  available = rownames(ref_b[["ADT"]]),
  assay_name = "ADT assay"
)

rna_features <- match_features_exact_or_normalized(
  targets = GALT_FIG2$rna_features,
  available = rownames(ref_b[["SCT"]]),
  assay_name = "SCT assay"
)


# =============================================================================
# 5. Map B_cell_subset to the paper's classification.2 vocabulary
# =============================================================================

source_labels <- as.character(ref_b$B_cell_subset)
paper_subset <- unname(USER_LABEL_TO_PAPER_SUBSET[normalize_key(source_labels)])

if (
  length(USER_CLUSTER_TO_PAPER_SUBSET) > 0L &&
    "seurat_clusters" %in% colnames(ref_b[[]])
) {
  cluster_id <- as.character(ref_b$seurat_clusters)
  override <- unname(USER_CLUSTER_TO_PAPER_SUBSET[cluster_id])
  paper_subset[!is.na(override)] <- override[!is.na(override)]
}

invalid_override <- setdiff(
  unique(unname(USER_CLUSTER_TO_PAPER_SUBSET)),
  GALT_FIG2$subset_order_top_to_bottom
)

if (length(invalid_override) > 0L) {
  stop(
    "USER_CLUSTER_TO_PAPER_SUBSET contains invalid paper labels: ",
    paste(invalid_override, collapse = ", ")
  )
}

mapping_audit <- unique(data.frame(
  B_cell_subset = source_labels,
  paper_subset = paper_subset,
  stringsAsFactors = FALSE
))

mapping_audit <- mapping_audit[
  order(mapping_audit$B_cell_subset, mapping_audit$paper_subset),
  ,
  drop = FALSE
]

write.table(
  mapping_audit,
  file = file.path(figures_dir, "Fig2BCD_subset_mapping_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

unmapped_labels <- sort(unique(source_labels[is.na(paper_subset)]))

if (length(unmapped_labels) > 0L && STRICT_PAPER_SUBSETS) {
  stop(
    "The following B_cell_subset labels cannot be assigned to the paper's ",
    "11 classes without biological guessing: ",
    paste(unmapped_labels, collapse = ", "),
    ". Review Fig2BCD_subset_mapping_audit.tsv and add justified entries to ",
    "USER_LABEL_TO_PAPER_SUBSET or USER_CLUSTER_TO_PAPER_SUBSET."
  )
}

keep_cells <- !is.na(paper_subset)

if (!all(keep_cells)) {
  warning(sum(!keep_cells), " cells with unmapped labels will be omitted.")
  ref_b <- subset(ref_b, cells = colnames(ref_b)[keep_cells])
  paper_subset <- paper_subset[keep_cells]
}

ref_b$paper_subset <- factor(
  paper_subset,
  levels = rev(GALT_FIG2$subset_order_top_to_bottom)
)

observed_subsets <- unique(as.character(ref_b$paper_subset))
missing_paper_subsets <- setdiff(
  GALT_FIG2$subset_order_top_to_bottom,
  observed_subsets
)

if (length(missing_paper_subsets) > 0L) {
  warning(
    "The reference contains no cells for these paper subsets: ",
    paste(missing_paper_subsets, collapse = ", "),
    ". Their rows/labels cannot be drawn."
  )
}


# =============================================================================
# 6. Figure 2B - WNN UMAP colored by B-cell subset
# =============================================================================

p_2b <- DimPlot(
  object = ref_b,
  reduction = umap_reduction,
  group.by = "paper_subset",
  cols = GALT_FIG2$BCELL_SUBSET,
  pt.size = 0.15,
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  raster = ncol(ref_b) > 50000
) +
  scale_color_manual(values = GALT_FIG2$BCELL_SUBSET, drop = FALSE) +
  scale_fill_manual(values = GALT_FIG2$BCELL_SUBSET, drop = FALSE) +
  labs(x = "UMAP1", y = "UMAP2", title = NULL) +
  coord_fixed() +
  GALT_FIG2$theme_umap +
  NoLegend()

save_panel(
  plot = p_2b,
  stem = "Fig2B_CITEseq_Bcell_UMAP",
  width = 4,
  height = 4
)


# =============================================================================
# 7. Figure 2C - average ADT expression and percent expressed
# =============================================================================

p_2c <- DotPlot(
  object = ref_b,
  assay = "ADT",
  features = unname(adt_features),
  group.by = "paper_subset",
  scale = TRUE,
  scale.min = 0,
  col.min = GALT_FIG2$dot_color_limits[1],
  col.max = GALT_FIG2$dot_color_limits[2],
  dot.min = 0
) +
  scale_x_discrete(labels = setNames(names(adt_features), unname(adt_features))) +
  scale_color_gradientn(
    colours = GALT_FIG2$CYTOBANK_GRADIENT,
    limits = GALT_FIG2$dot_color_limits,
    breaks = GALT_FIG2$dot_color_limits,
    labels = c("Low", "High"),
    oob = scales::squish,
    name = "Average ADT\nexpression"
  ) +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = GALT_FIG2$dot_size_breaks,
    range = c(0, 6),
    name = "Percent\nexpressed"
  ) +
  guides(
    color = guide_colorbar(
      order = 1,
      title.position = "top",
      barheight = grid::unit(1.1, "in"),
      barwidth = grid::unit(0.18, "in")
    ),
    size = guide_legend(order = 2, title.position = "top")
  ) +
  GALT_FIG2$theme_dotplot

save_panel(
  plot = p_2c,
  stem = "Fig2C_CITEseq_ADT_DotPlot",
  width = 6.4,
  height = 4.2
)


# =============================================================================
# 8. Figure 2D - average RNA expression and percent expressed
# =============================================================================

p_2d <- DotPlot(
  object = ref_b,
  assay = "SCT",
  features = unname(rna_features),
  group.by = "paper_subset",
  scale = TRUE,
  scale.min = 0,
  col.min = GALT_FIG2$dot_color_limits[1],
  col.max = GALT_FIG2$dot_color_limits[2],
  dot.min = 0
) +
  scale_x_discrete(labels = setNames(names(rna_features), unname(rna_features))) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue",
    limits = GALT_FIG2$dot_color_limits,
    breaks = GALT_FIG2$dot_color_limits,
    labels = c("Low", "High"),
    oob = scales::squish,
    name = "Average RNA\nexpression"
  ) +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = GALT_FIG2$dot_size_breaks,
    range = c(0, 6),
    name = "Percent\nexpressed"
  ) +
  guides(
    color = guide_colorbar(
      order = 1,
      title.position = "top",
      barheight = grid::unit(1.1, "in"),
      barwidth = grid::unit(0.18, "in")
    ),
    size = guide_legend(order = 2, title.position = "top")
  ) +
  GALT_FIG2$theme_dotplot

save_panel(
  plot = p_2d,
  stem = "Fig2D_CITEseq_RNA_DotPlot",
  width = 9.0,
  height = 4.2
)


# =============================================================================
# 9. Provenance and session information
# =============================================================================

donor_columns <- intersect(
  c("donor", "donor_id", "sample", "sample_id", "orig.ident", "status"),
  colnames(ref_b[[]])
)

provenance <- c(
  "GALT Figure 2B-D CITE-seq reproduction",
  paste0("Reference: ", reference_file),
  paste0("UMAP reduction: ", umap_reduction),
  "Grouping column supplied by reference: B_cell_subset",
  "Grouping used for plots: paper_subset (explicitly audited mapping)",
  paste0("Available sample metadata: ", paste(donor_columns, collapse = ", ")),
  paste0("Cells plotted: ", ncol(ref_b)),
  paste0("Subsets observed: ", paste(sort(observed_subsets), collapse = ", ")),
  paste0(
    "Subsets absent: ",
    if (length(missing_paper_subsets) == 0L) "none" else
      paste(missing_paper_subsets, collapse = ", ")
  ),
  "No statistical annotations are added because published panels 2B-D show none.",
  paste0("Generated: ", Sys.time())
)

writeLines(
  provenance,
  con = file.path(figures_dir, "Fig2BCD_provenance.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(figures_dir, "Fig2BCD_sessionInfo.txt")
)

message("Figure 2B-D reproduction completed: ", figures_dir)
