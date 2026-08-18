# Inferred/recovered plotting definitions for GALT paper Figure 2B-D
# DOI: 10.1126/sciimmunol.ady8948
#
# Evidence status
# - Subset names, marker order, row order, and plotting calls are recovered from
#   Fig1_04_Classification.R in the authors' public repository.
# - Subset colors are recovered from solid fills in the published PDF figure.
# - The ADT gradient is sampled from the published Figure 2C color bar because
#   config/globals.R is not present in the public repository.
# - The RNA gradient is made explicit as lightgrey -> blue, matching the Seurat
#   DotPlot default used by the authors (they supplied no replacement scale).

GALT_FIG2 <- list()

# Displayed top-to-bottom in Figure 2C and Figure 2D.
GALT_FIG2$subset_order_top_to_bottom <- c(
  "ActB3",
  "aNAV",
  "Centroblast",
  "Centrocyte",
  "CSM",
  "DN2",
  "IgM-only",
  "MZ",
  "Naive",
  "PB",
  "TS"
)

# Recovered from the solid subset fills in the attached published PDF.
# These are empirical PDF colors, not an invented categorical palette.
GALT_FIG2$BCELL_SUBSET <- c(
  "ActB3" = "#E4E1E3",
  "aNAV" = "#F8A19F",
  "Centroblast" = "#FBE426",
  "Centrocyte" = "#00FFFF",
  "CSM" = "#16CC32",
  "DN2" = "#1873CB",
  "IgM-only" = "#B10DA1",
  "MZ" = "#FEAF16",
  "Naive" = "#FE00FA",
  "PB" = "#F6222E",
  "TS" = "#325A9B"
)

# Low -> high. The continuous scale was visually recovered from the color bar
# in Figure 2C. Exact original stop positions remain uncertain because the
# defining globals.R file was not committed to the public repository.
GALT_FIG2$CYTOBANK_GRADIENT <- c(
  "#45135F",
  "#42536E",
  "#3E8C7B",
  "#96B947",
  "#FEE406",
  "#FAB30F",
  "#F58018",
  "#F14C1E",
  "#ED1C24"
)

# Figure 2C, left-to-right. The repository code uses sort(ANTIBODIES$ID).
GALT_FIG2$adt_features <- c(
  "ADT-CD10",
  "ADT-CD11C",
  "ADT-CD1C",
  "ADT-CD21",
  "ADT-CD27",
  "ADT-CD45RB",
  "ADT-CXCR5",
  "ADT-FCRL4",
  "ADT-IgD",
  "ADT-IgM"
)

# Figure 2D, left-to-right, copied from the author's plotting call.
GALT_FIG2$rna_features <- c(
  "TCL1A",
  "IGHD",
  "IGHM",
  "IGHA1",
  "IGHG1",
  "CD27",
  "CD1C",
  "CD24",
  "MALAT1",
  "MME",
  "BCL6",
  "SERPINA9",
  "MKI67",
  "CD38",
  "ITGAX",
  "FCRL4",
  "FCRL5",
  "ZEB2",
  "CR2"
)

GALT_FIG2$dot_color_limits <- c(-2.5, 2.5)
GALT_FIG2$dot_size_breaks <- c(0, 25, 50, 75, 100)

GALT_FIG2$theme_umap <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.45, color = "black"),
    axis.title = ggplot2::element_text(size = 10, color = "black"),
    legend.position = "none",
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

GALT_FIG2$theme_dotplot <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(color = "black"),
    axis.ticks = ggplot2::element_line(color = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.45, color = "black"),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )
