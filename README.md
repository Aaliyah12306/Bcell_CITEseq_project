# Bcell_CITEseq_project
# B-cell CITE-seq Annotation for Ileal Crohn’s Disease — Project README

Version: 1.1  
Working language: English for scientific content; commands in R.  
Main environment: R 4.4.1, Seurat v5, RStudio.

## 1. Project purpose

Use the GALT/CITE-seq reference framework from:

> Atlas of human gut-associated lymphoid tissue reveals immunomodulatory interactions of B cells  
> DOI: https://www.science.org/doi/10.1126/sciimmunol.ady8948  
> Code: https://github.com/jspencer-lab/galt_cosmX_imc  
> Public CITE-seq data: GSE327417

to:

1. Reproduce B-cell CITE-seq clustering from GSE327417.
2. Build strictly reference-grounded marker tables for ileal Crohn’s disease B-cell subsets and other intestinal lineages.
3. Transfer reference B-cell labels to the user’s existing Crohn’s disease scRNA-seq B-cell clusters.
4. Score each user cluster for likely B-cell subset identity.

## 2. AI operating rules / token efficiency

- Read only necessary files. Do not print whole Seurat objects or large tables in chat.
- Save all outputs as `.csv`, `.tsv`, `.rds`, `.md` files under the project directory.
- When producing marker tables, cite references as `PMID` or `DOI`, not full citations.
- Do not rely on model memory for marker genes or references. Use web search and verify each marker source.
- Keep one task per major session. Reference file paths instead of re-pasting content.
- Local data path:  
  `Bcell_CITEseq_project/Pitcher2026/raw_extracted/`

## 3. Primary reference: required reading and extraction

### Task 1 — Literature structure and multi-omics integration

Read the main text, Methods, supplementary files, and reference list of the paper above. Produce `summary/paper_structure.md` containing:

- Study design and sample types.
- Data modalities used: scRNA-seq, CITE-seq, spatial transcriptomics, imaging mass cytometry, etc.
- Key B-cell subsets identified.
- How CITE-seq RNA + ADT data were integrated.
- How spatial data validated B-cell interactions.
- Main immunomodulatory interaction claims.
- Reproducibility assets: GEO accession, GitHub repo, analysis scripts.

#### Multi-omics integration deep-dive

While the primary framework is the GALT paper, **supplement the integration strategy section with additional methodological references** from other fields. This is not to replace the GALT approach, but to learn transferable integration logic.

Search for and briefly summarize integration strategies from high-impact studies in areas such as:

- Multimodal single-cell integration: RNA + protein, RNA + ATAC, RNA + spatial.
- CITE-seq WNN-like approaches outside gut immunology.
- Spatial transcriptomics + scRNA-seq deconvolution or mapping.
- Imaging mass cytometry + single-cell reference mapping.
- Cross-modality label transfer and anchor-based integration.
- Multi-omics factor analysis approaches, if relevant.

For each additional reference, record in `summary/multiomics_integration_references.md`:

`study | field | modalities | integration method | key takeaway | PMID_or_DOI`

**Constraint:** The primary analytical backbone must remain the GALT paper. Additional references are only used to enrich the methodological discussion and improve the integration logic. Do not deviate from the GALT-based analysis workflow.

## 4. Marker table requirements

### Task 2 — Ileal Crohn’s disease B-cell subset markers

Produce:

`tables/B_cell_markers_ileal_CD.tsv`

Columns:

`subset | marker | type | evidence | PMID_or_DOI | notes`

`type` = `mRNA` or `protein/ADT`.

**Mandatory requirements:**

- Every single marker must have at least one supporting PMID or DOI.
- You are strictly forbidden from including any cell subset or marker that lacks a verifiable literature citation.
- Do not add markers from memory. If no citation is found, exclude the marker.
- Do not include "commonly used" markers without a specific reference.
- Prioritize single-cell or CITE-seq studies of Crohn's disease, ileum, or gut-associated lymphoid tissue.
- For recently reported B-cell states, the defining publication must be cited.
- If a subset is controversial or newly described, clearly state this in the `notes` column.

Cover classical and recent B-cell subsets:

- Naive B
- Transitional B
- Germinal center B
- Unswitched memory B
- Switched memory B
- Double-negative memory B
- Plasmablast
- Plasma cell
- IgA+ plasma cell
- Regulatory B cell / IL-10+ B cell
- Age-associated B cell / CD11c+ B cell
- Any recently reported inflammatory or gut-specific B-cell states in ileal Crohn’s disease

Search strategy:

- Use PubMed/Europe PMC.
- Prioritize journals with IF > 10: Nature, Science, Immunity, Cell, Nature Medicine, Nature Immunology, Gastroenterology, Gut, Journal of Experimental Medicine, Cell Host & Microbe.
- Prefer single-cell or CITE-seq studies of Crohn’s disease, ileum, or gut-associated lymphoid tissue.

### Task 3 — Non-B lineage markers

Produce:

`tables/lineage_markers_ileal.tsv`

Same column format.

**Mandatory requirements:**

- Every marker must have a PMID or DOI from the primary literature.
- Any marker without a verifiable citation must be removed.
- Do not include a cell subset unless all listed markers have been validated in cited references.
- If a classical subset has conflicting marker definitions across studies, list the markers with their respective citations, and note the discrepancy.

Cover at least:

- Stromal cells: fibroblasts（ADAMDEC1+ fibroblasts,RSPO3+ fibroblasts,FAP+ fibroblasts,etc), myofibroblasts, pericytes, smooth muscle cells(RERGL+ SMCs,etc)
- T cells: CD4+ T, CD8+ T, naive T, memory T, Th17, Treg, Tfh, IELs
- Endothelial cells: blood endothelial(venous, arterial, capillary), lymphatic endothelial, Tips
- Epithelial cells: stem, absorptive enterocyte, goblet, Paneth, enteroendocrine, tuft
- Myeloid cells: monocytes, macrophages, dendritic cells, neutrophils, mast cells

**Verification rule:**

Before finalizing the table, run a self-check: for each row, confirm that the PMID or DOI actually supports the marker–subset association. If not, delete the row.

## 5. Data and environment

- R: 4.4.1
- Seurat: v5
- Local raw data:  
  `Bcell_CITEseq_project/Pitcher2026/raw_extracted/`
- Expected GEO format: 10x Genomics-like files. Inspect first with:

```r
list.files("Bcell_CITEseq_project/Pitcher2026/raw_extracted", recursive = TRUE)
