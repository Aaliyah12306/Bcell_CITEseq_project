# Paper structure: Atlas of human gut-associated lymphoid tissue reveals immunomodulatory interactions of B cells

## Citation and scope

King HW et al. *Science Immunology* (2026), **Atlas of human gut-associated lymphoid tissue reveals immunomodulatory interactions of B cells**. DOI: [10.1126/sciimmunol.ady8948](https://doi.org/10.1126/sciimmunol.ady8948); PMID: [42247482](https://pubmed.ncbi.nlm.nih.gov/42247482/).

The study builds a normal human appendix/GALT B-cell reference and maps it into intact tissue to identify spatially constrained B-cell interactions. Its inflammatory comparison is **ulcerative colitis (UC) appendix**, not ileal Crohn's disease (CD). Consequently, it is a valuable reference for cell-state transfer into the present CITE-seq project, but disease-specific states must be checked against ileal-CD cohorts rather than assumed to transfer unchanged.

## Study design and samples

- **CITE-seq reference:** CD19-positive B cells isolated from appendices of three deceased adult organ donors. The experiment generated paired 5-prime gene-expression, antibody-derived-tag (ADT), and immunoglobulin V(D)J libraries using 10x Genomics chemistry.
- **Spatial transcriptomics:** seven formalin-fixed, paraffin-embedded normal adult appendices profiled with the CosMx 1000-plex Human Universal Cell Characterization Panel.
- **Imaging validation:** independent normal appendix sections were assessed by immunofluorescence, confocal microscopy, imaging mass cytometry (IMC), and RNAscope. Four normal and four severe-UC appendices were compared by IMC for disease-associated spatial and protein changes.
- **External reference:** an adult-appendix subset of a published human intestinal single-cell atlas was used to assign broad CosMx lineages before the study's B-cell labels were transferred.

The design separates three questions: high-resolution B-cell state definition by CITE-seq; localization and cell-neighborhood inference by spatial transcriptomics; and orthogonal confirmation of key proteins, transcripts, and UC-associated architectural changes by imaging.

## Modalities and multi-omics integration

### 1. Paired RNA and ADT processing

RNA and protein were normalized separately. Per-sample RNA data were processed with SCTransform, while ADTs were centered-log-ratio normalized. Batch correction was performed independently for the two modalities. For RNA, variable genes were selected after excluding immunoglobulin and dissociation-associated genes; reciprocal-PCA anchors were then used for SCT integration. ADTs were integrated using reciprocal-PCA anchors in protein-PC space.

### 2. Weighted-nearest-neighbor B-cell atlas

The integrated RNA and ADT representations were joined using Seurat weighted-nearest neighbors (WNN). The released analysis uses RNA PCs 1-30 and ADT PCs 1-8, followed by WSNN clustering at resolution 2.0. Initial clusters were manually merged or subclustered using paired transcript and surface-protein evidence. V(D)J libraries supplied immunoglobulin-isotype context but were not the primary clustering representation.

### 3. Broad spatial cell-type transfer

The CosMx data were quality controlled, including removal of fields/cells affected by dispersed background transcripts. A published adult-appendix scRNA-seq reference was restricted to relevant lineages and approximately 480 overlapping genes. Seurat anchor-based label transfer assigned broad immune, mesenchymal, endothelial, and other cell types to spatial cells.

### 4. B-cell subset transfer and transcript imputation

Spatially classified B cells were isolated and mapped to the CITE-seq reference using 238 genes shared between the CosMx panel and the CITE-seq integration feature set. The same anchors transferred B-cell labels and imputed the larger CITE-seq RNA feature space into CosMx B cells. The imputed values enabled ligand-receptor screening beyond the directly measured 1000-gene panel, but validation was therefore essential for the highest-priority inferred signals.

### 5. Spatial-neighborhood and interaction analysis

A Delaunay graph linked nearby cells, with long edges removed using a 125-micrometre cutoff. Cell-type proximity was compared with 1,000 shuffled spatial configurations. Candidate ligand-receptor pairs came from CellTalkDB and were retained when expressed across at least 10% of relevant edges; the analysis was then focused on interactions involving B-subset marker genes. Metascape was used for pathway interpretation.

## B-cell subsets resolved

| Subset | Defining evidence in the study | Interpretation |
|---|---|---|
| Transitional (TS) | CD10-positive surface phenotype | Immature/transitional compartment entering the mature pool |
| Naive | IgD-positive, CD27-low conventional naive profile | Baseline follicular/naive compartment |
| Activated naive (aNAV) | CD11c protein and FCRL5 RNA | Activated-naive state with a strong interaction profile |
| Germinal-center centroblast | BCL6 program with Ki67 positivity | Proliferating dark-zone-like GC cells |
| Germinal-center centrocyte | BCL6 program with low/absent Ki67 | Nonproliferating light-zone-like GC cells |
| Marginal-zone-like (MZ) | High CD27 and IgM; NOTCH2-associated biology | GALT MZ-like/unswitched-memory compartment |
| IgM-only memory | IGHM with lower IgD and CD1c than MZ cells | Separate unswitched-memory-like state in this atlas |
| Class-switched memory (CSM) | CD27 with non-IgM/IgD isotypes | Switched memory compartment |
| Double-negative 2 (DN2) | Low CD27, IgD, CXCR5, and CD21; high CD11c; FCRL4/FCRL5 | GALT-enriched age-associated/atypical-memory-like state |
| Activated B 3 (ActB3) | MALAT1-high paper-defined cluster | A reference-specific activated state; MALAT1 is not sufficiently specific for stand-alone annotation |
| Plasmablast (PB) | MZB1, SLAMF7, abundant immunoglobulin; reduced MS4A1/CD19 | Antibody-secreting effector compartment |

The most distinctive result is the prominence of aNAV and especially DN2 cells in the inferred interaction network. DN2 cells were enriched in the subepithelial dome/follicle-associated epithelium region and were positioned near T cells, myeloid cells, and mesenchymal cells. Their high FCRL4 expression supported an inhibitory/atypical-memory-like phenotype.

## Immunomodulatory interactions and biological claims

The analysis connected B-cell location with candidate communication programs. Notable examples included TNFRSF13B (TACI)-TNFSF13B (BAFF), CD86-CTLA4, GRN-TNFRSF1A/B, and EBI3 interactions with IL27RA/IL6ST. Antigen-presentation genes such as HLA molecules, B2M, and CD74 also contributed to the DN2 interaction profile. Other subset-linked signals included naive-B-cell LTA-TNFR signaling, consistent with lymphoid architecture; MZ-cell NOTCH2-DLL1 signaling; CCR7-CCL19/CCL21 in MZ, CSM, and IgM-only cells; and CCR6-CCL20 in MZ/IgM-only cells.

These are spatially informed candidate interactions, not direct demonstrations of ligand-receptor causality. The paper strengthens the claims by validating selected molecules and by showing that their protein/transcript patterns and cell neighborhoods change in inflamed UC appendices.

## Validation strategy

1. **Orthogonal spatial assays:** IMC measured multiplexed proteins, while RNAscope localized selected transcripts in intact appendix. Conventional immunofluorescence and confocal microscopy supplied additional validation in independent sections.
2. **Key DN2 markers and mediators:** DN2 identity was confirmed using CD11c with FCRL4. CD86 protein and GRN/EBI3 spatial signals were validated, supporting the principal immunomodulatory claims. SIGLEC10 was not experimentally validated because a suitable reagent was unavailable and should therefore remain an imputed candidate.
3. **Independent disease comparison:** a semisupervised random-forest pipeline labeled IMC cells in four severe-UC and four control appendices. UC tissue showed altered B/T organization, displacement of DN2 cells, and changes in CD86, GRN, and EBI3 patterns. This validates inflammation-associated remodeling in UC, but it is not a Crohn-specific validation.
4. **Reproducible computational assets:** raw/processed data are deposited under GEO accession [GSE327417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327417); analysis code is available at [jspencer-lab/GALT_CosMX_IMC](https://github.com/jspencer-lab/GALT_CosMX_IMC); and the associated archive is available as [Zenodo record 18405128](https://zenodo.org/records/18405128).

## Implications for ileal-Crohn CITE-seq annotation

- Use the paired RNA+ADT appendix atlas as a high-quality reference for GALT B-cell states, especially MZ-like, IgM-only, aNAV, GC, DN2, and plasmablast compartments.
- Preserve modality-specific evidence during transfer: surface CD11c/FCRL4/CD27/IgD can resolve states that overlap transcriptionally, while RNA provides BCL6/MKI67, FCRL5, isotype, and secretory programs.
- Treat label transfer as a starting point. Ileal CD can expand IgG plasma cells, inflammatory myeloid-stromal modules, CCR7-high CD4 T cells, and fibrotic fibroblast states not represented by healthy appendix CITE-seq.
- Validate DN2/ABC-like and regulatory-like calls with marker panels and spatial/functional evidence. CD86, GRN, and EBI3 describe an immunomodulatory DN2 state in this paper; they do not by themselves define classical IL-10-producing regulatory B cells.
