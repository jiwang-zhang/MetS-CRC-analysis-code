# Analysis code and intermediate results for the study of shared genetic candidates between metabolic syndrome and colorectal cancer

> **Note on the title.** The manuscript title at the time of this deposit still reads
> "…share GPAT3 and PLAU as candidate mediators…". The revised version of the manuscript
> replaces "mediators" with the term "shared candidates". The authoritative title is the one
> given on the associated article; this README deliberately avoids restating it.

This record contains the analysis code and the intermediate result files underlying the
findings reported in the manuscript above. All analyses were performed *in silico* on
publicly available data; **no wet-laboratory experiment was conducted, and no patient-level
identifiable data were used or deposited here**.

---

## Contents

| Directory | Contents |
|---|---|
| `01_MR/` | Mendelian randomisation: instrument clumping, harmonisation, forward (MetS → CRC) and reverse (CRC → MetS) analyses, heterogeneity and pleiotropy diagnostics |
| `02_DEG_and_ML/` | Differential expression for the CRC meta-cohort and the MetS dataset, co-DEG intersection, TNF/NF-κB prior-set construction, GO/KEGG enrichment, and the three machine-learning feature-selection algorithms (LASSO, random forest, SVM-RFE) |
| `03_single_cell_and_immune/` | Single-cell reanalysis, CIBERSORT immune deconvolution, linkET correlation analysis, and *in silico* virtual knockout (scTenifoldKnk) with enrichment |
| `04_druggability_and_docking/` | Drug-signature enrichment against DSigDB and molecular docking against the PLAU catalytic domain |
| `05_external_validation_GEPIA/` | External expression validation using GEPIA2 (TCGA COAD/READ versus GTEx normal tissue) |
| `07_results/` | Intermediate result files: MR estimates and pleiotropy tests, differential-expression tables for both contrasts, the 200-gene co-DEG list and prior-set intersections, hub–immune correlations, and the GPAT3 knockout footprint |
| `08_supplementary/` | Supplementary Tables S1–S4 as submitted |
| `09_manuscript/` | Markdown source of the manuscript **as originally submitted**. Revised wording will be deposited as a subsequent version of this record |

## Public data used (not redistributed here)

All input data are public and can be obtained from the following accessions:

| Resource | Accession | Used in |
|---|---|---|
| Gene Expression Omnibus | GSE21510, GSE32323, GSE39582 | CRC meta-cohort differential expression |
| Gene Expression Omnibus | GSE98895 | MetS contrast differential expression |
| Gene Expression Omnibus | GSE10950 | independent validation cohort |
| Gene Expression Omnibus | GSE132465 | single-cell reanalysis and virtual knockout |
| FinnGen | C3_COLORECTAL (DF13) | colorectal cancer outcome |
| van Walree et al. | common genetic factor GWAS of MetS components | MetS exposure |
| TCGA / GTEx | COAD, READ versus normal colon and rectum | external expression validation |
| MSigDB | HALLMARK_TNFA_SIGNALING_VIA_NFKB (release 2026.1.Hs) | TNF/NF-κB prior set |
| Enrichr / DSigDB | Drug Signatures Database | druggability screening |
| Protein Data Bank | 5ZA9 | molecular docking target |
| CIBERSORT | LM22 signature matrix | immune deconvolution |

## Software

Analyses were performed in R and Python. Principal packages used were `TwoSampleMR`,
`MendelianRandomization`, `MRPRESSO`, `limma`, `clusterProfiler`, `glmnet`, `randomForest`,
`e1071`, `Seurat`, `scTenifoldKnk`, `CIBERSORT` and `AutoDock Vina` (via the CB-Dock2
server). Exact package versions are recorded in each script header where they affect a
reported result; the R and Python versions used were 4.x and 3.x respectively.

## How to reproduce

Each directory preserves the sub-folder structure of the original working analysis, so the
order of operations is given by the numbered prefixes. Scripts assume the public data listed
above are available locally and paths at the top of each script are adjusted to point at
them. Running the scripts in numeric order within a directory reproduces the corresponding
intermediate result files deposited in `07_results/`.

## Licence

MIT — see `LICENSE`.

## Citation

If you use this code or these results, please cite the associated manuscript. The archival
DOI is displayed on the Zenodo record page for this repository.
