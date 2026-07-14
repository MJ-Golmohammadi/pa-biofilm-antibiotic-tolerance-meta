# Integrated Transcriptomic Meta-Analysis of Antibiotic Tolerance in *Pseudomonas aeruginosa* Biofilms

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

This repository contains the complete computational workflow, analysis scripts, processed datasets, and supplementary materials associated with the study:

> **Integrated transcriptomic meta-analysis reveals conserved molecular signatures underlying antibiotic tolerance in *Pseudomonas aeruginosa* biofilms**

The study integrates multiple publicly available transcriptomic datasets to identify conserved transcriptional signatures associated with antibiotic tolerance in *Pseudomonas aeruginosa* biofilms. The analytical workflow combines differential expression analysis, robust rank aggregation (RRA), functional enrichment, protein–protein interaction (PPI) network analysis, and network module detection to identify conserved adaptive mechanisms.

---

## Data Availability

All transcriptomic datasets analyzed in this study were obtained from the **NCBI Gene Expression Omnibus (GEO)**.

The accession numbers are reported in the manuscript.


## Analysis Workflow

The computational workflow consists of the following steps:

1. Download GEO datasets
2. Data preprocessing and quality assessment
3. Differential gene expression analysis
4. Robust Rank Aggregation (RRA)
5. Functional enrichment analysis
6. Protein–protein interaction network construction
7. Network module detection (MCODE)
8. Hub gene identification
9. Figure generation

---

## Software Requirements

The analysis was performed in **R (version 4.3.3)**.

Major R packages include:

- GEOquery
- affy
- limma
- RobustRankAggreg
- clusterProfiler
- igraph
- tidyverse
- ggplot2
- pheatmap

Network visualization and module detection were performed using:

- Cytoscape
- MCODE

---

## Reproducibility

All scripts are fully documented and can be executed sequentially to reproduce the analyses presented in the manuscript.

Processed intermediate files are provided to facilitate reproducibility and reduce computational time.

---

## Supplementary Materials

This repository includes:

- Supplementary Figures
- Supplementary Tables
- Processed datasets
- Analysis scripts

---

## Citation

If you use this repository, please cite the associated manuscript.

**Preprint**

*(Link will be added after publication on bioRxiv.)*

---

## Correspondence

**Mohammad Javad Golmohammadi**

Independent Researcher

Tehran, Iran

ORCID: https://orcid.org/0000-0002-9277-0023

Email: Mohammad.jg75@gmail.com

---

## License

This project is distributed under the MIT License.
