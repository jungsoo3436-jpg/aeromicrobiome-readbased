# aeromicrobiome-readbased

Analysis code and derived data for:

**Methane oxidation genes are only present in bioaerosols at sites with high methane concentrations**

Jungsoo Park, Nesil Han Martens, Sascha M.B. Krause, Anand Kumar,
Naomi Zimmerman, Nadine Borduas-Dedekind, S. Jane Fowler

Preprint: https://doi.org/10.64898/2026.10.13.758901 (posted 14 October 2026)

---

## What this study did

Aerosols and paired source environments (soil or water) were sampled at five
sites across the lower mainland of British Columbia spanning a 35-fold methane
range, from near-background wetlands (~2 ppm) to a dairy farm (26.8 ppm) and a
landfill (>=86 ppm, a lower bound above the analyser's calibration limit).
Ambient methane was measured concurrently at each site. Airborne
cell concentrations were quantified by flow cytometry, and shotgun
metagenomics with curated marker gene profiling resolved taxonomic composition
and functional gene representation across all 15 libraries (five composite
aerosol metagenomes and ten paired environmental metagenomes).

This repository covers the read-based and gene-catalog analyses: taxonomic
profiling, diversity, differential abundance, and curated marker gene panels.
The genome-resolved companion analysis lives in a separate repository,
[aeromicrobiome-mags](https://github.com/jungsoo3436-jpg/aeromicrobiome-mags).

## Repository layout

```
analysis/
  01_build_phyloseq.R                  Bracken table -> decontamination -> phyloseq -> rarefaction
  02_community_and_methane_figures.R   Figs 2, 3, 4, 5; Supp Figs 1, 2, 5; Tables S2, S3
  03_functional_marker_gene_panels.R   Supp Figs 6-10 (VOC/carbon, N, P, S, survival traits)
  04_archaea_virus_supplementary.R     Supp Fig 3 (archaeal and viral diversity and composition)
  05_domain_genus_composition.R        Supp Fig 4 (genus-level composition, four domains)
legacy/
  GSC_readbased_exploratory.R          Original exploratory session (see legacy/README.md)
data/
  bracken_with_taxonomy.tsv            Merged Kraken2/Bracken genus table
  metadata.txt                         Sample metadata
  phyloseq/                            Eight phyloseq objects, output of script 01
  gene_catalog/                        KO x sample proportional TPM, KO annotations
hpc/
  README.md                            Upstream steps run on the Fir cluster
figures/                               Figures and tables as published
```

## Reproducing the analysis

Scripts are numbered in dependency order and are run from the repository root,
each in a fresh R session. Paths are relative; nothing needs editing.

```r
setwd("/path/to/aeromicrobiome-readbased")
source("analysis/01_build_phyloseq.R")                 # optional, see below
source("analysis/02_community_and_methane_figures.R")
source("analysis/03_functional_marker_gene_panels.R")
source("analysis/04_archaea_virus_supplementary.R")
source("analysis/05_domain_genus_composition.R")
```

Script 01 regenerates the eight phyloseq objects from the Bracken table. Those
objects are already committed under `data/phyloseq/`, so scripts 02, 04 and 05
run without it. Rarefaction uses a fixed seed (24), and the NMDS ordinations use
seed 6, so repeated runs give identical output.

Scripts 02 and 03 read the KO by sample matrix in `data/gene_catalog/`, which
is the product of the assembly and gene-catalog steps described in
[`hpc/README.md`](hpc/README.md). Those steps need the raw reads from ENA and a
cluster; they are not reproducible from this repository alone.

### Software

Analyses were run in R 4.5.0, with DESeq2 1.50.2 for the differential abundance
analysis. Each script writes its own `sessionInfo` to `figures/`, which is the
authoritative record of the versions used.

```r
install.packages(c("tidyverse", "ggrepel", "patchwork", "ggh4x", "vegan"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "apeglm"))
```

## Figure and table mapping

| Manuscript item | Script | Output filename |
| --- | --- | --- |
| Figure 2A-B (flow cytometry) | 02 | `Fig2_flow_cytometry.pdf` |
| Figure 2C (TCC vs sequencing depth) | 02 | `Fig2C_TCC_vs_reads.pdf` |
| Figure 3A-E (diversity, NMDS, Procrustes) | 02 | `Fig3A`-`Fig3E_*.pdf`, plus combined layouts |
| Figure 4A-B (DESeq2 volcanoes) | 02 | `Fig4A_bacteria_volcano.pdf`, `Fig4B_eukaryota_volcano.pdf` |
| Figure 5A (C1 and trace gas heatmap) | 02 | `Fig5A_C1_metabolism_heatmap.pdf` |
| Figure 5B (pmoCAB vs ambient CH4) | 02 | `Fig5B_pmoCAB_vs_CH4.pdf` |
| Supplementary Fig. 1 (cells vs DNA and depth) | 02 | `SuppFig1_cells_DNA_reads.pdf` |
| Supplementary Fig. 2 (domain composition) | 02 | `SuppFig2_domain_composition.pdf` |
| Supplementary Fig. 3 (archaea and viruses) | 04 | `SuppFig3_archaea_virus.pdf` |
| Supplementary Fig. 4 (genus-level composition) | 05 | `SuppFig4_genus_composition_*.pdf` |
| Supplementary Fig. 5 (methanotroph genera) | 02 | `SuppFig5_methanotroph_genera_*.pdf` |
| Supplementary Figs 6-10 (marker gene panels) | 03 | `Supp_A_VOC_carbon.pdf` through `Supp_E_Survival.pdf` |
| Supplementary Table S1 (marker gene to KO map) | 03 | defined in the `tribble()` blocks of script 03 |
| Supplementary Table S2 (domain read proportions) | 02 | `TableS2_domain_read_proportions.csv` |
| Supplementary Table S3 (differentially abundant genera) | 02 | `TableS3_differentially_abundant_genera.csv` |

`figures/` holds the figures as they went into the manuscript, but not every
file in the table above: the Figure 3 panels, the Figure 4 volcanoes, the
Supplementary Fig. 3, 4 and 5 panels and Table S2 were saved elsewhere during
drafting and are not pre-committed. Running scripts 02 to 05 writes the
complete set into `figures/`.

`figures/Fig1_sampling_map.pdf` was prepared outside R.

Supplementary Table S1 is not written to a file. The marker genes and their KO
identifiers are the `tribble()` blocks at the head of each panel in script 03
and in Section 6 of script 02, which is where that table was compiled from.

## Methods summary

Read processing used KneadData 0.12.0 (Trimmomatic-based; sliding window 4:20,
minimum length 50 bp), with FastQC and MultiQC for quality assessment.

Taxonomic profiling used Kraken2 with the PlusPF database (the standard
database plus protozoa and fungi; downloaded August 2025) and Bracken at a read
length of 150 bp, analysed at genus level. Eukaryotic profiles are therefore
restricted to the fungi and protozoa represented in that database. Domains
were handled separately, with minimum relative abundance thresholds of 0.01%
(Bacteria, Archaea), 0.005% (Eukaryota) and 0.001% (Viruses) in at least one
sample. Alpha and beta diversity used counts rarefied to the per-domain minimum
library size (Bacteria 416,709; Archaea 965; Eukaryota 63,349; Viruses 511
reads). Community dissimilarity used Bray-Curtis distances visualised by NMDS,
tested with PERMANOVA (`adonis2`, 9,999 permutations) and `betadisper`.
Bacterial and eukaryotic community concordance was tested by Procrustes
analysis with PROTEST (999 permutations).

Differential abundance used DESeq2 on non-rarefied genus counts with a
site-blocked design (`~ Site + Compartment`), size factors by the positive
counts method, apeglm shrinkage, Benjamini-Hochberg correction, and thresholds
of adjusted p < 0.05 and |log2 fold change| > 2.

Gene-centric profiling aggregated Salmon gene-level TPM to KEGG Orthology
level, with functional categories and marker genes defined a priori. Category
abundances are sums of member KO proportional TPM. Heatmap tiles are within-site
z-scores for visualisation only; all statistical comparisons used
non-standardised proportional TPM.

See `hpc/README.md` for the upstream assembly, annotation and quantification
steps, and the Materials and Methods of the paper for full detail.

## Data availability

Raw metagenomic reads and the metagenome-assembled genomes are deposited
together in the European Nucleotide Archive under project accession
[PRJEB109287](https://www.ebi.ac.uk/ena/browser/view/PRJEB109287), released on
publication. Sample metadata are in the same project record.

This is the single deposit for both papers in the study. The MAGs in it are the
subject of the companion genome-resolved analysis,
[aeromicrobiome-mags](https://github.com/jungsoo3436-jpg/aeromicrobiome-mags);
the read-based analyses in this repository start from the same raw reads.

The derived tables in `data/` are included here so the figures can be
regenerated without re-running the cluster pipeline. See
[`data/README.md`](data/README.md) for what each file is and how it was made.

## Licence

Code is released under the MIT Licence (`LICENSE`). The derived data files in
`data/` and the figures in `figures/` are released under CC BY 4.0
(`LICENSE-data`).

## Funding

Government of Canada New Frontiers in Research Fund, NFRFE-2022-00163.

## Contact

Jungsoo Park, Department of Biological Sciences, Simon Fraser University
(jungsoop@sfu.ca); S. Jane Fowler (sjfowler@sfu.ca).
