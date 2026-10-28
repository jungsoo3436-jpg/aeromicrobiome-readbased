# Upstream processing on the cluster

Everything in `analysis/` starts from the derived tables in `data/`. This file
records how those tables were produced, so the chain from raw reads to figures
is documented end to end.

All steps ran on the Digital Research Alliance of Canada cluster **Fir**, under
Slurm. Raw reads are in ENA under
[PRJEB109287](https://www.ebi.ac.uk/ena/browser/view/PRJEB109287).

Input: 15 paired-end libraries (2 x 150 bp, Illumina NovaSeq 6000), comprising
five composite aerosol metagenomes and ten paired environmental metagenomes
(two per site), plus two negative control libraries that did not yield usable
data.

---

## 1. Quality control

| Step | Tool | Version | Key parameters |
| --- | --- | --- | --- |
| Raw read QC | FastQC | as cited | defaults |
| Trimming and filtering | KneadData | 0.12.0 | Trimmomatic-based; sliding window 4:20, minimum length 50 bp |
| Post-filter QC | FastQC, MultiQC | as cited | summarised across samples |

`multiqc_raw_pre-filtering.html` and `multiqc_report_post-filtering.html` in the
parent project directory hold the summary reports.

## 2. Read-based taxonomic profiling

| Step | Tool | Version | Key parameters |
| --- | --- | --- | --- |
| Classification | Kraken2 | as cited | standard database (bacteria, archaea, viruses, plasmids, human, UniVec), downloaded August 2025 |
| Abundance re-estimation | Bracken | as cited | read length 150 bp, genus level |

Per-sample Bracken reports were merged with their taxonomy into
`data/bracken_with_taxonomy.tsv`, which is the input to
`analysis/01_build_phyloseq.R`. Sample columns in that file carry a
`.bracken_num` suffix, which script 01 strips.

## 3. Assembly and gene catalog construction

| Step | Tool | Version | Key parameters |
| --- | --- | --- | --- |
| Co-assembly, per site | MEGAHIT | 1.2.9 | one composite aerosol plus two paired environmental samples per site; minimum contig length 1,000 bp |
| ORF prediction | Prodigal | 2.6.3 | metagenomic mode (`-p meta`) |
| Clustering | CD-HIT | 4.8.1 | 95% amino acid identity |

Yield across the five site-level co-assemblies:

- 4,621,648 contigs at >= 1,000 bp
- 16,088,825 predicted ORFs
- 15,177,915 non-redundant protein-coding genes in the final catalog

Note that these are **site-level co-assemblies**, each pooling the aerosol
sample with its two paired environmental samples. The genome-resolved companion
analysis
([aeromicrobiome-mags](https://github.com/jungsoo3436-jpg/aeromicrobiome-mags))
uses the same co-assembly design.

## 4. Gene quantification and functional annotation

| Step | Tool | Version | Key parameters |
| --- | --- | --- | --- |
| Quantification | Salmon | 1.10.1 | quasi-mapping against the gene catalog as a reference transcriptome; sequence- and GC-bias correction enabled; output TPM |
| Annotation | eggNOG-mapper | 2.1.12 | eggNOG 5.0 database, DIAMOND 2.1.8 alignment, default confidence thresholds |

KEGG Orthology identifiers were retained for annotations meeting eggNOG-mapper
default thresholds. Annotations assigned to eukaryotic or viral lineages by
eggNOG taxonomic scope were excluded, restricting functional analyses to
prokaryotic gene content.

TPM values were converted to within-sample relative abundances, then summed by
KO to give the KO by sample matrix. That matrix and the KO annotation table are
in `data/gene_catalog/` and are the inputs to scripts 02 (Figure 5A) and 03
(Supplementary Figs 5-9).

---

## Version gaps

The Materials and Methods cite Kraken2, Bracken, FastQC and MultiQC without
version numbers, and those versions were not recorded in a form that survived
to this repository. The versions of every R package used downstream are
recorded exactly, in the `sessionInfo_*.txt` files written to `figures/`.

## Job scripts

The Slurm submission scripts are not included here. Each step above was run as a
single-tool job with the parameters given in the tables, and the tool versions
were installed through `mamba`. The ENA deposit plus the parameters above are
sufficient to reproduce the derived tables in `data/`.
