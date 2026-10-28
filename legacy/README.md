# Legacy exploratory script

`GSC_readbased_exploratory.R` is the original interactive RStudio script used
while exploring the read-based data. It is kept here unmodified for provenance.

It is **not** the reproducible pipeline. Use `analysis/01` through `analysis/04`
instead, which are consolidated from it and from the gene-catalog script, carry
relative paths, and run start to finish in a fresh session.

This script is retained for two reasons:

1. It is where the phyloseq construction, rarefaction, DESeq2 volcano and
   methanotroph genus blocks now in `analysis/01` and `analysis/02` originally
   came from, so the consolidated code can be checked against it.
2. It contains exploratory analyses that are not in the paper: virus-specific
   NMDS, Chao1 and Pielou evenness panels, IndVal indicator genus analysis,
   per-domain genus-level taxonomy barplots, and a gene-centric methane panel.

Caveats if you run it: it uses `setwd()` with absolute paths, reassigns
`project_data` between domains so output file names were set by hand, and
depends on objects defined earlier in the session. It will not run unattended.
