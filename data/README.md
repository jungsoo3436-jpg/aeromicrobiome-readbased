# Derived data

These are the derived tables and objects the analysis scripts read. Raw
sequencing reads are not here; they are in ENA under
[PRJEB109287](https://www.ebi.ac.uk/ena/browser/view/PRJEB109287).

Everything in this folder is released under CC BY 4.0 (see `../LICENSE-data`).

## Sample naming

Fifteen libraries, named `<type>_<site><replicate>`:

- `A_` aerosol (one composite sample per site, pooled from nine impinger runs)
- `S_` soil, `W_` water (two paired environmental replicates per site)

Sites, ordered from high to low ambient methane:

| Code | Site | CH4 (ppm) | Paired environment |
| --- | --- | --- | --- |
| VL | Vancouver Landfill | 85.978 | Soil |
| DF | Dairy Farm | 26.759 | Soil |
| MK | Minnekhada Park | 2.624 | Soil |
| BB | Burns Bog | 2.597 | Water |
| RM | Richmond Marsh | 2.291 | Soil |

The aerosol library at Burns Bog is `A_BB1` and at the landfill `A_VL2`; the
numbering reflects which impinger pool yielded a usable library, not a
replicate. The two Simon Fraser University and Burnaby Refinery aerosol samples
yielded insufficient DNA for library construction and appear only in the flow
cytometry panels, where their values are hard-coded in
`analysis/02_community_and_methane_figures.R`.

## Files

### `bracken_with_taxonomy.tsv` (3.1 MB, 14,972 taxa)

Merged Kraken2/Bracken genus-level abundance estimates with taxonomy. Columns
`taxa`, `name`, `taxid` and the seven ranks `superkingdom` through `species`,
followed by one count column per library named `<sample>.bracken_num`. The
column set includes the two negative control libraries `AC_DF2` and `AC_RM2`,
which `01_build_phyloseq.R` drops.

Produced by step 2 of [`../hpc/README.md`](../hpc/README.md).

### `metadata.txt` (17 rows)

Sample metadata. `Sample_type` is Air, Soil or Water; `Sample_type2` collapses
this to Air versus Reference, which is the grouping used in the statistical
comparisons. `Control` marks the two negative control libraries as Negative.
`CH4` is the ambient methane concentration in ppm at the time of aerosol
sampling, measured by cavity ring-down spectroscopy, and is repeated across all
samples from a site. `Temp`, `Humidity` and `Wind_speed` are the field
conditions recorded at sampling.

### `phyloseq/` (8 objects)

Output of `01_build_phyloseq.R`, committed so that scripts 02 and 04 run without
re-running step 01.

- `psf_{bacteria,archaea,eukaryota,viruses}_nonrarefied.rds`
  Domain-split, decontaminated, abundance-filtered genus counts. Used for the
  DESeq2 differential abundance analysis, which models raw counts and must not
  receive rarefied input.
- `psf_{bacteria,archaea,eukaryota,viruses}_raref_min.rds`
  The same objects rarefied to the per-domain minimum library size (Bacteria
  416,709; Archaea 965; Eukaryota 63,349; Viruses 511 reads), seed 24. Used for
  alpha and beta diversity.

Abundance thresholds applied before rarefaction were 0.01% for Bacteria and
Archaea, 0.005% for Eukaryota and 0.001% for Viruses, in at least one sample.
The archaeal and viral rarefaction depths are low, which is why the diversity
estimates for those two domains are supplementary and are interpreted
cautiously.

### `gene_catalog/KOxSample_TPM.prop.tsv` (1.2 MB, 8,795 KOs x 15 samples)

KEGG Orthology by sample matrix. Values are within-sample relative abundances
of Salmon TPM, summed across all catalog genes assigned to the same KO. Rows are
KO identifiers, columns are the 15 libraries. These are the values used for all
statistical comparisons; the within-site z-scores in the heatmaps are computed
from this matrix at plotting time and are for visualisation only.

### `gene_catalog/KO_to_Function.clean3.tsv` (8.4 MB, 146,493 rows)

KO annotation table, no header, three tab-separated columns: KO identifier,
preferred gene name, description. Used to restrict the KO matrix to annotated
orthologs and to label the heatmap rows.

Both gene catalog files are products of steps 3 and 4 of
[`../hpc/README.md`](../hpc/README.md).

## What is not here

- Raw and quality-filtered reads: ENA PRJEB109287.
- Assembled contigs and the 15.2 M gene non-redundant protein catalog: too
  large to version here. Regenerable from the ENA reads with the parameters in
  `../hpc/README.md`.
- Flow cytometry event files. The per-run total and intact cell counts used in
  Figure 2 and Supplementary Fig. 1 are hard-coded in script 02, as reported in
  the paper.
