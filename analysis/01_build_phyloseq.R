###############################################################################
# Aerosol metagenome manuscript (read-based community analysis)
# Step 01: Build phyloseq objects from the Kraken2/Bracken genus table
#
# Author: Jungsoo Park (Simon Fraser University)
#
# What this does
#   1. Reads the merged Bracken genus table and splits taxonomy from counts
#   2. Builds a phyloseq object with the sample metadata
#   3. Drops the two negative control libraries (AC_DF2, AC_RM2)
#   4. Removes human, organellar and synthetic-construct taxa
#   5. Splits by domain and applies the per-domain abundance thresholds
#      reported in the Methods (Bacteria and Archaea 0.01%, Eukaryota 0.005%,
#      Viruses 0.001%, in at least one sample)
#   6. Rarefies each domain to its own minimum library size and saves both the
#      non-rarefied and rarefied objects
#
# Outputs (written to data/phyloseq/)
#   psf_{bacteria,archaea,eukaryota,viruses}_nonrarefied.rds
#   psf_{bacteria,archaea,eukaryota,viruses}_raref_min.rds
#
# These eight .rds files are already committed to the repository, so steps 02
# and 04 can be run without re-running this script. Run it to reproduce them
# from the Bracken table.
#
# Run from a FRESH R session, from the repository root.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(phyloseq)
})

DATA_DIR <- "data"
OUT_DIR  <- file.path(DATA_DIR, "phyloseq")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RNG_SEED <- 24   # rarefaction seed, fixed for reproducibility

# =============================================================================
# 1. Read the Bracken genus table
# =============================================================================

df <- read.delim(file.path(DATA_DIR, "bracken_with_taxonomy.tsv"),
                 sep = "\t", quote = "", check.names = FALSE,
                 stringsAsFactors = FALSE)

# "taxa" holds the stable Bracken identifier and is used as the row name
tax_cols  <- c("taxa", "name", "taxid", "superkingdom", "phylum", "class",
               "order", "family", "genus", "species")
rank_cols <- c("superkingdom", "phylum", "class", "order", "family",
               "genus", "species")
stopifnot(all(tax_cols %in% colnames(df)))

sample_cols <- setdiff(colnames(df), tax_cols)

counts <- as.matrix(df[, sample_cols, drop = FALSE])
mode(counts) <- "numeric"
rownames(counts) <- df$taxa
colnames(counts) <- sub("\\.bracken_num$", "", colnames(counts))

tax_mat <- as.matrix(df[, rank_cols])
rownames(tax_mat) <- df$taxa

cat("Bracken table:", nrow(counts), "taxa x", ncol(counts), "samples\n")

# =============================================================================
# 2. Metadata and phyloseq object
# =============================================================================

meta <- read.delim(file.path(DATA_DIR, "metadata.txt"), sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
if ("Sample" %in% colnames(meta)) rownames(meta) <- meta$Sample

missing_in_metadata <- setdiff(colnames(counts), rownames(meta))
if (length(missing_in_metadata)) {
  stop("Samples absent from metadata.txt: ",
       paste(missing_in_metadata, collapse = ", "))
}
meta <- meta[colnames(counts), , drop = FALSE]

ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE),
               tax_table(tax_mat),
               sample_data(meta))

# =============================================================================
# 3. Drop the negative control libraries
#    AC_DF2 and AC_RM2 are the field controls (Control == "Negative" in
#    metadata.txt). They are excluded from the community analyses, leaving the
#    15 libraries reported in the paper: five aerosol and ten paired
#    environmental metagenomes.
# =============================================================================

ac_samples <- grep("^AC", sample_names(ps), value = TRUE)
ps_noAC <- prune_samples(!sample_names(ps) %in% ac_samples, ps)

cat("Dropped", length(ac_samples), "negative control libraries (",
    paste(ac_samples, collapse = ", "), ");",
    nsamples(ps_noAC), "samples retained\n")
stopifnot(nsamples(ps_noAC) == 15)

# =============================================================================
# 4. Remove human, organellar and synthetic-construct taxa
# =============================================================================

tx <- as.data.frame(unclass(tax_table(ps_noAC)))
tx$taxa <- rownames(tx)

ab <- data.frame(taxa        = taxa_names(ps_noAC),
                 total_reads = taxa_sums(ps_noAC),
                 row.names   = taxa_names(ps_noAC))

cand <- tx %>%
  left_join(ab, by = "taxa") %>%
  mutate(any_field = apply(cbind(superkingdom, phylum, class, order,
                                 family, genus, species), 1,
                           function(z) paste(z, collapse = " "))) %>%
  filter(
    grepl("\\bHomo sapiens\\b|\\bHomo\\b|\\b9606\\b", any_field, ignore.case = TRUE) |
      grepl("mitochondri|chloroplast", any_field, ignore.case = TRUE) |
      grepl("phix|phiX174|lambda phage|control|vector|synthetic construct|plasmid",
            any_field, ignore.case = TRUE)
  ) %>%
  arrange(desc(total_reads)) %>%
  select(taxa, superkingdom, phylum, class, order, family, genus, species,
         total_reads)

print(head(cand, 30))

ps_clean <- prune_taxa(!taxa_names(ps_noAC) %in% cand$taxa, ps_noAC)
cat("Removed", ntaxa(ps_noAC) - ntaxa(ps_clean), "flagged taxa;",
    ntaxa(ps_clean), "retained\n")

# =============================================================================
# 5. Split by domain and apply per-domain abundance thresholds
# =============================================================================

relabund_mat <- as(otu_table(transform_sample_counts(ps_clean,
                                                     function(x) x / sum(x))),
                   "matrix")
dom_of <- tax_table(ps_clean)[, "superkingdom"]

# Threshold = minimum relative abundance in at least one sample (see Methods)
filter_domain <- function(domain, threshold) {
  rows <- relabund_mat[dom_of == domain, , drop = FALSE]
  keep <- apply(rows, 1, function(x) any(x >= threshold))
  prune_taxa(names(keep[keep]), ps_clean)
}

psf_bacteria  <- filter_domain("Bacteria",  1e-4)   # 0.01%
psf_archaea   <- filter_domain("Archaea",   1e-4)   # 0.01%
psf_eukaryota <- filter_domain("Eukaryota", 5e-5)   # 0.005%
psf_viruses   <- filter_domain("Viruses",   1e-5)   # 0.001%

for (nm in c("bacteria", "archaea", "eukaryota", "viruses")) {
  obj <- get(paste0("psf_", nm))
  cat(sprintf("%-10s: %5d taxa | reads = %d\n",
              nm, ntaxa(obj), sum(sample_sums(obj))))
}

saveRDS(psf_bacteria,  file.path(OUT_DIR, "psf_bacteria_nonrarefied.rds"))
saveRDS(psf_archaea,   file.path(OUT_DIR, "psf_archaea_nonrarefied.rds"))
saveRDS(psf_eukaryota, file.path(OUT_DIR, "psf_eukaryota_nonrarefied.rds"))
saveRDS(psf_viruses,   file.path(OUT_DIR, "psf_viruses_nonrarefied.rds"))

# =============================================================================
# 6. Rarefy each domain to its own minimum library size
#    Per-domain depths retain all 15 libraries while standardising effort
#    within each domain. Archaeal and viral depths are low, so the diversity
#    estimates for those domains are interpreted cautiously in the paper.
# =============================================================================

depths_min <- list(Bacteria  = 416709L,
                   Archaea   =    965L,
                   Eukaryota =  63349L,
                   Viruses   =    511L)

rare_one <- function(obj, depth, seed = RNG_SEED) {
  obj_sub <- prune_samples(sample_sums(obj) >= depth, obj)
  phyloseq::rarefy_even_depth(obj_sub, sample.size = depth, rngseed = seed,
                              trimOTUs = TRUE, verbose = FALSE)
}

set.seed(RNG_SEED)
psf_bacteria_raref_min  <- rare_one(psf_bacteria,  depths_min$Bacteria)
psf_archaea_raref_min   <- rare_one(psf_archaea,   depths_min$Archaea)
psf_eukaryota_raref_min <- rare_one(psf_eukaryota, depths_min$Eukaryota)
psf_viruses_raref_min   <- rare_one(psf_viruses,   depths_min$Viruses)

stopifnot(
  nsamples(psf_bacteria_raref_min)  == 15,
  nsamples(psf_archaea_raref_min)   == 15,
  nsamples(psf_eukaryota_raref_min) == 15,
  nsamples(psf_viruses_raref_min)   == 15,
  all(sample_sums(psf_bacteria_raref_min)  == depths_min$Bacteria),
  all(sample_sums(psf_archaea_raref_min)   == depths_min$Archaea),
  all(sample_sums(psf_eukaryota_raref_min) == depths_min$Eukaryota),
  all(sample_sums(psf_viruses_raref_min)   == depths_min$Viruses)
)

saveRDS(psf_bacteria_raref_min,  file.path(OUT_DIR, "psf_bacteria_raref_min.rds"))
saveRDS(psf_archaea_raref_min,   file.path(OUT_DIR, "psf_archaea_raref_min.rds"))
saveRDS(psf_eukaryota_raref_min, file.path(OUT_DIR, "psf_eukaryota_raref_min.rds"))
saveRDS(psf_viruses_raref_min,   file.path(OUT_DIR, "psf_viruses_raref_min.rds"))

cat("\nDone. Eight phyloseq objects written to", OUT_DIR, "\n")
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo_01.txt"))
