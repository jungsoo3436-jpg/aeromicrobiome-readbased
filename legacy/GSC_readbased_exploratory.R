library(tidyverse)
library(reshape2)
library(stringr)
library(ape)
library(phyloseq)
library(data.table)
library(viridis)
library(qualpalr)
library(ggplot2)
library(vegan)
library(dplyr)
library(tibble)
library(gridExtra)

# Set working directory
setwd("~/Desktop/Airmicrobiome_project/GSC_aeromicrobiome/R_work/")

# -------------------------------------------
# Step 1: Read in Bracken output file
# -------------------------------------------
# Path where you scp’d the file
fn <- "bracken_with_taxonomy.tsv"

# Read strictly as tab-delimited; keep original column names
df <- read.delim(fn, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)

# Inspect the first columns to confirm layout
head(df[ , 1:12])
ncol(df)

## Step 2: identify taxonomy columns vs sample columns
tax_cols <- c("taxa","name","taxid","superkingdom","phylum","class","order","family","genus","species")
stopifnot(all(tax_cols %in% colnames(df)))

sample_cols  <- setdiff(colnames(df), tax_cols)
sample_names <- sub("\\.bracken_num$", "", sample_cols)  # strip suffix

# Quick peek to confirm mapping
data.frame(original = head(sample_cols, 10),
           cleaned  = head(sample_names, 10))
length(sample_cols)

## Build counts matrix
counts <- as.matrix(df[, sample_cols])
mode(counts) <- "numeric"                           # ensure numeric
rownames(counts) <- df$taxa                         # use stable IDs as rownames
colnames(counts) <- sub("\\.bracken_num$", "", colnames(counts))

## Quick checks
dim(counts)
sum(is.na(counts))
head(rowSums(counts))
colSums(counts)

## Build taxonomy data.frame
tax <- df[, tax_cols]

## Keep rownames aligned with counts
rownames(tax) <- df$taxa

## Convert to matrix for phyloseq
tax_mat <- as.matrix(tax[, c("superkingdom","phylum","class","order","family","genus","species")])

## Quick peek
head(tax_mat)

## Step 5: load metadata and construct phyloseq

# 5a) Read metadata (expects a Sample column or rownames already set)
meta <- read.delim("metadata.txt", sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
if ("Sample" %in% colnames(meta)) rownames(meta) <- meta$Sample

# 5b) Check name alignment
missing_in_metadata <- setdiff(colnames(counts), rownames(meta))
missing_in_metadata

# 5c) Reorder metadata to match counts (drop any extras in metadata)
meta <- meta[colnames(counts), , drop = FALSE]

# 5d) Build phyloseq
library(phyloseq)
ps <- phyloseq(
  otu_table(counts, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta)
)
ps

## starting from your existing object `ps`

# which samples start with "AC"?
ac_samples <- grep("^AC", sample_names(ps), value = TRUE)
ac_samples

# drop them
ps_noAC <- prune_samples(!sample_names(ps) %in% ac_samples, ps)

# quick checks
length(sample_names(ps)); length(sample_names(ps_noAC))  # before vs after
setdiff(sample_names(ps), sample_names(ps_noAC))         # confirm which were removed
summary(sample_sums(ps_noAC))                            # library sizes now

## STEP A — Inspect likely contaminants in ps_noAC
library(dplyr)
tx <- as.data.frame(unclass(tax_table(ps_noAC)))
tx$taxa <- rownames(tx)

# add total abundance across samples
ab <- data.frame(taxa = taxa_names(ps_noAC),
                 total_reads = taxa_sums(ps_noAC),
                 row.names = taxa_names(ps_noAC))

cand <- tx %>%
  left_join(ab, by = "taxa") %>%
  mutate(any_field = apply(cbind(superkingdom, phylum, class, order, family, genus, species), 1,
                           function(z) paste(z, collapse = " "))) %>%
  # flag patterns: human, organelles, spike-ins/vectors
  filter(
    grepl("\\bHomo sapiens\\b|\\bHomo\\b|\\b9606\\b", any_field, ignore.case = TRUE) |
      grepl("mitochondri|chloroplast", any_field, ignore.case = TRUE) |
      grepl("phix|phiX174|lambda phage|control|vector|synthetic construct|plasmid", any_field, ignore.case = TRUE)
  ) %>%
  arrange(desc(total_reads)) %>%
  select(taxa, superkingdom, phylum, class, order, family, genus, species, total_reads)

# show the top 30 hits (if any)
head(cand, 30)

## STEP B — Drop contaminants

contam_taxa <- cand$taxa   # all flagged candidates from STEP A
ps_clean <- prune_taxa(!taxa_names(ps_noAC) %in% contam_taxa, ps_noAC)

# quick check
ntaxa(ps_noAC); ntaxa(ps_clean)    # before vs after
nsamples(ps_clean)                 # sample count should stay same

# Done! You can now inspect this object or continue with diversity analysis
View(ps)
# For specific table view
View(as.data.frame(unclass(sample_data(ps))))
View(as.data.frame(unclass(otu_table(ps))))
View(as.data.frame(unclass(tax_table(ps))))

## Relative abundance
relabund <- transform_sample_counts(ps_clean, function(x) x / sum(x))
relabund_mat <- as(otu_table(relabund), "matrix")

## Bacteria
keep_bac <- apply(relabund_mat[tax_table(ps_clean)[,"superkingdom"]=="Bacteria", ], 1,
                  function(x) any(x >= 0.0001))   # 0.01%
psf_bacteria <- prune_taxa(names(keep_bac[keep_bac]), ps_clean)

## Archaea
keep_arch <- apply(relabund_mat[tax_table(ps_clean)[,"superkingdom"]=="Archaea", ], 1,
                   function(x) any(x >= 0.0001))  # 0.01%
psf_archaea <- prune_taxa(names(keep_arch[keep_arch]), ps_clean)

## Eukaryota
keep_euk <- apply(relabund_mat[tax_table(ps_clean)[,"superkingdom"]=="Eukaryota", ], 1,
                  function(x) any(x >= 0.00005))  # 0.005%
psf_eukaryota <- prune_taxa(names(keep_euk[keep_euk]), ps_clean)

## Viruses
keep_vir <- apply(relabund_mat[tax_table(ps_clean)[,"superkingdom"]=="Viruses", ], 1,
                  function(x) any(x >= 0.00001))  # 0.001%
psf_viruses <- prune_taxa(names(keep_vir[keep_vir]), ps_clean)

cat("Bacteria :", ntaxa(psf_bacteria),  "| reads =", sum(sample_sums(psf_bacteria)), "\n")
cat("Archaea  :", ntaxa(psf_archaea),   "| reads =", sum(sample_sums(psf_archaea)),  "\n")
cat("Eukaryota:", ntaxa(psf_eukaryota), "| reads =", sum(sample_sums(psf_eukaryota)),"\n")
cat("Viruses  :", ntaxa(psf_viruses),   "| reads =", sum(sample_sums(psf_viruses)),  "\n")


## === SAVE FILTERED OBJECTS ===

## 1) Save the non-rarefied, domain-split, filtered objects
saveRDS(psf_bacteria,  "psf_bacteria_nonrarefied.rds")
saveRDS(psf_archaea,   "psf_archaea_nonrarefied.rds")
saveRDS(psf_eukaryota, "psf_eukaryota_nonrarefied.rds")
saveRDS(psf_viruses,   "psf_viruses_nonrarefied.rds")


## Rarefy each domain at the per-domain MIN depth 
#(performed at the minimum per-domian depth to retain all 14 samples while standardizing sequencing effort)
depths_min <- list(
  Bacteria  = 416709L,
  Archaea   = 965L,
  Eukaryota = 63349L,
  Viruses   = 511L
)

rare_one <- function(obj, depth, seed = 24){
  obj_sub <- prune_samples(sample_sums(obj) >= depth, obj)   # should be all 14 at min
  phyloseq::rarefy_even_depth(
    obj_sub, sample.size = depth, rngseed = seed, trimOTUs = TRUE, verbose = FALSE
  )
}

set.seed(24)
psf_bacteria_raref_min  <- rare_one(psf_bacteria,  depths_min$Bacteria)
psf_archaea_raref_min   <- rare_one(psf_archaea,   depths_min$Archaea)
psf_eukaryota_raref_min <- rare_one(psf_eukaryota, depths_min$Eukaryota)
psf_viruses_raref_min   <- rare_one(psf_viruses,   depths_min$Viruses)

# Quick sanity (should all be 14 samples and exactly the chosen depths)
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

# Save RDS (rarefied, min-depth)
saveRDS(psf_bacteria_raref_min,  "psf_bacteria_raref_min.rds")
saveRDS(psf_archaea_raref_min,   "psf_archaea_raref_min.rds")
saveRDS(psf_eukaryota_raref_min, "psf_eukaryota_raref_min.rds")
saveRDS(psf_viruses_raref_min,   "psf_viruses_raref_min.rds")


################################################################################
################################################################################
###START FROM HERE#
################################################################################
################################################################################
# Read RDS after all triming
project_data <- readRDS("psf_bacteria_nonrarefied.rds")
project_data <- readRDS("psf_archaea_nonrarefied.rds")
project_data <- readRDS("psf_eukaryota_nonrarefied.rds")
project_data <- readRDS("psf_viruses_nonrarefied.rds")


library(phyloseq)
library(dplyr)
library(ggplot2)
library(patchwork)

# =========================================================
# SETTINGS
# =========================================================
N_TOP     <- 15
RDS_FILES <- list(
  Bacteria  = "psf_bacteria_raref_min.rds",
  Archaea   = "psf_archaea_raref_min.rds",
  Eukaryota = "psf_eukaryota_raref_min.rds",
  Viruses   = "psf_viruses_raref_min.rds"
)

# Fixed sample order across all panels: aerosol first, then paired-env, grouped by site
SAMPLE_ORDER <- c(
  "A_VL2",
  "A_DF1",
  "A_MK1",
  "A_BB1",
  "A_RM1",
  "S_VL1", "S_VL2",
  "S_DF1", "S_DF2",
  "S_MK1", "S_MK2",
  "W_BB1", "W_BB2",
  "S_RM1", "S_RM2"
)

# Shared colour palette (up to 15 genera + Other)
BASE_COLS <- c(
  '#4E79A7','#F28E2B','#59A14F','#E15759','#76B7B2',
  '#EDC948','#B07AA1','#FF9DA7','#9C755F','#BAB0AC',
  '#D37295','#FABFD2','#8CD17D','#86BCB6','#499894',
  '#a9a9a9'  # Other
)

# =========================================================
# HELPER: build tidy relative abundance df for one domain
# =========================================================

build_df <- function(rds_path, domain_label) {
  
  ps <- readRDS(rds_path)
  
  # Collapse to genus
  ps_gen   <- tax_glom(ps, taxrank = "genus", NArm = TRUE)
  
  # Clean genus labels
  tx       <- as.data.frame(tax_table(ps_gen))
  tx$genus <- trimws(as.character(tx$genus))
  tx$genus[is.na(tx$genus) | tx$genus == ""] <- "Unclassified"
  tx$genus <- sub("\\s+sp\\.?.*$", "", tx$genus)
  tax_table(ps_gen) <- tax_table(as.matrix(tx))
  
  # Relative abundance
  df <- ps_gen %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    mutate(genus = as.character(genus))
  
  # Top N genera globally
  top_genera <- df %>%
    group_by(genus) %>%
    summarise(total = sum(Abundance), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = N_TOP) %>%
    pull(genus)
  
  # Collapse rest to Other, renormalize within sample
  df <- df %>%
    mutate(genus_plot = ifelse(genus %in% top_genera, genus, "Other")) %>%
    group_by(Sample, genus_plot) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup() %>%
    mutate(
      domain      = domain_label,
      genus_plot  = factor(genus_plot,
                           levels = c(top_genera[top_genera %in% unique(genus_plot)], "Other"))
    )
  
  return(list(df = df, top_genera = top_genera))
}

# =========================================================
# BUILD DATA FOR ALL DOMAINS
# =========================================================

domains   <- names(RDS_FILES)
all_data  <- lapply(domains, function(d) build_df(RDS_FILES[[d]], d))
names(all_data) <- domains

# =========================================================
# HELPER: make one barplot panel
# =========================================================

make_panel <- function(df, top_genera, domain_label,
                       show_x = FALSE, sample_order = SAMPLE_ORDER) {
  
  # Assign colours: top genera get fixed colours, Other gets grey
  genera_present <- levels(df$genus_plot)
  n_genera       <- length(genera_present[genera_present != "Other"])
  pal_cols       <- BASE_COLS[seq_len(n_genera)]
  pal            <- setNames(pal_cols, genera_present[genera_present != "Other"])
  if ("Other" %in% genera_present) pal["Other"] <- "#a9a9a9"
  
  # Fix sample order - only keep samples present in data
  samples_present <- intersect(sample_order, unique(df$Sample))
  df$Sample <- factor(df$Sample, levels = samples_present)
  
  # Sample type annotation for x-axis colour
  sample_type <- ifelse(grepl("^A_", samples_present), "Aerosol", "Paired-env")
  x_colours   <- ifelse(sample_type == "Aerosol", "grey10", "grey50")
  
  p <- ggplot(df, aes(x = Sample, y = Abundance, fill = genus_plot)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = pal, name = "Genus") +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0.005, 0.005)
    ) +
    # Vertical separator between aerosol and paired-env
    geom_vline(xintercept = 5.5, linetype = "dashed",
               linewidth = 0.4, color = "grey50") +
    labs(
      title = domain_label,
      x     = NULL,
      y     = "Relative abundance"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 10, color = "grey10"),
      axis.text.x      = if (show_x) {
        element_text(angle = 45, hjust = 1, size = 8, color = x_colours)
      } else {
        element_blank()
      },
      axis.ticks.x     = if (show_x) element_line() else element_blank(),
      axis.text.y      = element_text(size = 8),
      axis.title.y     = element_text(size = 9),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "grey70"),
      legend.text      = element_text(face = "italic", size = 7.5),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.key.size  = unit(0.55, "lines"),
      plot.background  = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  return(p)
}

# =========================================================
# BUILD AND SAVE EACH PANEL INDIVIDUALLY
# =========================================================

# All panels show x-axis labels
p_bac <- make_panel(all_data$Bacteria$df,  all_data$Bacteria$top_genera,  "Bacteria",  show_x = TRUE)
p_arc <- make_panel(all_data$Archaea$df,   all_data$Archaea$top_genera,   "Archaea",   show_x = TRUE)
p_euk <- make_panel(all_data$Eukaryota$df, all_data$Eukaryota$top_genera, "Eukaryota", show_x = TRUE)
p_vir <- make_panel(all_data$Viruses$df,   all_data$Viruses$top_genera,   "Viruses",   show_x = TRUE)

# Preview each in RStudio
print(p_bac)
print(p_arc)
print(p_euk)
print(p_vir)

# Save individually
pdf("Supp_Taxonomy_Bacteria.pdf",  width = 6, height = 3); print(p_bac); dev.off()
pdf("Supp_Taxonomy_Archaea.pdf",   width = 6, height = 3); print(p_arc); dev.off()
pdf("Supp_Taxonomy_Eukaryota.pdf", width = 6, height = 3); print(p_euk); dev.off()
pdf("Supp_Taxonomy_Viruses.pdf",   width = 6, height = 3); print(p_vir); dev.off()

cat("Saved: 4 individual PDF files.\n")

# ------------------------------
# Step 4: Assign color palette
# ------------------------------
cbPalette <- c(
  '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
  '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
  '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'  # 20 colors total
)

cbPalette <- c(
  '#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0',
  '#24796C', '#DAA51B', '#2F8AC4', '#764E9F', '#ED645A',
  '#CC3A8E', '#A5AA99', '#9630EE', '#B79B6E', '#5B1A18'  # 15 colors total
)


# ============================================================
# Methanotroph genera — relative abundance (whole community)
# Input  : phyloseq object (non-rarefied or rarefied)
# Output : PDF stacked bar plot faceted by Sample_type
# ============================================================
# Load your phyloseq object
project_data <- readRDS("psf_bacteria_nonrarefied.rds")

# Define methanotroph genera (update if needed)
methanotroph_genera <- c(
  # Alpha (Type II)
  "Methylocystis","Methylosinus","Methylocella","Methyloferula","Methylocapsa",
  # Gamma / Methylococcales (Type I & X)
  "Methylomonas", "Methylomicrobium","Methylosoma","Methylosphaera",
  "Methylococcus","Methylocaldum","Methylothermus","Methylovulum","Methyloprofundus",
  "Methylomarinum","Crenothrix",
  # Post-2018 splits (be safe and include)
  "Methylotuvimicrobium","Methylotuvicoccus","Methylotetracoccus","Methylotuvibacter",
  # Verrucomicrobiota
  "Methylacidiphilum","Methylacidimicrobium",
  # NC10 (bacterial, nitrite-dependent AOM)
  "Candidatus Methylomirabilis"
)


# ---- 1) Normalize whole community, then subset methanotrophs ----
ps_rel_all  <- transform_sample_counts(project_data, function(x) x / sum(x))
ps_meth_rel <- subset_taxa(ps_rel_all, genus %in% methanotroph_genera)

# If none present, stop early with a clear message
if (ntaxa(ps_meth_rel) == 0) {
  stop("No methanotroph genera found in this object after filtering.")
}

# ---- 2) Melt safely & tidy labels ----
df_meth <- psmelt(ps_meth_rel)

# psmelt may have renamed sample column; restore as 'Sample' if needed
if ("sample_Sample" %in% names(df_meth)) df_meth$Sample <- df_meth$sample_Sample

# Make genus a clean character; fill blanks as 'Unclassified'
df_meth$genus <- as.character(df_meth$genus)
df_meth$genus[df_meth$genus == "" | is.na(df_meth$genus)] <- "Unclassified"


# Map colors only to the genera that are present (extend if needed)
g_levels <- sort(unique(df_meth$genus))
if (length(base_cols) < length(g_levels)) {
  base_cols <- colorRampPalette(base_cols)(length(g_levels))
} else {
  base_cols <- base_cols[seq_along(g_levels)]
}
pal <- setNames(base_cols, g_levels)

library(dplyr)
df_meth_plot <- df_meth %>%
  group_by(Sample, Sample_type, genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# ---- 4) Plot & save ----
pdf("Methanotroph_genera_by_sampletype.pdf", width = 12, height = 8)
ggplot(df_meth_plot, aes(x = Sample, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity", width = 0.9, colour = NA) +
  facet_grid(~ Sample_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = "Genus") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),  # show tenths of a percent
    expand = expansion(mult = c(0, 0.02)),           # start at 0, tiny headroom
    breaks = scales::breaks_extended(n = 6)
  ) +
  labs(
    title = "Methanotroph genera across samples (normalized to whole community)",
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    strip.background = element_blank(),
    strip.text = element_text(color = 'black')
  ) +
  guides(fill = guide_legend(ncol = 1))
dev.off()

pdf("Methanotroph_genera_by_heatmap.pdf", width = 12, height = 8)
# Heatmap of methanotroph genera (relative abundance, whole community)
ggplot(df_meth, aes(x = Sample, y = genus, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(~ Sample_type, scales = "free_x", space = "free_x") +
  # perceptually-uniform palette; sqrt handles zeros and wide dynamic range
  scale_fill_viridis_c(
    option = "magma",
    trans   = "sqrt",
    breaks  = c(0, 0.0001, 0.0005, 0.001),   # 0.1% and 0.5% (fractions, since your data is 0–1 relative abundance)
    labels  = scales::percent_format(accuracy = 0.01),
    limits  = c(0, 0.001),       # cap at 0.5% for clarity
    na.value = "white",
    guide   = guide_colorbar(barheight = unit(4, "cm"))
  ) +
  labs(
    title = "Methanotroph genera across samples (relative to whole community)",
    x = "Sample",
    y = "Genus",
    fill = "Rel. abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid        = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(color = "black"),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y       = element_text(size = 14),
    panel.spacing.x   = unit(0.6, "lines")
  )
dev.off()
# NMDS ordination and plot setup for rarefied genus-level metagenomic data
# --------------------------------------------------------------
project_data <- readRDS("psf_bacteria_raref_min.rds")
project_data <- readRDS("psf_archaea_raref_min.rds")
project_data <- readRDS("psf_eukaryota_raref_min.rds")
project_data <- readRDS("psf_viruses_raref_min.rds")


# 1. Run NMDS on Bray-Curtis distance

set.seed(6)
NMDS.bray <- ordinate(
  physeq   = project_data,
  method   = "NMDS",
  distance = "bray",
  trymax   = 200,        # allow more tries for a stable solution
  autotransform = FALSE  # counts already rarefied; no extra transform
)

# Quick sanity: print stress (lower is better; <0.2 ok, <0.1 great)
NMDS.bray$stress

# Extract sample metadata as data.frame
NMDS_metadata <- as(sample_data(project_data), "data.frame")

# Extract NMDS coordinates
NMDS_coords <- as.data.frame(NMDS.bray$points)

# Make sure rownames line up
stopifnot(all(rownames(NMDS_coords) == rownames(NMDS_metadata)))

# Merge coordinates into metadata
NMDS_metadata$NMDS1 <- NMDS_coords$MDS1
NMDS_metadata$NMDS2 <- NMDS_coords$MDS2

head(NMDS_metadata)


# Example using "Site" for shape, "month" for color
pdf("NMDS_virus.pdf", width = 8, height = 8)

ggplot(NMDS_metadata, aes(x = NMDS1, y = NMDS2, color = Sample_type2, shape=Site)) +
  geom_point(size = 9, alpha = 0.75) +
  #facet_wrap(~ time, scales = "free") +
  theme_bw() +
  #scale_shape_manual(values = c(15, 16, 17, 18, 8)) +  # adjust as needed
  scale_color_manual(values = c("Air" = "#009E73", "Reference" = "#7A7A7A")) +            # use same color palette
  labs(title = "NMDS Ordination (Bray-Curtis)",
       x = "NMDS1", y = "NMDS2") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(color = 'black')
  ) +
  guides(
    color = guide_legend(ncol = 1),
    shape = guide_legend(ncol = 1)
  )

dev.off()




project_data <- readRDS("psf_bacteria_raref_min.rds")
project_data <- readRDS("psf_archaea_raref_min.rds")
project_data <- readRDS("psf_eukaryota_raref_min.rds")
project_data <- readRDS("psf_viruses_raref_min.rds")

## Collapse to genus (removes NA genera)
project_data <- tax_glom(project_data, taxrank = "genus", NArm = TRUE)

### Alpha Diversity with Richness and Evenness
# --------------------------------------------
# Step 1: Estimate alpha diversity
# --------------------------------------------
alpha_div <- estimate_richness(project_data, measures = c("Observed", "Shannon"))
alpha_div$Evenness <- with(alpha_div, Shannon / log(Observed))
alpha_div$Sample <- rownames(alpha_div)

# --------------------------------------------
# Step 2: Merge with metadata
# --------------------------------------------
metadata_df <- as(sample_data(project_data), "data.frame")
metadata_df$Sample <- rownames(metadata_df)
alpha_df <- left_join(metadata_df, alpha_div, by = "Sample")

# Simple t-test between 2 groups
ttest_shan <- t.test(Shannon ~ Sample_type2, data = alpha_df)
print(ttest_shan)
cat("t =", ttest_shan$statistic, "\n")
cat("p-value =", ttest_shan$p.value, "\n")

# --------------------------------------------
# Step 3: Choose group variable (e.g., Site, virus_dose, Media)
# --------------------------------------------
library(agricolae)
group_var <- "Sample_type2"  # Change to your variable of interest

# --------------------------------------------
# Step 4: ANOVA + Tukey HSD for Richness
# --------------------------------------------
anova_rich <- aov(Observed ~ Sample_type2, data = alpha_df)
tukey_rich <- HSD.test(anova_rich, trt = "Sample_type2", group = TRUE)
alpha_df$rich_group <- tukey_rich$groups[match(alpha_df$Sample_type2, rownames(tukey_rich$groups)), "groups"]

# --------------------------------------------
# Step 5: ANOVA + Tukey HSD for Shannon
# --------------------------------------------
anova_shan <- aov(Shannon ~ Sample_type2, data = alpha_df)
tukey_shan <- HSD.test(anova_shan, trt = "Sample_type2", group = TRUE)
alpha_df$shan_group <- tukey_even$groups[match(alpha_df$Sample_type2, rownames(tukey_even$groups)), "groups"]



# --------------------------------------------
# Step 2: Plot richness and evenness
# --------------------------------------------
pdf("Shannon_virus.pdf", width = 3, height = 5)
ggplot(alpha_df, aes(x = Sample_type2, y = Shannon, fill = Sample_type2)) +
  # Thin box
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5, color = "black") +
  # Jittered points (important)
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.9, shape = 21, color = "black") +
  # Clean theme
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black",
                                fill = NA,
                                linewidth = 0.8)
  ) +
  # Remove legend (since color obvious)
  guides(fill = "none") +
  # Axis labels
  labs(x = NULL, y = "Shannon diversity") +
  # Manual colors (consistent with NMDS)
  scale_fill_manual(values = c("Air" = "#1b9e77",
                               "Reference" = "grey60")) +
  # Add p-value annotation (adjust y position manually)
  annotate("text", x = 1.5, y = max(alpha_df$Shannon) * 1.05,
           label = "t-test P = 0.9943068",
           size = 5)
dev.off()

# Example 2: Evenness 
pdf("Evenness_virus.pdf", width = 2, height = 6)
ggplot(alpha_df, aes(x = Sample_type2, y = Evenness)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(width = 0.2, size = 2) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Alpha Diversity (Pielou's Evenness)", x = "Sample_type", y = "Evenness") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#E58606", "#5D69B1", "#52BCA3", "#764E9F", "#CC61B0"))
dev.off()

# --------------------------------------------
# Step 6: Richness plot with Tukey letters
# --------------------------------------------
p_rich <- ggplot(alpha_df, aes(x = Sample_type2, y = Chao1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = max, geom = "text", aes(label = rich_group), vjust = -0.5) +
  theme_bw() +
  #facet_wrap(~ Site, ncol=5, scales = "fixed") +
  labs(title = "Alpha Diversity (Chao1 Richness by Site)", y = "Chao1", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cbPalette[1:length(unique(alpha_df$Site))])

print(p_rich)
# --------------------------------------------
p_even <- ggplot(alpha_df, aes(x = Sample_type2, y = Evenness)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(width = 0.2, size = 2) +
  stat_summary(fun = max, geom = "text", aes(label = rich_group), vjust = -0.5) +
  theme_bw() +
  facet_wrap(~ Site, ncol=5, scales = "fixed") +
  labs(title = "Alpha Diversity (Evenness by Site)", y = "Evenness", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cbPalette[1:length(unique(alpha_df$Site))])

print(p_even)
# --------------------------------------------
# Step 8: Save plots
# --------------------------------------------
ggsave("AlphaDiversity_Chao1_Tukey.pdf", p_rich, width = 10, height = 6)
ggsave("AlphaDiversity_Evenness_Tukey.pdf", p_even, width = 10, height = 6)



project_data <- readRDS("psf_bacteria_raref_min.rds")
project_data <- readRDS("psf_archaea_raref_min.rds")
project_data <- readRDS("psf_eukaryota_raref_min.rds")
project_data <- readRDS("psf_viruses_raref_min.rds")

##PERMANOVA TEST##

#subsampling for targetting comparisons and match with NMDS results

#project_data <- project_data %>%
  #subset_samples(sample_type %in% c("kelp_seedstring", "kelp")) %>%
  #subset_samples(kelp_type != "outplanted")

project_bray <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))


#adonis(project_bray ~ Type, data=sample_df, method="bray")
adonis2(project_bray ~ Sample_type2,
        data= sample_df , permutations=9999, by = "margin")


## betadispersion test
beta.Type <- betadisper(project_bray, sample_df$Sample_type2) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
permutest(beta.Type) #do permutation test of above #this is more customizable, see documentation

####IndVal ANALYSIS####

library(phyloseq)
library(indicspecies)  # multipatt()
library(vegan)         # how()
library(dplyr)

## Pick one domain to start (example: bacteria, rarefied at min-depth)
project_data <- readRDS("psf_bacteria_nonrarefied.rds")
project_data <- readRDS("psf_archaea_nonrarefied.rds")
project_data <- readRDS("psf_eukaryota_nonrarefied.rds")
project_data <- readRDS("psf_viruses_nonrarefied.rds")


# Pull metadata and build grouping variable
meta <- as(sample_data(project_data), "data.frame")
stopifnot("Sample_type2" %in% names(meta))   # expects values like "Air" and "Soil"/"Reference"

grp <- factor(ifelse(meta$Sample_type2 == "Air", "Air", "Reference"),
              levels = c("Air","Reference"))

# sanity check
table(grp)

## Step 2 — Collapse to genus and prepare matrix for IndVal
ps_gen <- tax_glom(project_data, taxrank = "genus", NArm = TRUE)

# Build counts matrix with samples as rows, genera as columns
mat <- as(otu_table(ps_gen), "matrix")
if (!taxa_are_rows(ps_gen)) mat <- t(mat)
X <- t(mat)  # samples in rows

# Optional: tiny prevalence filter (present in ≥ 2 samples)
prev <- colSums(X > 0)
keep <- prev >= 1
Xf <- X[, keep, drop = FALSE]

cat("Genera kept:", ncol(Xf), "of", ncol(X), "\n")

## Step 3 — IndVal (genus level)
library(indicspecies)
library(vegan)  # includes permute dependency

set.seed(42)


#  Build control object with stratification—this is valid in current CRAN versions:
ctrl <- how(nperm = 999, blocks = meta$Site)

set.seed(42)
res <- multipatt(Xf, grp, func = "IndVal.g",
                 control = ctrl)

summary(res, indvalcomp = TRUE, alpha = 0.05)

# Also save a human-readable summary to file
txt <- capture.output(summary(res, indvalcomp = TRUE, alpha = 0.05))
writeLines(txt, "indval_genus_bacteria_Air_vs_Reference_summary.txt")

## ============================================================
## Heat map of significant Air indicator genera (IndVal, p < 0.05)
## Requires:
##   - project_data  : phyloseq object used for IndVal
##   - res           : result from indicspecies::multipatt on *genus*-level matrix
## Output:
##   - INDVAL_Air_indicator_genera_heatmap.pdf
## ============================================================


## ---- 0) Get Air-only indicators (single-group, p < 0.05) from res$sign ----
tab <- as.data.frame(res$sign)
tab$taxon <- rownames(tab)

# sanity: expect membership columns named s.Air and s.Reference
stopifnot(all(c("s.Air","s.Reference","stat","p.value") %in% colnames(tab)))

air_map <- tab %>%
  filter(!is.na(p.value),
         `s.Air` == 1, `s.Reference` == 0,         # Air-only
         p.value < 0.05) %>%                        # significant only
  select(taxon, stat, p.value)

message("Air indicator TAXA (rows of res$sign, p<0.05): ", n_distinct(air_map$taxon))

## ---- 1) Collapse to genus & create a taxon->genus lookup (clean 'Candidatus ') ----
ps_gen <- tax_glom(project_data, taxrank = "genus", NArm = TRUE)

lookup <- data.frame(
  taxon     = taxa_names(ps_gen),
  genus_raw = as.character(tax_table(ps_gen)[,"genus"]),
  stringsAsFactors = FALSE
)

air_taxa_to_genus <- air_map %>%
  left_join(lookup, by = "taxon") %>%
  mutate(genus_raw   = ifelse(is.na(genus_raw) | genus_raw == "", "Unclassified", genus_raw),
         genus_clean = sub("^Candidatus\\s+", "", genus_raw))

air_genus_set <- unique(air_taxa_to_genus$genus_clean)
message("Air indicator GENERA (unique, after cleaning labels): ", length(air_genus_set))

## ---- 2) Relative abundances for *all* genera, then keep only Air indicators ----
ps_gen_rel <- transform_sample_counts(ps_gen, function(x) x / sum(x))

df <- psmelt(ps_gen_rel)
if ("sample_Sample" %in% names(df)) df$Sample <- df$sample_Sample  # phyloseq quirk

df <- df %>%
  mutate(genus = as.character(genus),
         genus = ifelse(is.na(genus) | genus == "", "Unclassified", genus),
         genus = sub("^Candidatus\\s+", "", genus)) %>%             # match cleaned set
  filter(genus %in% air_genus_set) %>%                              # keep Air indicators only
  select(Sample, genus, Abundance)

## ---- 3) Attach facet column: prefer Sample_type, fallback to Sample_type2 ----
meta_ps <- as(sample_data(ps_gen_rel), "data.frame")
meta_ps$Sample <- rownames(meta_ps)

facet_col <- if ("Sample_type" %in% names(meta_ps)) "Sample_type" else "Sample_type2"
stopifnot(facet_col %in% names(meta_ps))

sample_map <- meta_ps %>%
  transmute(Sample, FacetType = .data[[facet_col]])

df <- df %>%
  left_join(sample_map, by = "Sample")

## ---- 4) Complete the Sample × Genus grid so zeros plot (dark purple), no NAs ----
df <- df %>%
  complete(Sample, genus, fill = list(Abundance = 0)) %>%
  left_join(sample_map, by = "Sample")

# normalize FacetType name in case a suffix appeared
if (!"FacetType" %in% names(df)) {
  ft_cols <- grep("^FacetType", names(df), value = TRUE)
  df$FacetType <- df[[ft_cols[length(ft_cols)]]]
}
stopifnot("FacetType" %in% names(df))

df$Abundance <- tidyr::replace_na(df$Abundance, 0)

## ---- 5) Drop genera that are zero across *all* samples (safety) ----
df <- df %>%
  group_by(genus) %>%
  filter(sum(Abundance, na.rm = TRUE) > 0) %>%
  ungroup()

## ---- 6) Order rows by total abundance (descending) ----
## compute totals and order genera
ord_tbl <- df %>%
  group_by(genus) %>%
  summarise(tot = sum(Abundance), .groups = "drop") %>%
  arrange(desc(tot))

# (optional) keep only top 10 genera
TOP_N <- 15   # change to Inf to keep all
if (is.finite(TOP_N)) {
  top_genera <- head(ord_tbl$genus, TOP_N)
  df <- df %>% filter(genus %in% top_genera)
  ord_tbl <- ord_tbl %>% filter(genus %in% top_genera)
}

df <- df %>%
  mutate(genus = factor(genus, levels = ord_tbl$genus))

## ---- 7) Make per-facet x order so each column facet shows its own samples nicely ----
df <- df %>%
  arrange(FacetType, Sample) %>%
  group_by(FacetType) %>%
  mutate(Sample_f = factor(Sample, levels = unique(Sample))) %>%
  ungroup()

## ---- 8) Plot (facets = column panels by FacetType; Air indicators only) ----
pdf("INDVAL_Other_indicator_genera_heatmap_bacteria.pdf", width = 5, height = 7)

ggplot(df, aes(x = Sample_f, y = genus, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(~ FacetType, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(
    option = "magma", trans = "sqrt",
    breaks = c(0, 1e-3, 5e-3, 1e-2),        # adjust to data range
    labels = scales::percent_format(accuracy = 0.01),
    limits = c(0, 0.005),                     # now covers up to 1%
    oob = scales::squish,
    guide = guide_colorbar(barheight = unit(4, "cm"))
  ) +
  labs(
    title = "Air indicator genera (IndVal; p < 0.05) — relative abundance",
    x = "Sample",
    y = "Genus",
    fill = "Rel. abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(color = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y      = element_text(size = 11),
    panel.spacing.x  = unit(0.6, "lines")
  )

dev.off()


project_data <- readRDS("psf_bacteria_raref_min.rds")



## ================================
## 0) Setup
## ================================
library(phyloseq)
library(DESeq2)
library(dplyr)
library(tidyr)
library(readr)

setwd("~/Desktop/Airmicrobiome_project/GSC_aeromicrobiome/R_work/")
## ---- choose your object (nonrarefied preferred) ----
# project_data <- readRDS("psf_bacteria_nonrarefied.rds")  # <- recommended
project_data <- readRDS("psf_bacteria_nonrarefied.rds")
project_data <- readRDS("psf_archaea_nonrarefied.rds")
project_data <- readRDS("psf_eukaryota_nonrarefied.rds")
project_data <- readRDS("psf_viruses_nonrarefied.rds")

## ---- grouping variable: Air vs Reference ----
meta <- as(sample_data(project_data), "data.frame")
grp_col <- if ("Sample_type2" %in% names(meta)) "Sample_type2" else "Sample_type"
stopifnot(grp_col %in% names(meta))
# Make a clean 2-level factor: Air vs Reference (everything not Air -> Reference)
meta$Group <- factor(ifelse(meta[[grp_col]] == "Air", "Air", "Reference"),
                     levels = c("Reference","Air"))
sample_data(project_data)$Group <- meta$Group

cat("Group counts:\n"); print(table(sample_data(project_data)$Group))

## ================================
## 1) Collapse to GENUS
## ================================
ps_gen <- tax_glom(project_data, taxrank = "genus", NArm = TRUE)

## Optional: drop all-zero genera (just in case)
ps_gen <- prune_taxa(taxa_sums(ps_gen) > 0, ps_gen)

## ================================
## 2) DESeq2 (Air vs Reference) at genus level
## ================================

# Build DESeqDataSet from the genus-collapsed phyloseq object
dds <- phyloseq_to_deseq2(ps_gen, ~ Group)   # counts are coerced to integer

# (A) Gentle prefilter: keep genera with ≥10 total counts across all samples
keep <- rowSums(counts(dds)) >= 15
dds  <- dds[keep, ]

# (B) Size-factor estimation that tolerates many zeros
dds <- estimateSizeFactors(dds, type = "poscounts")
summary(sizeFactors(dds))  # quick sanity check

# (C) Make sure "Reference" is the baseline level
dds$Group <- relevel(dds$Group, ref = "Reference")

# (D) Differential abundance (Wald test)
dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)

# (E) Get results: positive log2FC => enriched in Air
res <- results(dds, contrast = c("Group","Air","Reference"), alpha = 0.05)

# (F) LFC shrinkage (apeglm). Get the exact coef name then shrink.
resultsNames(dds)  # inspect; expect "Group_Air_vs_Reference"
sres <- lfcShrink(dds, coef = "Group_Air_vs_Reference", type = "apeglm")

# (G) Add genus names and tidy for export
tx <- as.data.frame(tax_table(ps_gen))
tx$taxon <- rownames(tx)

out_deseq <- as.data.frame(sres) |>
  tibble::rownames_to_column("taxon") |>
  dplyr::left_join(tx[, c("taxon","genus")], by = "taxon") |>
  dplyr::mutate(
    genus = dplyr::coalesce(as.character(genus), "Unclassified"),
    genus = sub("^Candidatus\\s+", "", genus)  # drop "Candidatus "
  ) |>
  dplyr::arrange(padj, dplyr::desc(log2FoldChange))

# Save to CSV
readr::write_csv(out_deseq, "DESeq2_virus_genus_Air_vs_Reference_all.csv")

################
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(patchwork)

# =========================================================
# SETTINGS
# =========================================================
PADJ_CUT  <- 0.05
LFC_CUT   <- 2       # conservative: 4x fold change
TOP_N     <- 15      # labels per panel

COL_AIR   <- "#3A6B8A"   # blue  - aerosol enriched
COL_ENV   <- "#C1703A"   # orange - env enriched
COL_NS    <- "grey82"    # not significant

# =========================================================
# HELPER: build tidy DESeq2 results from phyloseq object
# =========================================================

get_deseq_results <- function(rds_path) {
  
  ps   <- readRDS(rds_path)
  meta <- as(sample_data(ps), "data.frame")
  
  grp_col <- if ("Sample_type2" %in% names(meta)) "Sample_type2" else "Sample_type"
  
  meta$Compartment <- factor(
    ifelse(meta[[grp_col]] == "Air", "Air", "Reference"),
    levels = c("Reference", "Air")
  )
  meta$Site <- factor(sub("^[AWS]_([A-Z]+).*", "\\1", rownames(meta)))
  
  sample_data(ps)$Compartment <- meta$Compartment
  sample_data(ps)$Site        <- meta$Site
  
  # Collapse to genus
  ps_gen <- tax_glom(ps, taxrank = "genus", NArm = TRUE)
  ps_gen <- prune_taxa(taxa_sums(ps_gen) > 0, ps_gen)
  
  # DESeq2 site-blocked
  dds  <- phyloseq_to_deseq2(ps_gen, ~ Site + Compartment)
  dds  <- dds[rowSums(counts(dds)) >= 10, ]
  dds  <- estimateSizeFactors(dds, type = "poscounts")
  dds$Compartment <- relevel(dds$Compartment, ref = "Reference")
  dds  <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
  
  coef_name <- grep("Compartment_Air", resultsNames(dds), value = TRUE)
  sres <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # Add genus taxonomy
  tx <- as.data.frame(tax_table(ps_gen)) %>% rownames_to_column("taxon")
  
  out <- as.data.frame(sres) %>%
    rownames_to_column("taxon") %>%
    left_join(tx[, c("taxon", "genus")], by = "taxon") %>%
    mutate(
      genus = coalesce(as.character(genus), "Unclassified"),
      genus = sub("^Candidatus\\s+", "", genus),
      padj  = replace_na(padj, 1),
      log2FoldChange = replace_na(log2FoldChange, 0)
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  return(out)
}

# =========================================================
# HELPER: make volcano ggplot
# =========================================================

make_volcano <- function(df, title_label) {
  
  res <- df %>%
    mutate(
      sig  = padj < PADJ_CUT & abs(log2FoldChange) >= LFC_CUT,
      side = case_when(
        sig & log2FoldChange > 0 ~ "Enriched in aerosol",
        sig & log2FoldChange < 0 ~ "Enriched in paired-env",
        TRUE                     ~ "Not significant"
      ),
      side = factor(side, levels = c(
        "Enriched in aerosol", "Enriched in paired-env", "Not significant"
      ))
    )
  
  # Top N labels: top by padj within each direction
  label_df <- res %>%
    filter(sig) %>%
    arrange(padj) %>%
    slice_head(n = TOP_N)
  
  n_air <- sum(res$side == "Enriched in aerosol")
  n_env <- sum(res$side == "Enriched in paired-env")
  
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    
    # All points
    geom_point(
      aes(color = side),
      size  = 5,
      alpha = 0.6
    ) +
    
    # Labels for top hits only
    geom_text_repel(
      data               = label_df,
      aes(label          = genus),
      size               = 2.8,
      fontface           = "italic",
      color              = "grey15",
      max.overlaps       = Inf,
      min.segment.length = 0.2,
      segment.color      = "grey55",
      segment.size       = 0.3,
      box.padding        = 0.4,
      point.padding      = 0.2
    ) +
    
    # Threshold lines
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dashed", linewidth = 0.35, color = "grey50") +
    geom_hline(yintercept = -log10(PADJ_CUT),
               linetype = "dashed", linewidth = 0.35, color = "grey50") +
    
    scale_color_manual(
      values = c(
        "Enriched in aerosol"    = COL_AIR,
        "Enriched in paired-env" = COL_ENV,
        "Not significant"        = COL_NS
      ),
      name = NULL
    ) +
    
    labs(
      title    = title_label,
      subtitle = paste0(
        "Site-blocked DESeq2 (~ Site + Compartment)  |  ",
        "padj < ", PADJ_CUT, ", |LFC| >= ", LFC_CUT, "\n",
        n_air, " aerosol-enriched  ·  ", n_env, " paired-env-enriched"
      ),
      x = expression(log[2]~"fold change (Aerosol / Paired-env)"),
      y = expression(-log[10]~"(adjusted "*italic(p)*"-value)")
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 11, color = "grey10"),
      plot.subtitle      = element_text(size = 8, color = "grey45", lineheight = 1.3),
      axis.text          = element_text(size = 9),
      axis.title         = element_text(size = 10),
      legend.position    = "top",
      legend.text        = element_text(size = 9),
      legend.key.size    = unit(0.8, "lines"),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey70"),
      plot.background    = element_rect(fill = "white", color = NA)
    )
}

# =========================================================
# RUN BACTERIA
# =========================================================

cat("Running Bacteria...\n")
bac_df  <- get_deseq_results("psf_bacteria_nonrarefied.rds")
p_bac   <- make_volcano(bac_df, "Bacteria")
print(p_bac)   # preview in RStudio

# =========================================================
# RUN EUKARYOTA
# =========================================================

cat("Running Eukaryota...\n")
euk_df  <- get_deseq_results("psf_eukaryota_nonrarefied.rds")
p_euk   <- make_volcano(euk_df, "Eukaryota")
print(p_euk)   # preview in RStudio

# =========================================================
# COMBINED FIGURE
# =========================================================

combined <- p_bac + p_euk +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Differential genus abundance between aerosol and paired environmental samples",
    theme = theme(
      plot.title      = element_text(face = "bold", size = 12, color = "grey10"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

print(combined)   # preview in RStudio

# Save
pdf("DESeq2_genus_top_effects.pdf", width = 12, height = 6)
print(combined)
dev.off()

cat("Saved: DESeq2_genus_top_effects.pdf\n")

# Also save individually
pdf("DESeq2_bacteria_volcano.pdf", width = 5, height = 6)
print(p_bac)
dev.off()

pdf("DESeq2_eukaryota_volcano.pdf", width = 5, height = 6)
print(p_euk)
dev.off()

cat("Saved individual PDFs.\n")
# Set working directory
setwd("~/Desktop/Airmicrobiome_project/GSC_aeromicrobiome/R_work/gene_centric/")

cpml <- "methane_panel.cpm_long.tsv"
cpmm <- "methane_panel.cpm_matrix.tsv"

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(tibble)})

# wide is already tab-delimited — keep using read_tsv
wide <- read_tsv(cpmm, show_col_types = FALSE)

# long is whitespace-delimited — read with read_table()
long <- read_table(
  cpml,
  show_col_types = FALSE,
  col_types = cols(
    sample = col_character(),
    KO     = col_character(),
    gene   = col_character(),
    block  = col_character(),
    CPM    = col_double()
  )
)

# sanity checks
glimpse(long)
stopifnot(all(c("sample","KO","gene","block","CPM") %in% names(long)))
stopifnot(is.numeric(long$CPM), !any(is.na(long$CPM)))


## ---- heatmap_ko_by_sample ----
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(forcats)})

# 1) Set nice ordering for blocks (only keep ones actually present)
block_order <- c("AntiMarkers","FormateOx","H4F/one-C","H4MPT_detox","MDH","RuMP","Serine","pMMO","sMMO")
long <- long %>%
  mutate(block = factor(block, levels = intersect(block_order, unique(block))))

# 2) Order KOs within blocks and create a label "KO:gene" for the y-axis
long <- long %>%
  mutate(KO_label = paste0(KO, " : ", gene)) %>%
  arrange(block, KO, gene) %>%
  mutate(KO_label = factor(KO_label, levels = unique(KO_label)))

# 3) Order samples using the wide matrix column order
sample_order <- names(wide)[4:ncol(wide)]
long <- long %>% mutate(sample = factor(sample, levels = sample_order))

# 4) Log10 transform for heatmap
long <- long %>% mutate(log10_cpm = log10(CPM + 1))

# 5) Plot
p <- ggplot(long, aes(x = sample, y = KO_label, fill = log10_cpm)) +
  geom_tile() +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(name = "log10(CPM + 1)", low = "white", high = "steelblue") +
  labs(x = "Sample", y = "KO : gene", title = "Methane panel gene abundances (CPM)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(angle = 0)
  )

print(p)
ggsave("methane_panel.heatmap.pdf", p, width = 10, height = 9)

# Optional: quick peek of block totals per sample (table)
block_totals <- long %>% group_by(sample, block) %>% summarise(CPM = sum(CPM), .groups="drop")
print(block_totals %>% tidyr::pivot_wider(names_from = sample, values_from = CPM))

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ----------------------------------------
# Step 1: Load bacteria and eukaryota
# ----------------------------------------
bac_data <- readRDS("psf_bacteria_nonrarefied.rds")
euk_data <- readRDS("psf_eukaryota_nonrarefied.rds")

# ----------------------------------------
# Step 2: Collapse to genus
# ----------------------------------------
bac_data <- tax_glom(bac_data, taxrank = "genus", NArm = TRUE)
euk_data <- tax_glom(euk_data, taxrank = "genus", NArm = TRUE)

# ----------------------------------------
# Step 3: Make sure same samples in both
# ----------------------------------------
shared_samples <- intersect(sample_names(bac_data), sample_names(euk_data))
bac_data <- prune_samples(shared_samples, bac_data)
euk_data <- prune_samples(shared_samples, euk_data)
cat("Shared samples:", length(shared_samples), "\n")

# ----------------------------------------
# Step 4: Convert to relative abundance
# ----------------------------------------
bac_rel <- transform_sample_counts(bac_data, function(x) x / sum(x))
euk_rel <- transform_sample_counts(euk_data, function(x) x / sum(x))

# ----------------------------------------
# Step 5: Bray-Curtis distance
# ----------------------------------------
bac_bray <- phyloseq::distance(bac_rel, method = "bray")
euk_bray <- phyloseq::distance(euk_rel, method = "bray")

# ----------------------------------------
# Step 6: NMDS ordination
# ----------------------------------------
bac_nmds <- metaMDS(bac_bray, k = 2, trymax = 100)
euk_nmds <- metaMDS(euk_bray, k = 2, trymax = 100)

# ----------------------------------------
# Step 7: Procrustes + Protest
# ----------------------------------------
procrustes_result <- procrustes(bac_nmds, euk_nmds)
protest_result    <- protest(bac_nmds, euk_nmds, permutations = 999)
print(protest_result)

# ----------------------------------------
# Step 8: Extract coordinates for ggplot2
# ----------------------------------------
# Procrustes gives: X = reference (bacteria), Y = rotated (eukaryota)
proc_df <- data.frame(
  # Bacteria (reference) coordinates
  bac_x = procrustes_result$X[, 1],
  bac_y = procrustes_result$X[, 2],
  # Eukaryota (rotated) coordinates
  euk_x = procrustes_result$Yrot[, 1],
  euk_y = procrustes_result$Yrot[, 2],
  Sample = rownames(procrustes_result$X)
)

# Pull in metadata (Sample_type2)
meta_df <- as(sample_data(bac_data), "data.frame")
meta_df$Sample <- rownames(meta_df)
proc_df <- left_join(proc_df, meta_df[, c("Sample", "Sample_type2")], by = "Sample")

# ----------------------------------------
# Step 9: Build publication-quality ggplot
# ----------------------------------------

# Annotation text from protest
r_val <- round(protest_result$t0, 3)
p_val <- round(protest_result$signif, 3)
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", p_val))
annotation_text <- paste0("r = ", r_val, "\n", p_label)

# Colors — adjust to match your palette
sample_colors <- c("Air" = "#009E73", "Reference" = "#7A7A7A")

p <- ggplot(proc_df) +
  # Arrows from eukaryota to bacteria
  geom_segment(aes(x = euk_x, y = euk_y,
                   xend = bac_x, yend = bac_y,
                   color = Sample_type2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 0.7, alpha = 0.8) +
  # Eukaryota points (open circle)
  geom_point(aes(x = euk_x, y = euk_y, color = Sample_type2),
             shape = 1, size = 4, stroke = 1.2) +
  # Bacteria points (filled circle)
  geom_point(aes(x = bac_x, y = bac_y, color = Sample_type2),
             shape = 16, size = 4) +
  # Sample labels
  geom_text_repel(aes(x = bac_x, y = bac_y, label = Sample),
                  size = 2.8, max.overlaps = 20,
                  color = "grey30", segment.color = "grey70") +
  # Color scale
  scale_color_manual(values = sample_colors, name = "Sample Type") +
  # Annotation box (r and p value)
  annotate("text", x = Inf, y = Inf,
           label = annotation_text,
           hjust = 1.1, vjust = 1.5,
           size = 4, fontface = "italic", color = "grey20") +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  labs(
    title = "Procrustes Analysis: Bacteria vs Eukaryota",
    subtitle = "Filled = Bacteria  |  Open = Eukaryota  |  Arrows show residual distance",
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )

print(p)

# ----------------------------------------
# Step 10: Save high resolution
# ----------------------------------------
ggsave("Procrustes_Bac_vs_Euk_publication.pdf", plot = p,
       width = 5, height = 8, dpi = 400)
