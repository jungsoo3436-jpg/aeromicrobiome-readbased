###############################################################################
# Aerosol metagenome manuscript (read-based community analysis)
# Step 05: Genus-level composition across all four domains (Supplementary Fig. 4)
#
# Author: Jungsoo Park (Simon Fraser University)
#
# One stacked barplot panel per domain, showing the 15 most abundant genera
# across the dataset with everything else collapsed to "Other" and renormalised
# within each sample. This is the panel the Results cite when contrasting the
# broadly dispersed, low-abundance bacterial aerosol communities with the
# eukaryotic aerosols enriched in spore-forming and filamentous fungi.
#
# Input  : data/phyloseq/psf_{bacteria,archaea,eukaryota,viruses}_raref_min.rds
# Output : figures/SuppFig4_genus_composition_{Bacteria,Archaea,Eukaryota,Viruses}.pdf
#
# Samples are ordered aerosol first, then paired environmental, grouped by site
# with sites running from high to low ambient CH4. The dashed line separates the
# five aerosol samples from the ten paired environmental samples.
#
# Archaeal and viral panels rest on shallow rarefaction depths (965 and 511
# reads; see 01_build_phyloseq.R) and are shown for completeness only.
#
# Run from a FRESH R session, from the repository root.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(patchwork)
})

PS_DIR  <- "data/phyloseq"
OUT_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# SETTINGS
# =========================================================
N_TOP     <- 15
RDS_FILES <- list(
  Bacteria  = file.path(PS_DIR, "psf_bacteria_raref_min.rds"),
  Archaea   = file.path(PS_DIR, "psf_archaea_raref_min.rds"),
  Eukaryota = file.path(PS_DIR, "psf_eukaryota_raref_min.rds"),
  Viruses   = file.path(PS_DIR, "psf_viruses_raref_min.rds")
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

# Save individually
pdf(file.path(OUT_DIR, "SuppFig4_genus_composition_Bacteria.pdf"),  width = 6, height = 3); print(p_bac); dev.off()
pdf(file.path(OUT_DIR, "SuppFig4_genus_composition_Archaea.pdf"),   width = 6, height = 3); print(p_arc); dev.off()
pdf(file.path(OUT_DIR, "SuppFig4_genus_composition_Eukaryota.pdf"), width = 6, height = 3); print(p_euk); dev.off()
pdf(file.path(OUT_DIR, "SuppFig4_genus_composition_Viruses.pdf"),   width = 6, height = 3); print(p_vir); dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo_05.txt"))
cat("\nDone. Supplementary Fig. 4 panels written to", OUT_DIR, "\n")
