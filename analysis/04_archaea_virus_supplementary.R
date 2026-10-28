###############################################################################
# Aerosol metagenome manuscript (read-based community analysis)
# Step 04: Archaeal and viral diversity and composition (Supplementary Fig. 3)
#
# Author: Jungsoo Park (Simon Fraser University)
#
# Archaea and Viruses showed no significant aerosol vs paired-environment
# difference in either diversity or composition and were not examined further
# in the paper. This script runs that test for both domains and writes the
# supporting panels.
#
# The rarefaction depths for these two domains are low (Archaea 965 reads,
# Viruses 511 reads; see 01_build_phyloseq.R), so these estimates are
# interpreted cautiously. That caveat is the reason these panels are
# supplementary rather than main-text.
#
# Input  : data/phyloseq/psf_{archaea,viruses}_raref_min.rds
# Output : figures/SuppFig3_archaea_virus.pdf
#          figures/SuppFig3_archaea_virus_stats.tsv
#
# Settings match 02_community_and_methane_figures.R so the archaeal and viral
# tests are directly comparable with the bacterial and eukaryotic ones:
# genus level, rarefied counts, Bray-Curtis, NMDS seed 6, PERMANOVA with
# 9,999 permutations, dispersion by betadisper.
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

site_levels <- c("VL", "DF", "MK", "BB", "RM")   # high -> low ambient CH4
col_type    <- c(Air = "#009E73", Reference = "#7A7A7A")

N_PERM   <- 9999
NMDS_SEED <- 6

# Pull the aerosol vs reference grouping and the site out of the sample names,
# the same way 02_community_and_methane_figures.R does.
add_groups <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")
  grp  <- if ("Sample_type2" %in% names(meta)) "Sample_type2" else "Sample_type"
  sample_data(ps)$Compartment <- factor(
    ifelse(meta[[grp]] == "Air", "Air", "Reference"),
    levels = c("Reference", "Air")
  )
  sample_data(ps)$Site <- factor(
    sub("^[AWS]_([A-Z]+).*", "\\1", rownames(meta)), levels = site_levels
  )
  ps
}

analyse_domain <- function(rds, domain_label) {

  ps <- add_groups(readRDS(file.path(PS_DIR, rds)))
  ps <- tax_glom(ps, taxrank = "genus", NArm = TRUE)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  meta <- as(sample_data(ps), "data.frame")

  # ---- Shannon diversity ---------------------------------------------------
  div <- estimate_richness(ps, measures = "Shannon") %>%
    rownames_to_column("Sample") %>%
    left_join(meta %>% rownames_to_column("Sample"), by = "Sample")

  tt <- t.test(Shannon ~ Compartment, data = div)

  p_div <- ggplot(div, aes(Compartment, Shannon, fill = Compartment)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(aes(shape = Site), width = 0.14, size = 3.2, alpha = 0.9) +
    scale_fill_manual(values = col_type, guide = "none") +
    labs(title = domain_label, x = NULL, y = "Shannon diversity",
         subtitle = sprintf("two-sample t-test, p = %.3f", tt$p.value)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  # ---- NMDS on Bray-Curtis ------------------------------------------------
  set.seed(NMDS_SEED)
  ord <- ordinate(ps, method = "NMDS", distance = "bray",
                  trymax = 200, autotransform = FALSE)

  crd <- as.data.frame(ord$points)
  stopifnot(all(rownames(crd) == rownames(meta)))
  crd$Compartment <- meta$Compartment
  crd$Site        <- meta$Site

  dm  <- phyloseq::distance(ps, method = "bray")
  set.seed(NMDS_SEED)
  pmv <- vegan::adonis2(dm ~ Compartment, data = meta, permutations = N_PERM)
  bd  <- vegan::permutest(vegan::betadisper(dm, meta$Compartment),
                          permutations = N_PERM)

  p_ord <- ggplot(crd, aes(MDS1, MDS2, colour = Compartment, shape = Site)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_colour_manual(values = col_type, name = NULL) +
    labs(title = domain_label, x = "NMDS1", y = "NMDS2",
         subtitle = sprintf("stress = %.3f | PERMANOVA R2 = %.3f, p = %.4f",
                            ord$stress, pmv$R2[1], pmv$`Pr(>F)`[1])) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  stats <- tibble(
    domain          = domain_label,
    n_genera        = ntaxa(ps),
    shannon_p       = tt$p.value,
    nmds_stress     = ord$stress,
    permanova_F     = pmv$F[1],
    permanova_R2    = pmv$R2[1],
    permanova_p     = pmv$`Pr(>F)`[1],
    betadisper_F    = bd$statistic,
    betadisper_p    = bd$tab[1, "Pr(>F)"]
  )

  list(div = p_div, ord = p_ord, stats = stats)
}

arc <- analyse_domain("psf_archaea_raref_min.rds", "Archaea")
vir <- analyse_domain("psf_viruses_raref_min.rds", "Viruses")

stats <- bind_rows(arc$stats, vir$stats)
print(stats)
write_tsv(stats, file.path(OUT_DIR, "SuppFig3_archaea_virus_stats.tsv"))

supp2 <- (arc$div + arc$ord) / (vir$div + vir$ord) +
  plot_layout(widths = c(1, 1.5), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Archaeal and viral diversity and composition in aerosol and paired environmental samples",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

ggsave(file.path(OUT_DIR, "SuppFig3_archaea_virus.pdf"),
       plot = supp2, width = 11, height = 8)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo_04.txt"))
cat("\nDone. Supplementary Fig. 3 panels and statistics written to", OUT_DIR, "\n")
