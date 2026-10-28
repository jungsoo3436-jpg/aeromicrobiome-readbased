###############################################################################
# Aerosol microbiome manuscript: analyses and figures (revision, 2026)
# Author: Jungsoo Park (SFU)
#
# Sections
#   0. Setup: paths, packages, shared settings
#   1. Supplementary Fig. 1: aerosol cell counts vs DNA yield and sequencing depth
#   2. Figure 2: flow cytometry (TCC, ICC, intact fraction)
#   3. Figure 3: community diversity/composition (Shannon, NMDS, Procrustes) + stats
#   4. Figure 4 numbers: site-blocked DESeq2 summary counts
#   5. Supplementary Table S2: domain read proportions
#   6. Figure 5A: C1 metabolism marker gene heatmap
#   7. Supplementary Fig. 2: domain composition per sample
#   8. Supplementary Table S3: differentially abundant genera
#   9. Figure 5B: pmoCAB log2 fold change vs ambient CH4
#  10. Figure 4: DESeq2 volcano plots (bacteria, eukaryota)
#  11. Supplementary Fig. 5: methanotroph genera in read-based profiles
#  12. Session info (record package versions for reproducibility)
#
# Sections 10 and 11 were originally run from a separate exploratory script and
# are consolidated here so that one script reproduces every read-based figure.
#
# Notes
# - Run from a FRESH R session (quit RStudio fully, do not restore .RData).
# - vegan and DESeq2 functions are called with "pkg::" to avoid masking conflicts.
# - Edit only the paths in Section 0 to run on another machine.
###############################################################################


# =============================================================================
# 0. SETUP
# =============================================================================

# ---- Paths (relative to the repository root) ---------------------------------
R_WORK  <- "data"                    # Bracken genus table and sample metadata
PS_DIR  <- "data/phyloseq"           # phyloseq objects from 01_build_phyloseq.R
KO_DIR  <- "data/gene_catalog"       # KO x sample proportional TPM matrix
OUT_DIR <- "figures"                 # all outputs go here
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Packages ----------------------------------------------------------------
# Install once if needed:
# install.packages(c("tidyverse", "ggrepel", "patchwork", "ggh4x", "vegan"))
# BiocManager::install(c("phyloseq", "DESeq2", "apeglm"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(ggh4x)
  library(phyloseq)
})

# ---- Shared settings ---------------------------------------------------------
site_levels <- c("VL", "DF", "MK", "BB", "RM")                       # ordered high -> low CH4
site_cols   <- c(VL = "#D55E00", DF = "#E69F00", MK = "#009E73",     # Okabe-Ito, colour-blind safe
                 BB = "#0072B2", RM = "#CC79A7")

fmt_p    <- function(p) if (p < 0.001) sprintf("%.4f", p) else sprintf("%.3f", p)  # no 2e-04 style
save_pdf <- function(p, file, w, h) { pdf(file.path(OUT_DIR, file), width = w, height = h); print(p); dev.off() }



# =============================================================================
# SUPPLEMENTARY FIG 1: aerosol cell counts vs DNA yield and sequencing depth
#   - Three panels in a 2 x 2 grid: TCC vs DNA, ICC vs DNA (top), ICC vs reads
#     (TCC vs reads is shown separately as Fig. 2C)
#   - SFU and CH: aerosol library failed -> no reads; DNA panels only
#   - SFU DNA below Qubit detection -> plotted at 0.1 ng (Spearman ranks unaffected)
#   - Requires Section 0 (packages, OUT_DIR)
# =============================================================================

supp1_df <- tibble(
  Site   = c("BB", "MK", "VL", "RM", "DF", "CH", "SFU"),
  TCC    = c(56003.39, 6134.244, 247268.2, 2725.61, 1507207, 5158.39, 1694.56),
  ICC    = c(1431.894, 4090.109, 139702.5, 38.11605, 84575.77, 1066.91, 302.24),
  dna    = c(13.77, 6.75, 17.64, 40.45, 74.1, 4.83, 0.1),
  reads  = c(6430565, 2630772, 10094145, 1866743, 98902286, NA, NA),
  Status = c(rep("Library recovered", 5), rep("Library not recovered", 2))
)

supp1_long <- supp1_df %>%
  pivot_longer(c(TCC, ICC), names_to = "Metric", values_to = "Cells") %>%
  pivot_longer(c(dna, reads), names_to = "Predictor", values_to = "X") %>%
  filter(!is.na(X), paste(Metric, Predictor) != "TCC reads") %>%   # TCC vs reads = Fig. 2C
  mutate(logCells = log10(Cells),
         logX     = log10(X),
         Panel    = factor(paste(Metric, Predictor),
                           levels = c("TCC dna", "ICC dna", "ICC reads"),
                           labels = c("TCC vs DNA yield", "ICC vs DNA yield",
                                      "ICC vs sequencing reads")))

supp1_stats <- supp1_long %>%
  group_by(Panel) %>%
  summarise(n     = n(),
            rho   = cor.test(Cells, X, method = "spearman", exact = TRUE)$estimate,
            p     = cor.test(Cells, X, method = "spearman", exact = TRUE)$p.value,
            x_pos = min(logX) + 0.06 * (max(logX) - min(logX)),
            y_pos = Inf,
            .groups = "drop") %>%
  mutate(label = sprintf("rho == %.2f*','~~italic(p) == %.3f*','~~n == %d", rho, p, n))

print(supp1_stats)

supp1 <- ggplot(supp1_long, aes(logX, logCells)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2,
              colour = scales::alpha("blue", 0.5)) +              # visual guide only
  geom_point(aes(shape = Status), size = 6, alpha = 0.7, stroke = 1.4) +
  scale_shape_manual(values = c("Library recovered" = 16,
                                "Library not recovered" = 2)) +
  geom_text_repel(aes(label = Site), size = 6, box.padding = 0.6,
                  min.segment.length = 0, max.overlaps = Inf) +
  geom_text(data = supp1_stats, aes(x_pos, y_pos, label = label),
            parse = TRUE, hjust = 0, vjust = 1.3, size = 6, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  facet_wrap(~ Panel, scales = "free", ncol = 2) +
  labs(x = expression(log[10]*"(DNA yield, ng) or "*log[10]*"(sequencing reads)"),
       y = expression(log[10]*"(cells m"^{-3}*")"),
       shape = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 7))) +
  theme_bw(base_size = 19) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 19),
        panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", linewidth = 0.8),
        panel.spacing    = unit(1.2, "lines"),
        axis.title       = element_text(size = 20),
        axis.text        = element_text(size = 17, colour = "black"),
        axis.ticks       = element_line(colour = "black"),
        legend.position  = "bottom",
        legend.text      = element_text(size = 18),
        legend.key.size  = unit(1.4, "lines"))

ggsave(file.path(OUT_DIR, "SuppFig1_cells_DNA_reads.pdf"),
       plot = supp1, width = 11, height = 10)


library(tidyverse)
library(ggrepel)

df <- tibble(
  Site  = c("BB", "MK", "VL", "RM", "DF"),
  TCC   = c(56003.39, 6134.244, 247268.2, 2725.61, 1507207),
  reads = c(6430565, 2630772, 10094145, 1866743, 98902286)
)

ct <- cor.test(df$TCC, df$reads, method = "spearman", exact = TRUE)

lab <- sprintf(
  "rho == %.2f*','~~italic(p) == %.3f*','~~n == %d",
  ct$estimate, ct$p.value, nrow(df)
)

print(ct)

p <- ggplot(df, aes(x = log10(reads), y = log10(TCC))) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1.1,
    colour = scales::alpha("blue", 0.5)
  ) +
  geom_point(
    size = 5,
    alpha = 0.7
  ) +
  geom_text_repel(
    aes(label = Site),
    size = 6,
    box.padding = 0.6,
    point.padding = 0.7,
    force = 1.5,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    vjust = 1.4,
    parse = TRUE,
    size = 5.5,
    label = lab
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.25))
  ) +
  labs(
    x = expression(log[10]*"(sequencing reads)"),
    y = expression(log[10]*"(TCC, cells m"^{-3}*")")
  ) +
  theme_bw(base_size = 17) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.7),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 16, colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

ggsave(
  file.path(OUT_DIR, "Fig2C_TCC_vs_reads.pdf"),
  plot = p,
  width = 5,
  height = 7
)




# =============================================================================
# 2. FIGURE 2: flow cytometry
# =============================================================================

add_fcm <- function(site, sample, type, vol, tcc, icc)
  tibble(Site = site, Sample = sample, Type = type, cone_vol_mL = vol,
         Rep = 1:3, TCC_150uL = tcc, ICC_150uL = icc)

fcm_raw <- bind_rows(
  # Vancouver Landfill
  add_fcm("VL","S1","air",13.14,c(729,746,775),c(446,462,453)),
  add_fcm("VL","S3","air",13.51,c(585,591,520),c(173,176,164)),
  add_fcm("VL","S6","air",13.27,c(844,845,854),c(177,173,174)),
  add_fcm("VL","S8","air",13.43,c(47417,47158,48139),c(27128,27484,27281)),
  add_fcm("VL","NC1","control",15.07,c(34,2,1),c(0,0,0)),
  add_fcm("VL","NC2","control",15.16,c(1,0,1),c(0,3,2)),
  
  # Dairy Farm
  add_fcm("DF","S2","air",11.66,c(74501,73309,74362),c(3302,3135,3228)),
  add_fcm("DF","S6","air",11.65,c(76377,77309,77397),c(3364,3456,3306)),
  add_fcm("DF","S7","air",11.96,c(102094,100278,100442),c(9877,5295,5185)),
  add_fcm("DF","S9","air",11.71,c(92939,96006,93214),c(6114,5933,5966)),
  add_fcm("DF","NC1","control",15.80,c(17,12,17),c(0,0,1)),
  add_fcm("DF","NC2","control",15.32,c(9,9,14),c(1,0,0)),
  
  # Minnekhada
  add_fcm("MK","S2","air",13.22,c(359,329,342),c(311,349,330)),
  add_fcm("MK","S4","air",13.45,c(322,316,300),c(321,282,261)),
  add_fcm("MK","S6","air",13.18,c(308,318,312),c(108,153,109)),
  add_fcm("MK","S9","air",13.08,c(296,278,274),c(99,87,85)),
  add_fcm("MK","NC1","control",15.16,c(1,0,0),c(0,0,0)),
  add_fcm("MK","NC2","control",15.27,c(0,0,1),c(1,0,0)),
  
  # Burns Bog
  add_fcm("BB","S1","air",10.89,c(1305,1290,1320),c(6.667,3.333,3.333)),
  add_fcm("BB","S3","air",10.13,c(11477,11321,11058),c(159,146,169)),
  add_fcm("BB","S6","air",10.37,c(1495,1408,1509),c(138,166,152)),
  add_fcm("BB","S8","air",10.45,c(713,716,703),c(64,61,59)),
  add_fcm("BB","NC1","control",15.36,c(0,0,0),c(0,0,0)),
  add_fcm("BB","TC","control",15.42583333,c(1,0,0),c(0,0,1)),
  
  # Richmond Marsh
  add_fcm("RM","S3","air",11.50,c(161,166,193),c(4,3,4)),
  add_fcm("RM","S4","air",10.97,c(125,117,94),c(2,3,3)),
  add_fcm("RM","S8","air",11.95,c(165,142,174),c(2,0,2)),
  add_fcm("RM","S9","air",11.67,c(208,193,170),c(1,3,0)),
  add_fcm("RM","NC1","control",15.01,c(0,0,1),c(1,2,0)),
  add_fcm("RM","NC2","control",15.32,c(1,0,1),c(0,0,0)),
  
  # SFU
  add_fcm("SFU","S2","air",12.03,c(82,91,98),c(13,14,16)),
  add_fcm("SFU","S4","air",12.41,c(74,70,96),c(15,13,17)),
  add_fcm("SFU","S5","air",12.44,c(126,116,95),c(30,23,25)),
  add_fcm("SFU","S9","air",11.17,c(89,98,108),c(10,12,14)),
  add_fcm("SFU","NC1","control",15.21,c(3,2,2),c(2,1,2)),
  add_fcm("SFU","NC2","control",15.23,c(5,5,4),c(3,2,2)),
  
  # CH
  add_fcm("CH","S1","air",12.30583333,c(301,317,332),c(74,72,73)),
  add_fcm("CH","S2","air",11.50583333,c(308,315,291),c(55,56,55)),
  add_fcm("CH","S3","air",11.94583333,c(256,263,272),c(60,59,54)),
  add_fcm("CH","S4","air",11.45583333,c(298,297,289),c(61,63,48)),
  add_fcm("CH","NC1","control",15.42583333,c(1,4,1),c(1,2,1)),
  add_fcm("CH","TC","control",15.42583333,c(0,1,0),c(1,0,0))
)

# Left-to-right site order
fcm_sites <- c("VL","DF","MK","BB","RM","SFU","CH")
stopifnot(setequal(unique(fcm_raw$Site), fcm_sites))

to_m3 <- function(count150, vol_mL) count150 * (1000 / 150) * vol_mL / 4.5

fcm_raw <- fcm_raw %>%
  mutate(Site = factor(Site, levels = fcm_sites),
         TCC = to_m3(TCC_150uL, cone_vol_mL),
         ICC = to_m3(ICC_150uL, cone_vol_mL))

# Negative-control maximum per site
ctrl_max <- fcm_raw %>%
  filter(Type == "control") %>%
  group_by(Site) %>%
  summarise(ctrl_max = max(c(TCC, ICC)), .groups = "drop")

# Collection means
coll <- fcm_raw %>%
  filter(Type == "air") %>%
  group_by(Site, Sample) %>%
  summarise(TCC_mean = mean(TCC), ICC_mean = mean(ICC), .groups = "drop") %>%
  left_join(ctrl_max, by = "Site") %>%
  mutate(ICC_pct = 100 * ICC_mean / TCC_mean,
         near_det = ICC_mean <= ctrl_max) %>%
  arrange(Site, Sample) %>%
  group_by(Site) %>%
  mutate(coll_idx = factor(row_number(), levels = 1:4),
         x = as.numeric(Site) + (row_number() - 2.5) * 0.18) %>%
  ungroup()

# Checks
print(coll %>%
        select(Site, Sample, TCC_mean, ICC_mean, ICC_pct, ctrl_max, near_det),
      n = Inf)

print(coll %>%
        group_by(Site) %>%
        summarise(TCC = mean(TCC_mean),
                  ICC = mean(ICC_mean),
                  ICC_pct = mean(ICC_pct),
                  .groups = "drop"))

fcm_raw %>%
  filter(Type == "air") %>%
  pivot_longer(c(TCC_150uL, ICC_150uL),
               names_to = "Metric", values_to = "events") %>%
  group_by(Site, Sample, Metric) %>%
  summarise(mean_events = mean(events),
            cv = 100 * sd(events) / mean(events),
            .groups = "drop") %>%
  mutate(event_class = if_else(mean_events >= 10, ">= 10 events", "< 10 events")) %>%
  group_by(event_class) %>%
  summarise(n = n(), median_cv = median(cv), max_cv = max(cv), .groups = "drop") %>%
  print()

pa_fcm <- coll %>%
  pivot_longer(c(TCC_mean, ICC_mean),
               names_to = "Metric", values_to = "cells") %>%
  mutate(Metric = factor(sub("_mean", "", Metric), levels = c("TCC","ICC")))

site_sum <- coll %>%
  group_by(Site) %>%
  summarise(m = mean(ICC_pct), s = sd(ICC_pct), .groups = "drop") %>%
  mutate(xs = as.numeric(Site))

shades <- c("1"="#6baed6", "2"="#3182bd", "3"="#08519c", "4"="#08306b")

x_scale <- scale_x_continuous(
  breaks = 1:7, labels = fcm_sites,
  limits = c(0.5, 7.5), expand = c(0,0)
)

# Separates environmental sites from SFU/CH
failed_line <- geom_vline(
  xintercept = 5.5, linetype = "dashed",
  colour = "grey55", linewidth = 0.4
)

fcm_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))

# Panel A
fig2A <- ggplot() +
  geom_segment(data = coll,
               aes(x, xend = x, y = ICC_mean, yend = TCC_mean, colour = coll_idx),
               linewidth = 0.5) +
  geom_point(data = pa_fcm,
             aes(x, cells, colour = coll_idx, shape = Metric),
             fill = "white", size = 5, alpha = 0.7, stroke = 1) +
  failed_line +
  annotate("text", x = 6.5, y = 10^6.7,
           label = "Aerosol library\nnot recovered",
           size = 4, colour = "grey35", lineheight = 0.9) +
  x_scale +
  scale_y_log10(limits = c(10,1e7),
                breaks = 10^(1:7),
                labels = scales::label_log()) +
  scale_shape_manual(values = c(TCC=16, ICC=21),
                     breaks = c("TCC","ICC"),
                     labels = c("TCC (total cells)", "ICC (intact cells)"),
                     name = NULL) +
  scale_colour_manual(values = shades, guide = "none") +
  guides(shape = guide_legend(
    override.aes = list(colour = "grey20", size = 3.5))) +
  labs(x = NULL, y = expression("Cells m"^{-3}*" air")) +
  fcm_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Panel B
fig2B <- ggplot() +
  geom_errorbar(data = site_sum,
                aes(xs, ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.25, linewidth = 0.5, colour = "grey20") +
  geom_segment(data = site_sum,
               aes(x = xs - 0.25, xend = xs + 0.25, y = m, yend = m),
               linewidth = 0.9, colour = "grey20") +
  geom_point(data = coll,
             aes(x, ICC_pct, colour = coll_idx),
             size = 5, alpha = 0.9) +
  geom_text(data = site_sum,
            aes(xs, pmin(m + s, 100) + 9,
                label = sprintf("%.0f%%", m)),
            size = 5, colour = "grey20") +
  failed_line +
  x_scale +
  scale_y_continuous(limits = c(0,110),
                     breaks = seq(0,100,25),
                     labels = function(x) paste0(x,"%")) +
  scale_colour_manual(values = shades, guide = "none") +
  labs(x = NULL, y = "ICC / TCC") +
  fcm_theme +
  theme(axis.text.x = element_text(face = "bold", size = 14))

# Combine
fig2 <- (fig2A / fig2B) +
  plot_layout(heights = c(1.6,1), guides = "collect") &
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.background = element_rect(
          colour = "black", fill = "white", linewidth = 0.5),
        legend.margin = margin(4,8,4,8))

save_pdf(fig2, "Fig2_flow_cytometry.pdf", 5.5, 7)

# =============================================================================
# 3. FIGURE 3: community diversity and composition (genus level)
#    - A–D: rarefied tables (Shannon Welch t-test; NMDS + PERMANOVA; betadisper)
#    - E:   non-rarefied relative abundance (Procrustes/PROTEST; residual Wilcoxon)
#    - One legend for all panels: colour = site; shape = sample type;
#      filled = bacteria, open = eukaryotes
# =============================================================================

key_levels <- c("Bacteria, Aerosol", "Bacteria, Paired environment",
                "Eukaryota, Aerosol", "Eukaryota, Paired environment")
key_shapes <- setNames(c(17, 16, 2, 1), key_levels)
pt_alpha <- 0.75; pt_stroke <- 1.2

fig3_theme <- theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", linewidth = 0.6),
        axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

common_scales <- list(   # identical scales -> patchwork merges into ONE legend
  scale_colour_manual(values = site_cols, limits = site_levels, name = "Site"),
  scale_shape_manual(values = key_shapes, limits = key_levels, name = NULL),
  guides(colour = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 16, size = 4, alpha = 1)),
         shape  = guide_legend(nrow = 2, byrow = TRUE, order = 2,
                               override.aes = list(size = 3.5, alpha = 1, stroke = pt_stroke))))

prep <- function(rds, domain) {   # load, collapse to genus, add Site/Type/Key
  ps   <- tax_glom(readRDS(file.path(PS_DIR, rds)), taxrank = "genus", NArm = TRUE)
  meta <- data.frame(sample_data(ps))
  meta$Sample <- sample_names(ps)
  meta$Site   <- factor(sub("^[AWS]_([A-Z]+).*", "\\1", meta$Sample), levels = site_levels)
  meta$Type   <- factor(if_else(meta$Sample_type2 == "Air", "Aerosol", "Paired environment"),
                        levels = c("Aerosol", "Paired environment"))
  meta$Key    <- factor(paste(domain, meta$Type, sep = ", "), levels = key_levels)
  list(ps = ps, meta = meta)
}

make_shannon <- function(d) {
  a <- estimate_richness(d$ps, measures = "Shannon") %>% rownames_to_column("Sample") %>%
    left_join(select(d$meta, Sample, Site, Type, Key), by = "Sample")
  tt <- t.test(Shannon ~ Type, data = a)                                 # Welch t-test
  ggplot(a, aes(Type, Shannon)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, colour = "grey30") +
    geom_point(aes(colour = Site, shape = Key), size = 3.5, alpha = pt_alpha, stroke = pt_stroke,
               position = position_jitter(width = 0.12, seed = 1)) +
    annotate("text", x = 1.5, y = max(a$Shannon) * 1.06, parse = TRUE, size = 4.5,
             label = paste0("italic(p)~'= ", fmt_p(tt$p.value), "'")) +
    common_scales + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "Shannon diversity") + fig3_theme
}

make_nmds <- function(d, domain) {
  bray <- phyloseq::distance(d$ps, method = "bray")
  set.seed(1); perm <- vegan::adonis2(bray ~ Type, data = d$meta, permutations = 9999)
  set.seed(1); bd_p <- vegan::permutest(vegan::betadisper(bray, d$meta$Type), permutations = 9999)
  set.seed(6); ord  <- vegan::metaMDS(bray, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE)
  cat(sprintf("%s | PERMANOVA F(%d,%d)=%.2f R2=%.3f p=%s | betadisper F=%.2f p=%.3f | stress=%.3f\n",
              domain, perm$Df[1], perm$Df[2], perm$F[1], perm$R2[1], fmt_p(perm$`Pr(>F)`[1]),
              bd_p$tab$F[1], bd_p$tab$`Pr(>F)`[1], ord$stress))
  pts <- as.data.frame(ord$points) %>% rownames_to_column("Sample") %>%
    left_join(select(d$meta, Sample, Site, Type, Key), by = "Sample")
  ggplot(pts, aes(MDS1, MDS2, colour = Site, shape = Key)) +
    geom_point(size = 4.5, alpha = pt_alpha, stroke = pt_stroke) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, parse = TRUE, size = 4.2,
             label = paste0("italic(R)^2~'= ", sprintf("%.3f", perm$R2[1]),
                            ",'~italic(p)~'= ", fmt_p(perm$`Pr(>F)`[1]), "'")) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8, size = 3.8,
             label = paste0("Stress = ", round(ord$stress, 2))) +
    common_scales + scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(x = "NMDS1", y = "NMDS2") + fig3_theme
}

# Panels A–D (rarefied)
bac <- prep("psf_bacteria_raref_min.rds",  "Bacteria")
euk <- prep("psf_eukaryota_raref_min.rds", "Eukaryota")
fig3A <- make_shannon(bac); fig3B <- make_nmds(bac, "Bacteria")
fig3C <- make_shannon(euk); fig3D <- make_nmds(euk, "Eukaryota")

# Panel E (non-rarefied relative abundance)
bac_nr <- prep("psf_bacteria_nonrarefied.rds",  "Bacteria")
euk_nr <- prep("psf_eukaryota_nonrarefied.rds", "Eukaryota")
shared <- intersect(sample_names(bac_nr$ps), sample_names(euk_nr$ps))
rel_bray <- function(ps) phyloseq::distance(transform_sample_counts(prune_samples(shared, ps),
                                                                    function(x) x / sum(x)), method = "bray")
set.seed(1)
bac_nmds <- vegan::metaMDS(rel_bray(bac_nr$ps), k = 2, trymax = 100, trace = FALSE)
euk_nmds <- vegan::metaMDS(rel_bray(euk_nr$ps), k = 2, trymax = 100, trace = FALSE)
set.seed(1); prot <- vegan::protest(bac_nmds, euk_nmds, permutations = 999)
proc <- vegan::procrustes(bac_nmds, euk_nmds)

proc_df <- tibble(Sample = rownames(proc$X), bac_x = proc$X[, 1], bac_y = proc$X[, 2],
                  euk_x = proc$Yrot[, 1], euk_y = proc$Yrot[, 2],
                  residual = as.numeric(residuals(proc))) %>%
  left_join(select(bac_nr$meta, Sample, Site, Type), by = "Sample")
cat(sprintf("PROTEST r=%.3f p=%s\n", prot$t0, fmt_p(prot$signif)))
print(wilcox.test(residual ~ Type, data = proc_df, exact = TRUE))       # residuals: aerosol vs paired env

pts_E <- bind_rows(
  transmute(proc_df, Site, Type, x = bac_x, y = bac_y, Domain = "Bacteria"),
  transmute(proc_df, Site, Type, x = euk_x, y = euk_y, Domain = "Eukaryota")) %>%
  mutate(Key = factor(paste(Domain, Type, sep = ", "), levels = key_levels))

fig3E <- ggplot() +
  geom_segment(data = proc_df, aes(euk_x, euk_y, xend = bac_x, yend = bac_y, colour = Site),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               linewidth = 0.6, alpha = pt_alpha, show.legend = FALSE) +
  geom_point(data = pts_E, aes(x, y, colour = Site, shape = Key), size = 3.5, alpha = pt_alpha, stroke = pt_stroke) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, parse = TRUE, size = 4.2,
           label = paste0("italic(r)~'= ", sprintf("%.3f", prot$t0), ",'~italic(p)~'= ", fmt_p(prot$signif), "'")) +
  common_scales + scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "Dimension 1", y = "Dimension 2") + fig3_theme

# Save: panels without legends (for PowerPoint), A–D preview, combined with one legend
no_leg <- theme(legend.position = "none")
save_pdf(fig3A + no_leg, "Fig3A_bacteria_shannon.pdf",  3.6, 4.2)
save_pdf(fig3B + no_leg, "Fig3B_bacteria_NMDS.pdf",     5.0, 4.2)
save_pdf(fig3C + no_leg, "Fig3C_eukaryota_shannon.pdf", 3.6, 4.2)
save_pdf(fig3D + no_leg, "Fig3D_eukaryota_NMDS.pdf",    5.0, 4.2)
save_pdf(fig3E + no_leg, "Fig3E_procrustes.pdf",        5.6, 8)
save_pdf(((fig3A + fig3B + plot_layout(widths = c(1, 1.4))) /
            (fig3C + fig3D + plot_layout(widths = c(1, 1.4)))) & no_leg, "Fig3_AD_preview.pdf", 7.5, 8.5)
save_pdf(fig3A + fig3B + fig3C + fig3D + fig3E +
           plot_layout(design = "ABE\nCDE", widths = c(1, 1.4, 1.6), guides = "collect") &
           theme(legend.position = "bottom", legend.box = "vertical"), "Fig3_combined.pdf", 14, 9.5)


# =============================================================================
# 4. FIGURE 4 NUMBERS: site-blocked DESeq2 (~ Site + Compartment)
#    - Genus level, non-rarefied; genera with < 10 total counts removed
#    - poscounts size factors; apeglm LFC shrinkage
#    - Significant: padj < 0.05 and |log2FC| >= 2
# =============================================================================

get_deseq <- function(rds) {
  ps   <- readRDS(file.path(PS_DIR, rds))
  meta <- as(sample_data(ps), "data.frame")
  sample_data(ps)$Compartment <- factor(ifelse(meta$Sample_type2 == "Air", "Air", "Reference"),
                                        levels = c("Reference", "Air"))
  sample_data(ps)$Site <- factor(sub("^[AWS]_([A-Z]+).*", "\\1", rownames(meta)))
  ps_gen <- tax_glom(ps, taxrank = "genus", NArm = TRUE)
  ps_gen <- prune_taxa(taxa_sums(ps_gen) > 0, ps_gen)
  dds <- phyloseq_to_deseq2(ps_gen, ~ Site + Compartment)
  dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
  coef <- grep("Compartment_Air", DESeq2::resultsNames(dds), value = TRUE)
  res  <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm")
  tx   <- as.data.frame(tax_table(ps_gen)) %>% rownames_to_column("taxon")
  as.data.frame(res) %>% rownames_to_column("taxon") %>%
    left_join(tx[, c("taxon", "genus")], by = "taxon") %>%
    mutate(padj = replace_na(padj, 1), log2FoldChange = replace_na(log2FoldChange, 0))
}

count_sig <- function(df, label) {
  sig <- df$padj < 0.05 & abs(df$log2FoldChange) >= 2
  tibble(domain = label, tested = nrow(df), significant = sum(sig),
         aerosol_enriched = sum(sig & df$log2FoldChange > 0),
         environment_enriched = sum(sig & df$log2FoldChange < 0),
         pct_not_significant = round(100 * mean(!sig), 1))
}

bac_de <- get_deseq("psf_bacteria_nonrarefied.rds")
euk_de <- get_deseq("psf_eukaryota_nonrarefied.rds")
print(bind_rows(count_sig(bac_de, "Bacteria"), count_sig(euk_de, "Eukaryota")))


# =============================================================================
# 5. SUPPLEMENTARY TABLE S2: domain read proportions (% of classified reads)
#    - Bracken table; human reads and "AC" samples removed (as in main pipeline)
# =============================================================================

bracken  <- read.delim(file.path(R_WORK, "bracken_with_taxonomy.tsv"), sep = "\t", quote = "",
                       check.names = FALSE, stringsAsFactors = FALSE)
tax_cols <- c("taxa","name","taxid","superkingdom","phylum","class","order","family","genus","species")

dom <- bracken %>%
  filter(!grepl("\\bHomo\\b", genus)) %>%
  select(superkingdom, all_of(setdiff(colnames(bracken), tax_cols))) %>%
  pivot_longer(-superkingdom, names_to = "Sample", values_to = "reads") %>%
  mutate(Sample = sub("\\.bracken_num$", "", Sample)) %>%
  filter(!grepl("^AC", Sample), superkingdom %in% c("Bacteria", "Archaea", "Eukaryota", "Viruses")) %>%
  group_by(Sample, superkingdom) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  group_by(Sample) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup() %>%
  mutate(Type = if_else(grepl("^A_", Sample), "Aerosol", "Paired environment"))

dom_wide <- dom %>% select(Sample, Type, superkingdom, pct) %>%
  pivot_wider(names_from = superkingdom, values_from = pct) %>% arrange(Type, Sample)
print(dom_wide, n = Inf)
write_csv(dom_wide, file.path(OUT_DIR, "TableS2_domain_read_proportions.csv"))


# =============================================================================
# 6. FIGURE 5A: C1 metabolism marker gene heatmap
#    - KO x sample matrix of proportional TPM (Salmon -> eggNOG KO, prokaryotic)
#    - Colour = within-site z-score per gene; grey = not detected
#    - Genes never detected in any sample are dropped from the plot
#    - Note: K10944–K10946 are shared by pmoCAB and amoCAB in KEGG
# =============================================================================

ko_mat  <- read.delim(file.path(KO_DIR, "KOxSample_TPM.prop.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
ko_mat[] <- lapply(ko_mat, function(x) as.numeric(as.character(x)))
ko_func <- read.delim(file.path(KO_DIR, "KO_to_Function.clean3.tsv"), header = FALSE, sep = "\t",
                      quote = "\"", col.names = c("KO", "Preferred", "Description"))
ko_mat  <- ko_mat[rownames(ko_mat) %in% ko_func$KO, ]

c1_map <- tribble(
  ~Category, ~Gene, ~KO,
  "Methane oxidation","pmoC","K10946", "Methane oxidation","pmoA","K10944", "Methane oxidation","pmoB","K10945",
  "Methane oxidation","mmoX","K16157", "Methane oxidation","mmoY","K16158", "Methane oxidation","mmoZ","K16159",
  "Methane oxidation","mmoD","K16162", "Methane oxidation","mmoC","K16161",
  "Methanol oxidation","xoxF","K10823", "Methanol oxidation","mxaF","K14028", "Methanol oxidation","mxaI","K14029",
  "Formaldehyde assimilation & detoxification","fae","K10798",
  "Formaldehyde assimilation & detoxification","fch","K01455",
  "Formaldehyde assimilation & detoxification","mch","K08691",
  "Formaldehyde assimilation & detoxification","folD","K01491",
  "Formaldehyde assimilation & detoxification","frmA","K00122",
  "Formaldehyde assimilation & detoxification","frmB","K01468",
  "Formate oxidation","fdhA","K00123", "Formate oxidation","fdhB","K00124", "Formate oxidation","fdhC","K00125",
  "C1 assimilation (RuMP pathway)","hps","K13812", "C1 assimilation (RuMP pathway)","phi","K02331",
  "C1 assimilation (Serine cycle)","glyA","K00600", "C1 assimilation (Serine cycle)","mcl","K01638",
  "C1 assimilation (Serine cycle)","sga","K01783",  "C1 assimilation (Serine cycle)","mtkA","K00614",
  "C1 assimilation (Serine cycle)","mtkB","K00615",
  "CO oxidation","coxL","K03518", "CO oxidation","coxM","K03519", "CO oxidation","coxS","K03520",
  "Calvin cycle","rbcL","K01601", "Calvin cycle","rbcS","K01602", "Calvin cycle","prkB","K00855",
  "Hydrogen metabolism","hoxK","K06282", "Hydrogen metabolism","hoxG","K06281", "Hydrogen metabolism","hoxZ","K06283",
  "Hydrogen metabolism","hupS","K18016", "Hydrogen metabolism","hupL","K18017", "Hydrogen metabolism","hupC","K18018",
  "Hydrogen metabolism","hoxF","K06284", "Hydrogen metabolism","hoxU","K06285", "Hydrogen metabolism","hoxY","K06286",
  "Hydrogen metabolism","hoxH","K06287", "Hydrogen metabolism","hoxB","K06288", "Hydrogen metabolism","hoxC","K06289",
  "Hydrogen metabolism","hoxJ","K06290", "Hydrogen metabolism","hoxA","K06291", "Hydrogen metabolism","hhyS","K22622",
  "Hydrogen metabolism","hhyL","K22623", "Hydrogen metabolism","hypA","K04651", "Hydrogen metabolism","hypB","K04652",
  "Hydrogen metabolism","hypC","K04653", "Hydrogen metabolism","hypD","K04654", "Hydrogen metabolism","hypE","K04655",
  "Hydrogen metabolism","hypF","K04656", "Hydrogen metabolism","hupD","K06280", "Hydrogen metabolism","hydA","K00533")

category_levels <- unique(c1_map$Category)                  # order as listed above
gene_levels     <- c(unique(c1_map$Gene)[1:14], "mtdA", unique(c1_map$Gene)[15:57])  # mtdA: no KO (placeholder)

infer_site <- function(s)                                    # base R grepl (avoids stringi issues)
  case_when(grepl("_VL", s) ~ "VL", grepl("_DF", s) ~ "DF", grepl("_MK", s) ~ "MK",
            grepl("_BB", s) ~ "BB", grepl("_RM", s) ~ "RM", TRUE ~ NA_character_)

ko_long <- ko_mat %>% rownames_to_column("KO") %>%
  pivot_longer(-KO, names_to = "Sample", values_to = "prop") %>%
  mutate(Site = factor(infer_site(Sample), levels = site_levels)) %>% filter(!is.na(Site))

# x-axis labels: Aerosol, then Soil 1/2 (or Water 1/2 at BB)
sample_tbl <- ko_long %>% distinct(Site, Sample) %>%
  mutate(is_air = grepl("^A_", Sample),
         env_type = case_when(grepl("^W_", Sample) ~ "Water", grepl("^S_", Sample) ~ "Soil")) %>%
  arrange(Site, desc(is_air), Sample) %>% group_by(Site) %>%
  mutate(Label = if_else(is_air, "Aerosol", paste(env_type, cumsum(!is_air)))) %>% ungroup() %>%
  mutate(Sample_in_site = factor(Sample, levels = unique(Sample)))
sample_labels <- setNames(sample_tbl$Label, sample_tbl$Sample)

gene_long <- c1_map %>% left_join(ko_long, by = "KO", relationship = "one-to-many") %>%
  group_by(Category, Gene, Site, Sample) %>% summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop")

gene_long <- bind_rows(distinct(c1_map, Category, Gene),
                       tibble(Category = "Formaldehyde assimilation & detoxification", Gene = "mtdA")) %>%
  tidyr::crossing(select(sample_tbl, Site, Sample)) %>%
  left_join(gene_long, by = c("Category", "Gene", "Site", "Sample")) %>%
  mutate(prop = replace_na(prop, 0)) %>%
  left_join(select(sample_tbl, Sample, Sample_in_site), by = "Sample") %>%
  mutate(Category = factor(Category, levels = category_levels),
         GeneLabel = factor(Gene, levels = rev(gene_levels))) %>%
  group_by(Gene) %>% filter(any(prop > 0)) %>% ungroup() %>%           # drop never-detected genes
  group_by(Site, Gene) %>%
  mutate(z = if (sd(prop) == 0) 0 else (prop - mean(prop)) / sd(prop),
         z_plot = ifelse(prop == 0, NA, z)) %>% ungroup()

vmax <- max(abs(gene_long$z_plot), na.rm = TRUE)

# ---- x positions: samples grouped by site, small gap between sites ----------
gap <- 0.6
x_tbl <- sample_tbl %>%
  arrange(Sample_in_site) %>%
  group_by(Site) %>%
  mutate(idx = row_number()) %>%
  ungroup() %>%
  mutate(x = (as.integer(Site) - 1) * (3 + gap) + idx)

site_centers <- x_tbl %>% group_by(Site) %>% summarise(xc = mean(x), .groups = "drop")
x_min <- min(x_tbl$x) - 0.5
x_max <- max(x_tbl$x) + 0.5

# ---- y positions: header row per category, then its genes (top to bottom) ---
genes_kept <- gene_long %>%
  distinct(Category, Gene) %>%
  mutate(Category = factor(Category, levels = category_levels),
         Gene     = factor(Gene, levels = gene_levels)) %>%
  arrange(Category, Gene) %>%
  mutate(across(c(Category, Gene), as.character))

rows <- bind_rows(lapply(unique(genes_kept$Category), function(cat)
  tibble(Category = cat,
         Row  = c(paste0("HDR_", cat), genes_kept$Gene[genes_kept$Category == cat]),
         Type = c("header", rep("gene", sum(genes_kept$Category == cat))))))

hdr_h <- 1.8                                 # header row height (gene row = 1)
rows <- rows %>%
  mutate(h   = if_else(Type == "header", hdr_h, 1),
         top = sum(h) - cumsum(h) + h,       # top edge of each row (first row at top)
         y   = top - h / 2,                  # row centre
         ymin = top - h, ymax = top)
hdr  <- filter(rows, Type == "header")
gene_y <- filter(rows, Type == "gene") %>% select(Gene = Row, y)

plot_df <- gene_long %>%
  mutate(Gene = as.character(Gene)) %>%
  left_join(select(x_tbl, Sample, x), by = "Sample") %>%
  left_join(gene_y, by = "Gene")

# ---- Plot ------------------------------------------------------------------
fig5A <- ggplot() +
  # header bars across the heatmap
  geom_rect(data = hdr, aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax),
            fill = "grey92", colour = NA) +
  geom_text(data = hdr, aes(x = x_min + 0.15, y = y, label = Category),
            hjust = 0, size = 12 / .pt, fontface = "bold") +
  # heatmap tiles
  geom_tile(data = plot_df, aes(x, y, fill = z_plot), colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-vmax, vmax), name = "z-score (within site)", na.value = "grey85",
                       guide = guide_colorbar(barwidth = 10, barheight = 0.8,
                                              title.position = "top", title.hjust = 0.5)) +
  scale_x_continuous(breaks = x_tbl$x, labels = x_tbl$Label, expand = c(0, 0),
                     sec.axis = dup_axis(breaks = site_centers$xc,
                                         labels = as.character(site_centers$Site))) +
  scale_y_continuous(breaks = gene_y$y,
                     labels = parse(text = paste0("italic('", gene_y$Gene, "')")),
                     expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, colour = "black"),
        axis.text.x.top    = element_text(size = 16, face = "bold", colour = "black"),
        axis.ticks.x.top   = element_blank(),
        axis.text.y  = element_text(size = 13, colour = "black"),
        axis.ticks   = element_line(colour = "grey70", linewidth = 0.3),
        panel.grid   = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.4),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 11),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(OUT_DIR, "Fig5A_C1_metabolism_heatmap.pdf"), plot = fig5A,
       width = 200, height = 260, units = "mm")

# =============================================================================
# 7. SESSION INFO (record package versions; report in Methods if needed)
# =============================================================================



# =============================================================================
# SUPPLEMENTARY FIG. 2: domain composition (% of classified reads) per sample
#   - Uses `dom` from Section 5 (Bracken; human reads removed)
#   - Bars: one per sample; facets: sites (high -> low CH4)
# =============================================================================

dom_plot <- dom %>%
  mutate(Site   = factor(sub("^[AWS]_([A-Z]+).*", "\\1", Sample), levels = site_levels),
         is_air = grepl("^A_", Sample),
         env    = case_when(grepl("^W_", Sample) ~ "Water", grepl("^S_", Sample) ~ "Soil")) %>%
  group_by(Site) %>%
  arrange(desc(is_air), Sample, .by_group = TRUE) %>%
  mutate(Label = if_else(is_air, "Aerosol",
                         paste(env, match(Sample, unique(Sample[!is_air]))))) %>%
  ungroup() %>%
  mutate(Label = factor(Label, levels = c("Aerosol", "Soil 1", "Soil 2", "Water 1", "Water 2")),
         superkingdom = factor(superkingdom,
                               levels = c("Viruses", "Archaea", "Eukaryota", "Bacteria")))

# Summary for the Results text (send this output)
dom_plot %>%
  group_by(Type, superkingdom) %>%
  summarise(mean = mean(pct), min = min(pct), max = max(pct), .groups = "drop") %>%
  print(n = Inf)

supp2 <- ggplot(dom_plot, aes(Label, pct, fill = superkingdom)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
  facet_grid(~ Site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(Bacteria = "#4E79A7", Eukaryota = "#F28E2B",
                               Archaea = "#59A14F", Viruses = "#B07AA1"),
                    breaks = c("Bacteria", "Eukaryota", "Archaea", "Viruses"),
                    name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Classified reads (%)") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 16))

ggsave(file.path(OUT_DIR, "SuppFig2_domain_composition.pdf"),
       plot = supp2, width = 11, height = 6)




# =============================================================================
# SUPPLEMENTARY TABLE S3: all differentially abundant genera (Fig. 4)
#   - From site-blocked DESeq2 results (Section 4)
#   - Significant: padj < 0.05 and |shrunken log2FC| >= 2
# =============================================================================

tableS3 <- bind_rows(
  bac_de %>% mutate(Domain = "Bacteria"),
  euk_de %>% mutate(Domain = "Eukaryota")
) %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 2) %>%
  mutate(Direction = if_else(log2FoldChange > 0, "Enriched in aerosols", "Depleted in aerosols")) %>%
  transmute(Domain, Genus = genus, Direction,
            log2FC = round(log2FoldChange, 2),
            padj   = signif(padj, 3),
            baseMean = round(baseMean, 1)) %>%
  arrange(Domain, Direction, padj)

# Counts should match the text: Bacteria 141 enriched / 7 depleted; Eukaryota 5 / 2
print(count(tableS3, Domain, Direction))

write_csv(tableS3, file.path(OUT_DIR, "TableS3_differentially_abundant_genera.csv"))


# =============================================================================
# 6b. FIGURE 5B: pmoCAB log2(aerosol / paired environment) vs ambient CH4
#   - Sum of pmoA, pmoB, pmoC proportional TPM per sample
#   - Paired environment = mean of the two soil/water replicates per site
#   - Sites without pmoCAB in the aerosol are shown as "n.d." (no ratio computed)
#   - Uses ko_long (Section 6) and site_cols (Section 0)
#   - Note: these KOs are shared by pmoCAB and amoCAB in KEGG
# =============================================================================

pmo_kos <- c(pmoA = "K10944", pmoB = "K10945", pmoC = "K10946")
ch4     <- tibble(Site    = c("VL", "DF", "MK", "BB", "RM"),
                  CH4_ppm = c(85.98, 26.76, 2.62, 2.60, 2.29))

pmo_df <- ko_long %>%
  filter(KO %in% pmo_kos) %>%
  group_by(Site, Sample) %>%
  summarise(pmo = sum(prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = if_else(grepl("^A_", Sample), "Aerosol", "Environment")) %>%
  group_by(Site, Type) %>%
  summarise(pmo = mean(pmo), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = pmo) %>%
  mutate(Site = as.character(Site)) %>%
  left_join(ch4, by = "Site") %>%
  mutate(detected = Aerosol > 0,
         log2FC   = if_else(detected, log2(Aerosol / Environment), NA_real_))

print(pmo_df)   # check: VL and DF detected; BB, MK, RM not detected

# "n.d." row placed below the detected values
nd_y     <- floor(min(pmo_df$log2FC, na.rm = TRUE)) - 1.5
y_breaks <- pretty(c(min(pmo_df$log2FC, na.rm = TRUE), 0))
y_breaks <- y_breaks[y_breaks > nd_y + 0.75]
pmo_df   <- mutate(pmo_df, y = if_else(detected, log2FC, nd_y))

fig5B <- ggplot(pmo_df, aes(CH4_ppm, y, colour = Site)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = nd_y + 0.75, colour = "grey80") +          # separates n.d. row
  geom_point(aes(shape = detected), size = 3.5, stroke = 1.2) +
  geom_text_repel(aes(label = Site), size = 10 / .pt, show.legend = FALSE,
                  box.padding = 0.4, min.segment.length = 0) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = site_cols, guide = "none") +
  scale_x_log10() +
  scale_y_continuous(breaks = c(nd_y, y_breaks),
                     labels = c("n.d.", y_breaks)) +
  labs(x = expression(CH[4]*" (ppm)"),
       y = expression(log[2]*"(aerosol / paired environment)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.6),
        axis.text  = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 11))

# Small size for the inset next to the hydrogen metabolism rows of Fig. 5A
ggsave(file.path(OUT_DIR, "Fig5B_pmoCAB_vs_CH4.pdf"), plot = fig5B,
       width = 80, height = 100, units = "mm")


# =============================================================================
# FIGURE 4: differential genus abundance, aerosol vs paired environment
#   - Site-blocked DESeq2 (~ Site + Compartment) on NON-rarefied genus counts
#   - Genera with < 10 total counts dropped; size factors by "poscounts"
#   - apeglm shrinkage; Benjamini-Hochberg adjustment
#   - Significant: padj < 0.05 AND |log2FC| >= 2; top 15 by padj labelled
#   - Requires Section 0; DESeq2 and apeglm must be installed
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tibble)
})


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
bac_df  <- get_deseq_results(file.path(PS_DIR, "psf_bacteria_nonrarefied.rds"))
p_bac   <- make_volcano(bac_df, "Bacteria")

# =========================================================
# RUN EUKARYOTA
# =========================================================

cat("Running Eukaryota...\n")
euk_df  <- get_deseq_results(file.path(PS_DIR, "psf_eukaryota_nonrarefied.rds"))
p_euk   <- make_volcano(euk_df, "Eukaryota")

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


# Save
pdf(file.path(OUT_DIR, "Fig4_volcano_combined.pdf"), width = 12, height = 6)
print(combined)
dev.off()

cat("Saved: Fig4_volcano_combined.pdf\n")

# Also save individually
pdf(file.path(OUT_DIR, "Fig4A_bacteria_volcano.pdf"), width = 5, height = 6)
print(p_bac)
dev.off()

pdf(file.path(OUT_DIR, "Fig4B_eukaryota_volcano.pdf"), width = 5, height = 6)
print(p_euk)
dev.off()

cat("Saved individual PDFs.\n")


# =============================================================================
# SUPPLEMENTARY FIG. 5: methanotroph genera in read-based taxonomic profiles
#   - Supports the Results statement that pmoCAB detection in aerosols at VL
#     and DF coincides with canonical methanotrophic genera in the Bracken
#     profiles at those sites
#   - Genus list is defined a priori from the described aerobic methanotrophs
#     (Alphaproteobacteria Type II, Methylococcales Type I and X, the post-2018
#     Methylotuvi* splits, Verrucomicrobiota, and NC10 Methylomirabilis)
#   - Relative abundance is computed against the WHOLE bacterial community
#     before subsetting, so values are comparable across samples
#   - Requires Section 0
# =============================================================================

methanotroph_genera <- c(
  # Alphaproteobacteria (Type II)
  "Methylocystis", "Methylosinus", "Methylocella", "Methyloferula", "Methylocapsa",
  # Gammaproteobacteria / Methylococcales (Type I and X)
  "Methylomonas", "Methylomicrobium", "Methylosoma", "Methylosphaera",
  "Methylococcus", "Methylocaldum", "Methylothermus", "Methylovulum",
  "Methyloprofundus", "Methylomarinum", "Crenothrix",
  # Post-2018 genus splits within Methylococcaceae
  "Methylotuvimicrobium", "Methylotuvicoccus", "Methylotetracoccus",
  "Methylotuvibacter",
  # Verrucomicrobiota (acidophilic methanotrophs)
  "Methylacidiphilum", "Methylacidimicrobium",
  # NC10 (nitrite-dependent anaerobic methane oxidation)
  "Candidatus Methylomirabilis"
)

ps_bac_all  <- readRDS(file.path(PS_DIR, "psf_bacteria_nonrarefied.rds"))
ps_rel_all  <- transform_sample_counts(ps_bac_all, function(x) x / sum(x))
ps_meth_rel <- subset_taxa(ps_rel_all, genus %in% methanotroph_genera)

if (ntaxa(ps_meth_rel) == 0) {
  warning("No methanotroph genera detected; Supplementary Fig. 5 not written.")
} else {

  df_meth <- psmelt(ps_meth_rel)
  if ("sample_Sample" %in% names(df_meth)) df_meth$Sample <- df_meth$sample_Sample

  df_meth$genus <- as.character(df_meth$genus)
  df_meth$genus[df_meth$genus == "" | is.na(df_meth$genus)] <- "Unclassified"

  # Order samples aerosol-first within each site, sites high -> low CH4
  df_meth <- df_meth %>%
    mutate(Site   = factor(sub("^[AWS]_([A-Z]+).*", "\\1", Sample),
                           levels = site_levels),
           is_air = grepl("^A_", Sample)) %>%
    group_by(Site) %>%
    arrange(desc(is_air), Sample, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))

  g_levels <- sort(unique(df_meth$genus))
  pal <- setNames(
    colorRampPalette(c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1",
                       "#E15759", "#76B7B2", "#EDC948"))(length(g_levels)),
    g_levels
  )

  df_meth_plot <- df_meth %>%
    group_by(Sample, Site, Sample_type, genus) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

  supp4A <- ggplot(df_meth_plot, aes(Sample, Abundance, fill = genus)) +
    geom_col(width = 0.9, colour = NA) +
    facet_grid(~ Site, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, name = "Genus") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.02)),
                       breaks = scales::breaks_extended(n = 6)) +
    labs(x = NULL, y = "Relative abundance of bacterial community") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text = element_text(face = "bold", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                     colour = "black")) +
    guides(fill = guide_legend(ncol = 1))

  supp4B <- ggplot(df_meth, aes(Sample, genus, fill = Abundance)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(~ Site, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(option = "magma", trans = "sqrt",
                         labels = scales::percent_format(accuracy = 0.01),
                         name = "Relative\nabundance") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text = element_text(face = "bold", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                     colour = "black"),
          axis.text.y = element_text(face = "italic", size = 10,
                                     colour = "black"))

  ggsave(file.path(OUT_DIR, "SuppFig5_methanotroph_genera_barplot.pdf"),
         plot = supp4A, width = 12, height = 7)
  ggsave(file.path(OUT_DIR, "SuppFig5_methanotroph_genera_heatmap.pdf"),
         plot = supp4B, width = 12, height = 7)

  # Per-site totals quoted in the Results
  df_meth_plot %>%
    group_by(Site, Sample) %>%
    summarise(total_methanotroph = sum(Abundance), .groups = "drop") %>%
    arrange(Site, Sample) %>%
    print(n = Inf)
}

# =============================================================================
# SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo_02.txt"))
cat("\nDone. Figures and tables written to", OUT_DIR, "\n")
