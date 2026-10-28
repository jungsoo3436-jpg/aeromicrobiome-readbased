###############################################################################
# Aerosol metagenome manuscript (read-based community analysis)
# Step 03: Curated marker gene panels for the supplementary functional figures
#
# Author: Jungsoo Park (Simon Fraser University)
#
# Produces the five supplementary marker gene heatmaps (Supplementary Figs 6-10):
#   Supp_A_VOC_carbon.pdf   VOC and complex carbon degradation
#   Supp_B_Nitrogen.pdf     Nitrogen cycling
#   Supp_C_Phosphorus.pdf   Phosphorus cycling
#   Supp_D_Sulfur.pdf       Sulfur cycling
#   Supp_E_Survival.pdf     Survival and persistence traits
#
# Input  : data/gene_catalog/KOxSample_TPM.prop.tsv   KO x sample proportional TPM
#          data/gene_catalog/KO_to_Function.clean3.tsv  KO -> gene name, description
# Output : figures/
#
# Marker genes for every category are defined a priori in the tribble() blocks
# below and are the source of Supplementary Table S1. KOs absent from every
# sample are dropped before plotting but remain listed in that table.
#
# Abundances are proportional TPM. Tiles are within-site z-scores, computed per
# gene across the three samples of a site, for visualisation only; grey tiles
# mark non-detection. Samples are ordered aerosol first within each site, and
# sites run from high to low ambient CH4 (VL, DF, MK, BB, RM).
#
# The C1 and trace gas panel (main Figure 5A) is produced by
# 02_community_and_methane_figures.R, not here.
#
# Run from a FRESH R session, from the repository root.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggh4x)
})

KO_DIR  <- "data/gene_catalog"
OUT_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# LOAD KO x SAMPLE MATRIX (proportional TPM)
# =========================================================
ko_mat <- read.delim(
  file.path(KO_DIR, "KOxSample_TPM.prop.tsv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ko_mat[] <- lapply(ko_mat, function(x) as.numeric(as.character(x)))

ko_func <- read.delim(
  file.path(KO_DIR, "KO_to_Function.clean3.tsv"),
  header = FALSE,
  sep = "\t",
  quote = "\"",
  stringsAsFactors = FALSE,
  col.names = c("KO", "Preferred", "Description")
)

# Restrict to KOs with a functional annotation (prokaryotic scope)
ko_mat <- ko_mat[rownames(ko_mat) %in% ko_func$KO, ]
stopifnot(all(rownames(ko_mat) %in% ko_func$KO))

cat("KO matrix:", nrow(ko_mat), "KOs x", ncol(ko_mat), "samples\n")

# =========================================================
# SHARED HELPERS
# =========================================================

infer_site <- function(s) {
  case_when(
    str_detect(s, "_VL") ~ "VL",
    str_detect(s, "_DF") ~ "DF",
    str_detect(s, "_MK") ~ "MK",
    str_detect(s, "_BB") ~ "BB",
    str_detect(s, "_RM") ~ "RM",
    TRUE ~ NA_character_
  )
}

site_order <- c("VL","DF","MK","BB","RM")

make_sample_order <- function(ko_long) {
  ko_long %>%
    distinct(Site, Sample) %>%
    group_by(Site) %>%
    summarise(levels = list({
      s <- unique(Sample)
      air <- s[grepl("^A_", s)]
      other <- sort(s[!grepl("^A_", s)])
      c(air, other)
    }), .groups = "drop") %>%
    unnest(levels) %>%
    group_by(Site) %>%
    mutate(Sample_in_site = factor(levels, levels = unique(levels))) %>%
    rename(Sample = levels) %>%
    ungroup() %>%
    select(Site, Sample, Sample_in_site)
}

build_gene_long <- function(ko_map, ko_mat, ko_long, sample_order_tbl,
                            category_levels, gene_levels,
                            extra_rows = NULL) {
  
  # ----------------------------------------------------------
  # Key fix: when same KO appears in multiple genes/categories
  # (e.g. narG and nxrA both mapped to K00370),
  # we join by KO then filter to keep only the intended gene
  # ----------------------------------------------------------
  ko_map_clean <- ko_map %>%
    filter(KO %in% rownames(ko_mat)) %>%
    distinct(Category, Gene, KO)
  
  gene_long_raw <- ko_map_clean %>%
    left_join(ko_long, by = "KO", relationship = "many-to-many") %>%
    group_by(Category, Gene, Site, Sample) %>%
    summarise(prop = sum(prop, na.rm = TRUE), .groups = "drop")
  
  all_rows <- bind_rows(
    ko_map %>% distinct(Category, Gene),
    extra_rows
  ) %>% distinct()
  
  expand_grid(
    Category = all_rows$Category,
    Gene     = all_rows$Gene,
    Site     = site_order,
    Sample   = ko_long %>% distinct(Site, Sample) %>% pull(Sample) %>% unique()
  ) %>%
    mutate(Site = infer_site(Sample)) %>%
    filter(!is.na(Site)) %>%
    mutate(Site = factor(Site, levels = site_order)) %>%
    left_join(gene_long_raw, by = c("Category","Gene","Site","Sample")) %>%
    mutate(prop = replace_na(prop, 0)) %>%
    left_join(sample_order_tbl, by = c("Site","Sample")) %>%
    mutate(
      Category  = factor(Category, levels = category_levels),
      Gene      = factor(Gene, levels = gene_levels),
      GeneLabel = as.character(Gene)
    ) %>%
    filter(Sample %in% colnames(ko_mat)) %>%
    { . ->> tmp
      tmp %>%
        group_by(Category, Gene) %>%
        summarise(any_present = any(prop > 0), .groups = "drop") %>%
        filter(any_present) %>%
        semi_join(tmp, ., by = c("Category","Gene"))
    } %>%
    group_by(Site, Gene) %>%
    mutate(
      z      = if (sd(prop, na.rm = TRUE) == 0) 0
      else (prop - mean(prop, na.rm = TRUE)) / sd(prop, na.rm = TRUE),
      z_plot = ifelse(prop == 0, NA, z)
    ) %>%
    ungroup()
}

make_heatmap <- function(gene_long, gene_levels_present, title_label,
                         vmax = NULL, base_size = 9) {
  if (is.null(vmax)) {
    vmax <- max(abs(gene_long$z_plot), na.rm = TRUE)
    if (!is.finite(vmax) || vmax == 0) vmax <- 1
  }
  
  present_genes <- gene_long %>%
    filter(!is.na(z_plot)) %>%
    pull(GeneLabel) %>%
    unique()
  
  ordered_genes <- rev(gene_levels_present[gene_levels_present %in% present_genes])
  gene_long$GeneLabel <- factor(gene_long$GeneLabel, levels = ordered_genes)
  gene_long <- gene_long %>% filter(GeneLabel %in% ordered_genes)
  
  ggplot(gene_long, aes(x = Sample_in_site, y = GeneLabel, fill = z_plot)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-vmax, vmax),
      name     = "z-score\n(within site)",
      na.value = "grey85",
      guide    = guide_colorbar(
        barwidth       = 6,
        barheight      = 0.5,
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    facet_grid(
      Category ~ Site,
      scales = "free",
      space  = "free"
    ) +
    labs(title = title_label) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title         = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                        size = base_size - 1, color = "grey20"),
      axis.text.y        = element_text(size = base_size - 1, face = "italic",
                                        color = "grey10"),
      axis.ticks         = element_line(color = "grey70", linewidth = 0.3),
      panel.grid         = element_blank(),
      panel.border       = element_rect(color = "grey80", linewidth = 0.4),
      strip.text.x       = element_text(face = "bold", size = base_size,
                                        color = "grey10"),
      strip.text.y       = element_text(face = "bold", size = base_size - 2,
                                        angle = 0, hjust = 0, color = "grey10"),
      strip.background.x = element_rect(fill = "grey95", color = "grey70"),
      strip.background.y = element_rect(fill = "grey95", color = "grey70"),
      legend.position    = "bottom",
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 2),
      plot.background    = element_rect(fill = "white", color = NA)
    )
}

# =========================================================
# LOAD DATA
# =========================================================
ko_long <- ko_mat %>%
  as.data.frame() %>%
  rownames_to_column("KO") %>%
  pivot_longer(-KO, names_to = "Sample", values_to = "prop") %>%
  mutate(prop = as.numeric(prop), Site = infer_site(Sample)) %>%
  filter(!is.na(Site)) %>%
  mutate(Site = factor(Site, levels = site_order))

sample_order_tbl <- make_sample_order(ko_long)

# =========================================================
# PANEL A — VOC & complex carbon degradation
# =========================================================
voc_ko_map <- tribble(
  ~Category, ~Gene, ~KO,
  "Alcohol & aldehyde oxidation",    "adh",    "K13953",
  "Alcohol & aldehyde oxidation",    "aldh",   "K00128",
  "Organic acid metabolism",         "acs",    "K01895",
  "Organic acid metabolism",         "ldh",    "K00016",
  "Organic acid metabolism",         "poxB",   "K00156",
  "Aromatic compound degradation",   "benA",   "K04110",
  "Aromatic compound degradation",   "benB",   "K04111",
  "Aromatic compound degradation",   "benC",   "K04112",
  "Aromatic compound degradation",   "catA",   "K03381",
  "Aromatic compound degradation",   "catB",   "K03382",
  "Aromatic compound degradation",   "catC",   "K03383",
  "Aromatic compound degradation",   "pcaG",   "K00448",
  "Aromatic compound degradation",   "pcaH",   "K00449",
  "Alkane & hydrocarbon degradation","alkB",   "K00496",
  "Alkane & hydrocarbon degradation","ladA",   "K16244",
  "Alkane & hydrocarbon degradation","cyp153", "K15446",
  "Terpene degradation",             "atuA",   "K06443",
  "Terpene degradation",             "atuB",   "K06444",
  "Terpene degradation",             "limA",   "K07539",
  "Carbohydrate degradation",        "celA",   "K01179",
  "Carbohydrate degradation",        "xynA",   "K01181",
  "Carbohydrate degradation",        "amyA",   "K01176"
) %>% distinct(Category, Gene, KO)

voc_cat_levels <- c(
  "Alcohol & aldehyde oxidation",
  "Organic acid metabolism",
  "Aromatic compound degradation",
  "Alkane & hydrocarbon degradation",
  "Terpene degradation",
  "Carbohydrate degradation"
)
voc_gene_levels <- c(
  "adh","aldh",
  "acs","ldh","poxB",
  "benA","benB","benC","catA","catB","catC","pcaG","pcaH",
  "alkB","ladA","cyp153",
  "atuA","atuB","limA",
  "celA","xynA","amyA"
)

gl_voc <- build_gene_long(voc_ko_map, ko_mat, ko_long, sample_order_tbl,
                          voc_cat_levels, voc_gene_levels)
p_A <- make_heatmap(gl_voc, voc_gene_levels, "(A) VOC & complex carbon degradation")

# =========================================================
# PANEL B — Nitrogen cycling
# Fix: narG/narH deliberately share KOs with nxrA/nxrB
# We keep both as separate gene labels — abundance will be
# identical but they represent different organisms/contexts
# =========================================================
n_ko_map <- tribble(
  ~Category, ~Gene, ~KO,
  "Nitrogen fixation",               "nifH",  "K02588",
  "Nitrogen fixation",               "nifD",  "K02586",
  "Nitrogen fixation",               "nifK",  "K02591",
  "Ammonia oxidation",               "amoA",  "K10944",
  "Ammonia oxidation",               "amoB",  "K10945",
  "Ammonia oxidation",               "amoC",  "K10946",
  "Ammonia oxidation",               "hao",   "K10535",
  "Nitrite oxidation",               "nxrA",  "K00370",
  "Nitrite oxidation",               "nxrB",  "K00371",
  "Denitrification",                 "napA",  "K02567",
  "Denitrification",                 "napB",  "K02568",
  "Denitrification",                 "narI",  "K00374",
  "Denitrification",                 "nirK",  "K00368",
  "Denitrification",                 "nirS",  "K15864",
  "Denitrification",                 "norB",  "K04561",
  "Denitrification",                 "norC",  "K02305",
  "Denitrification",                 "nosZ",  "K00376",
  "DNRA",                            "nrfA",  "K03385",
  "DNRA",                            "nrfH",  "K15876",
  "Assimilatory nitrate reduction",  "nasA",  "K00372",
  "Assimilatory nitrate reduction",  "nirA",  "K00366",
  "Anaerobic ammonium oxidation",    "hzsA",  "K20932",
  "Anaerobic ammonium oxidation",    "hzsB",  "K20933",
  "Anaerobic ammonium oxidation",    "hzsC",  "K20934",
  "Anaerobic ammonium oxidation",    "hzo",   "K20935",
  "Urea hydrolysis & transport",     "ureA",  "K01428",
  "Urea hydrolysis & transport",     "ureB",  "K01429",
  "Urea hydrolysis & transport",     "ureC",  "K01430",
  "Urea hydrolysis & transport",     "urtA",  "K11959",
  "Urea hydrolysis & transport",     "urtB",  "K11960",
  "Urea hydrolysis & transport",     "urtC",  "K11961",
  "Urea hydrolysis & transport",     "urtD",  "K11962",
  "Urea hydrolysis & transport",     "urtE",  "K11963"
) %>% distinct(Category, Gene, KO)

n_cat_levels <- c(
  "Nitrogen fixation",
  "Ammonia oxidation",
  "Nitrite oxidation",
  "Denitrification",
  "DNRA",
  "Assimilatory nitrate reduction",
  "Anaerobic ammonium oxidation",
  "Urea hydrolysis & transport"
)
n_gene_levels <- c(
  "nifH","nifD","nifK",
  "amoA","amoB","amoC","hao",
  "nxrA","nxrB",
  "napA","napB","narI","nirK","nirS","norB","norC","nosZ",
  "nrfA","nrfH",
  "nasA","nirA",
  "hzsA","hzsB","hzsC","hzo",
  "ureA","ureB","ureC","urtA","urtB","urtC","urtD","urtE"
)

gl_n <- build_gene_long(n_ko_map, ko_mat, ko_long, sample_order_tbl,
                        n_cat_levels, n_gene_levels)
p_B <- make_heatmap(gl_n, n_gene_levels, "(B) Nitrogen cycling")

# =========================================================
# PANEL C — Phosphorus cycling
# =========================================================
p_ko_map <- tribble(
  ~Category, ~Gene, ~KO,
  "Phosphate uptake (Pst/Pit)",      "pstS",  "K02040",
  "Phosphate uptake (Pst/Pit)",      "pstC",  "K02038",
  "Phosphate uptake (Pst/Pit)",      "pstA",  "K02037",
  "Phosphate uptake (Pst/Pit)",      "pstB",  "K02036",
  "Phosphate uptake (Pst/Pit)",      "pitA",  "K01083",
  "Phosphate sensing (Pho regulon)", "phoB",  "K07657",
  "Phosphate sensing (Pho regulon)", "phoR",  "K07636",
  "Organic P mineralization",        "phoA",  "K01077",
  "Organic P mineralization",        "phoD",  "K01113",
  "Organic P mineralization",        "phoX",  "K01114",
  "Phosphonate utilization",         "phnC",  "K02044",
  "Phosphonate utilization",         "phnD",  "K02043",
  "Phosphonate utilization",         "phnE",  "K02042",
  "Phosphonate utilization",         "phnG",  "K06163",
  "Phosphonate utilization",         "phnH",  "K06164",
  "Phosphonate utilization",         "phnI",  "K06165",
  "Phosphonate utilization",         "phnJ",  "K06166",
  "Phosphonate utilization",         "phnK",  "K06167",
  "Phosphonate utilization",         "phnL",  "K06168",
  "Phosphonate utilization",         "phnM",  "K06169",
  "Polyphosphate metabolism",        "ppk1",  "K00937",
  "Polyphosphate metabolism",        "ppk2",  "K01093",
  "Polyphosphate metabolism",        "ppx",   "K01095"
) %>% distinct(Category, Gene, KO)

p_cat_levels <- c(
  "Phosphate uptake (Pst/Pit)",
  "Phosphate sensing (Pho regulon)",
  "Organic P mineralization",
  "Phosphonate utilization",
  "Polyphosphate metabolism"
)
p_gene_levels <- c(
  "pstS","pstC","pstA","pstB","pitA",
  "phoB","phoR",
  "phoA","phoD","phoX",
  "phnC","phnD","phnE","phnG","phnH","phnI","phnJ","phnK","phnL","phnM",
  "ppk1","ppk2","ppx"
)

gl_p <- build_gene_long(p_ko_map, ko_mat, ko_long, sample_order_tbl,
                        p_cat_levels, p_gene_levels)
p_C <- make_heatmap(gl_p, p_gene_levels, "(C) Phosphorus cycling")

# =========================================================
# PANEL D — Sulfur cycling (with oxidative Dsr genes added)
# =========================================================
s_ko_map <- tribble(
  ~Category, ~Gene, ~KO,
  "Dissimilatory sulfate reduction", "sat",   "K00958",
  "Dissimilatory sulfate reduction", "aprA",  "K00394",
  "Dissimilatory sulfate reduction", "aprB",  "K00395",
  "Dissimilatory sulfate reduction", "dsrA",  "K11180",
  "Dissimilatory sulfate reduction", "dsrB",  "K11181",
  "Sulfur oxidation",                "soxA",  "K17222",
  "Sulfur oxidation",                "soxX",  "K17223",
  "Sulfur oxidation",                "soxY",  "K17224",
  "Sulfur oxidation",                "soxZ",  "K17225",
  "Sulfur oxidation",                "soxB",  "K17226",
  "Sulfur oxidation",                "soxC",  "K17227",
  "Sulfur oxidation",                "soxD",  "K17228",
  "Sulfur oxidation",                "sqr",   "K17218",
  "Sulfur oxidation",                "fccA",  "K15228",
  "Sulfur oxidation",                "fccB",  "K15229",
  "Sulfur oxidation",                "soeA",  "K21307",
  "Sulfur oxidation",                "soeB",  "K21308",
  "Sulfur oxidation",                "soeC",  "K21309",
  "Sulfur oxidation",                "dsrC",  "K17217",  # oxidative Dsr
  "Sulfur oxidation",                "dsrE",  "K17219",  # oxidative Dsr
  "Sulfur oxidation",                "dsrF",  "K17220",  # oxidative Dsr
  "Sulfur oxidation",                "dsrH",  "K17221",  # oxidative Dsr
  "Thiosulfate disproportionation",  "phsA",  "K08352",
  "Thiosulfate disproportionation",  "phsB",  "K08353",
  "Thiosulfate disproportionation",  "phsC",  "K08354",
  "Assimilatory sulfate reduction",  "cysD",  "K00955",
  "Assimilatory sulfate reduction",  "cysN",  "K00956",
  "Assimilatory sulfate reduction",  "cysC",  "K00957",
  "Assimilatory sulfate reduction",  "cysH",  "K00392",
  "Assimilatory sulfate reduction",  "cysI",  "K00380",
  "Assimilatory sulfate reduction",  "cysJ",  "K00381",
  "Organosulfur utilization",        "tauA",  "K07239",
  "Organosulfur utilization",        "tauB",  "K07240",
  "Organosulfur utilization",        "tauC",  "K07241",
  "Organosulfur utilization",        "tauD",  "K03118",
  "Organosulfur utilization",        "ssuA",  "K12339",
  "Organosulfur utilization",        "ssuB",  "K12340",
  "Organosulfur utilization",        "ssuC",  "K12341",
  "Organosulfur utilization",        "ssuD",  "K12342",
  "Organosulfur utilization",        "ssuE",  "K12343",
  "Organosulfur utilization",        "msuD",  "K18275",
  "Organosulfur utilization",        "dmdA",  "K14468",
  "Sulfur storage & transfer",       "sseA",  "K00658",
  "Sulfur storage & transfer",       "glpE",  "K01011"
) %>% distinct(Category, Gene, KO)

s_cat_levels <- c(
  "Dissimilatory sulfate reduction",
  "Sulfur oxidation",
  "Thiosulfate disproportionation",
  "Assimilatory sulfate reduction",
  "Organosulfur utilization",
  "Sulfur storage & transfer"
)
s_gene_levels <- c(
  "sat","aprA","aprB","dsrA","dsrB",
  "soxA","soxX","soxY","soxZ","soxB","soxC","soxD","sqr","fccA","fccB",
  "soeA","soeB","soeC","dsrC","dsrE","dsrF","dsrH",
  "phsA","phsB","phsC",
  "cysD","cysN","cysC","cysH","cysI","cysJ",
  "tauA","tauB","tauC","tauD","ssuA","ssuB","ssuC","ssuD","ssuE","msuD","dmdA",
  "sseA","glpE"
)

gl_s <- build_gene_long(s_ko_map, ko_mat, ko_long, sample_order_tbl,
                        s_cat_levels, s_gene_levels)
p_D <- make_heatmap(gl_s, s_gene_levels, "(D) Sulfur cycling")

# =========================================================
# PANEL E — Survival & persistence traits
# =========================================================
surv_ko_map <- tribble(
  ~Category, ~Gene, ~KO,
  "DNA repair & genome stability",        "recA",      "K03553",
  "DNA repair & genome stability",        "uvrA",      "K03701",
  "DNA repair & genome stability",        "uvrB",      "K03702",
  "DNA repair & genome stability",        "uvrC",      "K03703",
  "DNA repair & genome stability",        "mutS",      "K02335",
  "DNA repair & genome stability",        "uvrD",      "K03555",
  "Oxidative stress defense",             "sodA",      "K04564",
  "Oxidative stress defense",             "katG",      "K03781",
  "Oxidative stress defense",             "ahpC",      "K03386",
  "Oxidative stress defense",             "trxB",      "K03671",
  "Carbon & energy storage",              "ppk1",      "K00937",
  "Carbon & energy storage",              "ppk2",      "K01093",
  "Carbon & energy storage",              "ppx",       "K01524",
  "Carbon & energy storage",              "phaA",      "K00626",
  "Carbon & energy storage",              "phaB",      "K00025",
  "Carbon & energy storage",              "phaC",      "K03821",
  "Sporulation & dormancy",               "spo0A",     "K07699",
  "UV protection & photolyase",           "crtB",      "K02291",
  "UV protection & photolyase",           "phr",       "K01669",
  "Biofilm & EPS production",             "bcsA",      "K00694",
  "Desiccation & osmotic stress",         "otsA",      "K00697",
  "Acid stress tolerance",                "gadA_gadB", "K01580",
  "Acid stress tolerance",                "gadC",      "K01581",
  "Mechanosensing & membrane adaptation", "mscL",      "K03823",
  "Mechanosensing & membrane adaptation", "mscS",      "K03824",
  "Trace gas oxidation",                  "coxL",      "K03518",
  "Trace gas oxidation",                  "coxM",      "K03519",
  "Trace gas oxidation",                  "coxS",      "K03520"
) %>% distinct(Category, Gene, KO)

surv_cat_levels <- c(
  "DNA repair & genome stability",
  "Oxidative stress defense",
  "Carbon & energy storage",
  "Sporulation & dormancy",
  "UV protection & photolyase",
  "Biofilm & EPS production",
  "Desiccation & osmotic stress",
  "Acid stress tolerance",
  "Mechanosensing & membrane adaptation",
  "Trace gas oxidation"
)
surv_gene_levels <- c(
  "recA","uvrA","uvrB","uvrC","mutS","uvrD",
  "sodA","katG","ahpC","trxB",
  "ppk1","ppk2","ppx","phaA","phaB","phaC",
  "spo0A",
  "crtB","phr",
  "bcsA",
  "otsA",
  "gadA_gadB","gadC",
  "mscL","mscS",
  "coxL","coxM","coxS"
)

gl_surv <- build_gene_long(surv_ko_map, ko_mat, ko_long, sample_order_tbl,
                           surv_cat_levels, surv_gene_levels)
p_E <- make_heatmap(gl_surv, surv_gene_levels, "(E) Survival & persistence traits")

# =========================================================
# COMBINE & SAVE
# =========================================================
combined <- (p_A + p_B) / (p_C + p_D) / p_E +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Supplementary Figure | Functional gene potential across aerosol and paired environmental samples",
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )

# =========================================================
# SAVE — separately for manual PPT assembly
# =========================================================

pdf(file.path(OUT_DIR, "Supp_A_VOC_carbon.pdf"), width = 6, height = 5)
print(p_A)
dev.off()

pdf(file.path(OUT_DIR, "Supp_B_Nitrogen.pdf"), width = 6, height = 6)
print(p_B)
dev.off()

pdf(file.path(OUT_DIR, "Supp_C_Phosphorus.pdf"), width = 5, height = 6)
print(p_C)
dev.off()

pdf(file.path(OUT_DIR, "Supp_D_Sulfur.pdf"), width = 5, height = 5)
print(p_D)
dev.off()

pdf(file.path(OUT_DIR, "Supp_E_Survival.pdf"), width = 5, height = 6)
print(p_E)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_03.txt"))
