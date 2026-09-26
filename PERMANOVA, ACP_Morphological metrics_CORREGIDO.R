# ============================================================
# MORPHOLOGICAL ANALYSIS
# Stegastes acapulcoensis vs S. flavilatus
#
# Variables expressed as SL / morphological measurement
# PCA performed on centered, unscaled variables
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(grid)


# ============================================================
# 1. IMPORT DATA
# ============================================================

morph <- read.csv(
  "Morphological metrics.csv",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# Inspect data
# ------------------------------------------------------------

str(morph)
head(morph)
table(morph$species)


# ------------------------------------------------------------
# Keep variables used in the analysis
# ------------------------------------------------------------

morph <- morph %>%
  select(
    individual,
    species,
    ED,
    PO,
    VCP,
    DCF,
    DL,
    HD,
    BW
  ) %>%
  mutate(
    species = factor(species)
  )


# Morphological variables
morph_vars <- c(
  "ED",
  "PO",
  "VCP",
  "DCF",
  "DL",
  "HD",
  "BW"
)


# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------

# Number of individuals per species
table(morph$species)

# Missing values
colSums(is.na(morph[, morph_vars]))

# Descriptive summary
summary(morph[, morph_vars])


# ============================================================
# 2. PERMANOVA
# Unscaled morphological proportions
# Euclidean distance
# ============================================================

morph_raw <- morph[, morph_vars]


# Euclidean distance matrix
d_morph_raw <- dist(
  morph_raw,
  method = "euclidean"
)


# PERMANOVA
set.seed(123)

perm_morph_raw <- adonis2(
  d_morph_raw ~ species,
  data = morph,
  permutations = 9999
)

print(perm_morph_raw)


# ============================================================
# 3. PERMDISP
# Homogeneity of multivariate dispersion
# ============================================================

disp_morph_raw <- betadisper(
  d_morph_raw,
  group = morph$species
)


set.seed(123)

permdisp_morph_raw <- permutest(
  disp_morph_raw,
  permutations = 9999
)

print(permdisp_morph_raw)


# ============================================================
# 4. PRINCIPAL COMPONENT ANALYSIS (PCA)
# Centered but NOT scaled
# ============================================================

pca_morph_raw <- prcomp(
  morph[, morph_vars],
  center = TRUE,
  scale. = FALSE
)


# PCA summary
summary(pca_morph_raw)


# ============================================================
# 5. PCA VARIANCE EXPLAINED
# ============================================================

eigenvalues_raw <- pca_morph_raw$sdev^2


variance_raw <- (
  eigenvalues_raw /
    sum(eigenvalues_raw)
) * 100


pca_variance_raw <- data.frame(
  PC = paste0(
    "PC",
    seq_along(eigenvalues_raw)
  ),
  Eigenvalue = eigenvalues_raw,
  Variance_percent = variance_raw,
  Cumulative_percent = cumsum(variance_raw)
)


print(pca_variance_raw)


# Variance explained by PC1 and PC2
var_PC1_raw <- variance_raw[1]
var_PC2_raw <- variance_raw[2]


cat(
  "\nPC1 =",
  round(var_PC1_raw, 1),
  "%\n"
)

cat(
  "PC2 =",
  round(var_PC2_raw, 1),
  "%\n"
)

cat(
  "PC1 + PC2 =",
  round(var_PC1_raw + var_PC2_raw, 1),
  "%\n"
)


# ============================================================
# 6. PCA SCORES
# ============================================================

scores_raw <- as.data.frame(
  pca_morph_raw$x
) %>%
  mutate(
    species = morph$species,
    individual = morph$individual
  )


head(scores_raw)


# ============================================================
# 7. PCA LOADINGS
# ============================================================

loadings_raw <- as.data.frame(
  pca_morph_raw$rotation
) %>%
  rownames_to_column(
    "Variable"
  )


print(loadings_raw)


# ------------------------------------------------------------
# Loadings for PC1 and PC2 only
# ------------------------------------------------------------

loadings_PC12_raw <- loadings_raw %>%
  select(
    Variable,
    PC1,
    PC2
  )


print(loadings_PC12_raw)


# ============================================================
# 8. CONVEX HULLS
# ============================================================

hulls_raw <- scores_raw %>%
  group_by(species) %>%
  slice(
    chull(
      PC1,
      PC2
    )
  ) %>%
  ungroup()


# ============================================================
# 9. LOADINGS FOR MORPHOSPACE
# ============================================================

# Multiplication factor ONLY for graphical representation
# This does not modify PCA loadings or PCA results.

arrow_mult_raw <- 4


loadings_plot_raw <- loadings_PC12_raw %>%
  mutate(
    PC1_plot = PC1 * arrow_mult_raw,
    PC2_plot = PC2 * arrow_mult_raw
  )


# ============================================================
# 10. MORPHOSPACE
# Convex hulls + morphological vectors
# ============================================================


# ------------------------------------------------------------
# Species colors
# ------------------------------------------------------------

col_species <- c(
  "S. acapulcoensis" = "#C16540",
  "S. flavilatus"    = "#00C1C8"
)


# ------------------------------------------------------------
# Species names in italics
# ------------------------------------------------------------

long_labels <- c(
  expression(
    italic("Stegastes acapulcoensis")
  ),
  expression(
    italic("Stegastes flavilatus")
  )
)


# ------------------------------------------------------------
# Morphospace
# ------------------------------------------------------------

pca_hull_biplot <- ggplot(
  scores_raw,
  aes(
    x = PC1,
    y = PC2,
    color = species
  )
) +
  
  # Convex hulls
  geom_polygon(
    data = hulls_raw,
    aes(
      x = PC1,
      y = PC2,
      group = species,
      fill = species,
      color = species
    ),
    alpha = 0.18,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Individuals
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  # Horizontal zero line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  
  # Vertical zero line
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  
  # Morphological vectors
  geom_segment(
    data = loadings_plot_raw,
    aes(
      x = 0,
      y = 0,
      xend = PC1_plot,
      yend = PC2_plot
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = grid::unit(
        0.18,
        "cm"
      )
    ),
    linewidth = 0.65,
    color = "black"
  ) +
  
  # Variable labels
  ggrepel::geom_text_repel(
    data = loadings_plot_raw,
    aes(
      x = PC1_plot,
      y = PC2_plot,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    family = "Arial",
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 0.4,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  # Species colors
  scale_color_manual(
    values = col_species,
    labels = long_labels
  ) +
  
  scale_fill_manual(
    values = col_species,
    guide = "none"
  ) +
  
  # Axis labels
  labs(
    title = NULL,
    x = paste0(
      "PC1 (",
      round(var_PC1_raw, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_PC2_raw, 1),
      "%)"
    ),
    color = NULL
  ) +
  
  # Theme
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    text = element_text(
      size = 11,
      family = "Arial"
    ),
    
    axis.title = element_text(
      size = 11,
      family = "Arial"
    ),
    
    axis.text = element_text(
      size = 11,
      family = "Arial",
      color = "black"
    ),
    
    legend.position = "top",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 11,
      family = "Arial"
    )
  )


# Display morphospace
print(pca_hull_biplot)

############################################################
# 26) EXPORT HIGH-RESOLUTION FIGURE
############################################################

ggsave(
  filename = paste0(
    "F:/Doctorado/TESIS version 3/",
    "capitulo 2.- corregido/figuras/PCA.tiff"
  ),
  plot = pca_hull_biplot,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6,
  height = 6,
  units = "in"
)


# ============================================================
# 11. PCA LOADINGS — NUMERICAL VALUES
# ============================================================

print(
  loadings_PC12_raw
)

# ============================================================
# 12. GRAPH OF PCA LOADINGS
# PC1 and PC2
# ============================================================

# ------------------------------------------------------------
# Prepare loadings
# ------------------------------------------------------------

loadings_long_raw <- loadings_PC12_raw %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "Loading"
  )


# ------------------------------------------------------------
# Color
# ------------------------------------------------------------

col_loadings <- "#006D61"


# ------------------------------------------------------------
# Loadings plot
# ------------------------------------------------------------

p_loadings <- ggplot(
  loadings_long_raw,
  aes(
    x = reorder(Variable, Loading),
    y = Loading
  )
) +
  
  # Bars
  geom_col(
    width = 0.70,
    fill = col_loadings
  ) +
  
  # Horizontal bars
  coord_flip() +
  
  # Separate PC1 and PC2
  facet_wrap(
    ~ PC,
    ncol = 2,
    scales = "free_y"
  ) +
  
  # Zero reference
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  
  # Labels
  labs(
    x = "Morphological variable",
    y = "Loading"
  ) +
  
  # Same style as morphospace
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    text = element_text(
      size = 11,
      family = "Arial"
    ),
    
    axis.title = element_text(
      size = 11,
      family = "Arial"
    ),
    
    axis.text = element_text(
      size = 11,
      family = "Arial",
      color = "black"
    ),
    
    strip.text = element_text(
      size = 11,
      family = "Arial"
    ),
    
    strip.background = element_blank(),
    
    panel.spacing = grid::unit(
      1,
      "cm"
    )
  )


# Display figure
print(p_loadings)

############################################################
#EXPORT HIGH-RESOLUTION FIGURE
############################################################

ggsave(
  filename = paste0(
    "F:/Doctorado/TESIS version 3/",
    "capitulo 2.- corregido/figuras/loadings.tiff"
  ),
  plot = p_loadings,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6,
  height = 6,
  units = "in"
)
# ============================================================
# 13. MEAN MORPHOLOGICAL PROPORTIONS ± SD BY SPECIES
# ============================================================

morph_summary_raw <- morph %>%
  group_by(species) %>%
  summarise(
    across(
      all_of(morph_vars),
      list(
        Mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        SD = ~ sd(
          .x,
          na.rm = TRUE
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


print(
  morph_summary_raw,
  width = Inf
)


# ============================================================
# 14. MORPHOLOGICAL SUMMARY — LONG FORMAT
# ============================================================

morph_table_raw <- morph %>%
  pivot_longer(
    cols = all_of(
      morph_vars
    ),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  
  group_by(
    species,
    Variable
  ) %>%
  
  summarise(
    n = sum(
      !is.na(Value)
    ),
    Mean = mean(
      Value,
      na.rm = TRUE
    ),
    SD = sd(
      Value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  mutate(
    `Mean ± SD` = sprintf(
      "%.2f ± %.2f",
      Mean,
      SD
    )
  )


print(
  morph_table_raw,
  n = Inf
)


# ============================================================
# 15. FINAL TABLE
# Variables as rows and species as columns
# ============================================================

morph_table_wide_raw <- morph_table_raw %>%
  
  select(
    species,
    Variable,
    `Mean ± SD`
  ) %>%
  
  pivot_wider(
    names_from = species,
    values_from = `Mean ± SD`
  )


# Preserve morphological variable order
morph_table_wide_raw <- morph_table_wide_raw %>%
  mutate(
    Variable = factor(
      Variable,
      levels = morph_vars
    )
  ) %>%
  arrange(
    Variable
  )


print(
  morph_table_wide_raw,
  n = Inf
)