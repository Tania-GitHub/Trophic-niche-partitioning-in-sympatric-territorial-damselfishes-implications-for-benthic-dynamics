# ============================================================
# DIETARY ITEM COMPOSITION
# PERMANOVA + PERMDISP + PCoA + SIMPER
# ============================================================

library(vegan)
library(ape)
library(ggplot2)
library(dplyr)


# ============================================================
# 1. IMPORT DATA
# ============================================================

datos <- read.csv(
  "ITEM.csv",
  header = TRUE
)

# Metadata
meta <- datos[, 1:2]

# Dietary items
mat <- datos[, -c(1, 2)]


# ============================================================
# 2. DATA PREPARATION
# ============================================================

# Remove microplastics
mat <- mat[
  ,
  colnames(mat) != "Microplasticos"
]

# Convert dietary variables to numeric
mat <- as.data.frame(
  lapply(
    mat,
    as.numeric
  )
)

# Replace NA with 0
mat[is.na(mat)] <- 0


# ============================================================
# 3. PRESENCE-ABSENCE MATRIX
# ============================================================

mat_pa <- decostand(
  mat,
  method = "pa"
)

# Check
head(mat_pa)


# ============================================================
# 4. JACCARD DISSIMILARITY
# ============================================================

dist_jac <- vegdist(
  mat_pa,
  method = "jaccard"
)


# ============================================================
# 5. PERMANOVA
# ============================================================

set.seed(123)

permanova <- adonis2(
  dist_jac ~ Species,
  data = meta,
  permutations = 9999
)

print(permanova)


# ============================================================
# 7. PCoA
# ============================================================

pcoa_res <- pcoa(
  dist_jac
)

# Coordinates of first two axes
pcoa_axes <- as.data.frame(
  pcoa_res$vectors[, 1:2]
)

colnames(pcoa_axes) <- c(
  "PCoA1",
  "PCoA2"
)

# Add metadata
pcoa_df <- cbind(
  meta,
  pcoa_axes
)

# Relative eigenvalues (%)
eig <- pcoa_res$values$Relative_eig * 100

# Check variance represented
eig[1:5]


# ============================================================
# 8. CONVEX HULLS
# ============================================================

pcoa_hulls <- pcoa_df %>%
  group_by(Species) %>%
  slice(
    chull(PCoA1, PCoA2)
  ) %>%
  ungroup()


# ============================================================
# 9. PCoA PLOT
# ============================================================

composicion <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Species
  )
) +
  
  # Convex hulls
  geom_polygon(
    data = pcoa_hulls,
    aes(
      x = PCoA1,
      y = PCoA2,
      group = Species,
      fill = Species,
      color = Species
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
  
  # Species colors
  scale_color_manual(
    values = c(
      "S. acapulcoensis" = "#C16540",
      "S. flavilatus"    = "#00C1C8"
    ),
    labels = c(
      "Stegastes acapulcoensis",
      "Stegastes flavilatus"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "S. acapulcoensis" = "#C16540",
      "S. flavilatus"    = "#00C1C8"
    ),
    guide = "none"
  ) +
  
  # Axis labels
  labs(
    x = paste0(
      "PCoA1 (",
      round(eig[1], 2),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(eig[2], 2),
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
      family = "Arial",
      size = 11
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 11
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      family = "Arial",
      size = 11,
      face = "italic"
    ),
    
    legend.position = "top"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )


print(composicion)


# ============================================================
# 10. SAVE PCoA FIGURE
# ============================================================

ggsave(
  filename = paste0(
    "F:/Doctorado/TESIS version 3/",
    "capitulo 2.- corregido/figuras/",
    "composicion.tiff"
  ),
  plot = composicion,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6,
  height = 6,
  units = "in"
)


# ============================================================
# 11. SIMPER
# Presence-absence dietary data
# NOTE: vegan::simper() uses Bray-Curtis dissimilarity
# ============================================================

set.seed(123)

simper_res <- simper(
  mat_pa,
  meta$Species,
  permutations = 9999
)

# Full SIMPER summary
sim <- summary(
  simper_res
)

sim


# ============================================================
# 12. EXTRACT SPECIES COMPARISON
# ============================================================

# There are only two species, therefore one comparison
sim_df <- sim[[1]]

# Convert row names (dietary items) to a column
sim_df$Item <- rownames(
  sim_df
)

# Sort dietary items from highest to lowest contribution
sim_df <- sim_df %>%
  arrange(
    desc(average)
  )

# Calculate cumulative contribution
sim_df <- sim_df %>%
  mutate(
    cumulative = cumsum(average) /
      sum(average)
  )


# ============================================================
# 13. ITEMS CONTRIBUTING TO 70% OF DISSIMILARITY
# ============================================================

sim_70 <- sim_df %>%
  filter(
    cumulative <= 0.70
  )

sim_70


# ============================================================
# 14. OVERALL AVERAGE DISSIMILARITY
# ============================================================

simper_res[[1]]$overall

