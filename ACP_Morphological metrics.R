## =========================================
## PCA (raw metrics) + biplot
## (without ellipse legend, axes from -5 to 5.5)
## Morphological metrics.csv
## =========================================

# Packages
pkgs <- c("tidyverse", "factoextra", "ggrepel", "showtext", "ragg")

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(pkgs, library, character.only = TRUE))

# 1) Data
df <- read.csv(
  "Morphological metrics.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

vars <- c("2", "3", "4", "5", "6", "7")

stopifnot(all(c("individual", "species", vars) %in% names(df)))

df <- df |>
  mutate(
    species = factor(
      species,
      levels = c("S. acapulcoensis", "S. flavilatus")
    )
  ) |>
  drop_na(all_of(vars))

# 2) PCA
pca <- prcomp(
  df[, vars],
  center = TRUE,
  scale. = TRUE
)

# Individual scores
scores_ind <- as.data.frame(pca$x) |>
  mutate(
    individual = df$individual,
    species = df$species
  )

print(scores_ind)

# Explained variance
ve <- (pca$sdev^2) / sum(pca$sdev^2)

lab_pc1 <- paste0("PC1 (", round(100 * ve[1], 1), "%)")
lab_pc2 <- paste0("PC2 (", round(100 * ve[2], 1), "%)")

print(pca)

# Variable loadings
round(pca$rotation[, 1:2], 3)

# 3) Colors and labels
col_species <- c(
  "S. acapulcoensis" = "#C16540",
  "S. flavilatus"    = "#00C1C8"
)

long_labels <- c(
  expression(italic("Stegastes acapulcoensis")),
  expression(italic("Stegastes flavilatus"))
)

# 4) Coordinates for ellipses
ind <- get_pca_ind(pca)

ind_df <- as.data.frame(ind$coord) |>
  mutate(species = df$species)

# 5) Biplot
raw_biplot <- fviz_pca_biplot(
  pca,
  geom        = "point",
  habillage   = df$species,
  addEllipses = FALSE,
  label       = "var",
  repel       = TRUE,
  col.var     = "black",
  pointshape  = 16,
  pointsize   = 3,
  palette     = col_species
) +
  
  stat_ellipse(
    data = ind_df,
    aes(
      x = Dim.1,
      y = Dim.2,
      color = species,
      fill = species
    ),
    type = "t",
    level = 0.95,
    linewidth = 0.8,
    alpha = 0.18,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = col_species,
    labels = long_labels
  ) +
  
  scale_fill_manual(
    values = col_species,
    guide = "none"
  ) +
  
  labs(
    title = NULL,
    x = lab_pc1,
    y = lab_pc2
  ) +
  
  xlim(-5, 5.5) +
  ylim(-5, 5.5) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    text = element_text(size = 11, family = "Arial"),
    
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
  ) +
  
  guides(
    shape = "none",
    color = guide_legend(
      override.aes = list(size = 3)
    )
  )

print(raw_biplot)



# 6) Export high-resolution figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/PCA.tiff",
  plot = raw_biplot,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6,
  height = 6,
  units = "in"
)
