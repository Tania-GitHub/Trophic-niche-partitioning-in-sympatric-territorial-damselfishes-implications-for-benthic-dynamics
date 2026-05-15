# Load packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)
library(ggrepel)
library(grid)

# --------------------------
# Parameters
# --------------------------
eps <- 0 # Presence threshold for %FO

# --------------------------
# 1) Read and clean data
# --------------------------
df0 <- read_csv(
  "Grouped_items.csv",
  show_col_types = FALSE,
  na = c("", "NA")
)

names(df0) <- str_trim(names(df0))

df0 <- df0 %>% 
  mutate(Species = str_trim(Species)) %>%
  filter(!(is.na(Species) | Species == ""))

# Identify item columns and remove Microplastics if present
cols_items <- setdiff(names(df0), c("Stomach", "Species", "Microplastics"))

# Convert item columns to numeric
df0 <- df0 %>% 
  mutate(across(all_of(cols_items), ~ suppressWarnings(as.numeric(.x))))

# Keep only stomachs with content
df <- df0 %>%
  mutate(total_row = rowSums(across(all_of(cols_items)), na.rm = TRUE)) %>%
  filter(total_row > 0)

# --------------------------
# 2) %FO and %W function
# --------------------------
calc_FO_W_raw <- function(d_species, items, eps = 0) {
  
  n_stomachs <- nrow(d_species)
  
  FO <- colSums(as.matrix(d_species[, items]) > eps, na.rm = TRUE) / 
    n_stomachs * 100
  
  item_sums <- colSums(d_species[, items, drop = FALSE], na.rm = TRUE)
  
  total_sum <- sum(item_sums, na.rm = TRUE)
  
  W <- item_sums / total_sum * 100
  
  tibble(
    Item = names(FO),
    FO_pct = as.numeric(FO),
    W_pct = as.numeric(W),
    N_stomachs = n_stomachs
  )
}

# --------------------------
# 3) Calculate %FO, %W, IRI_mod and %IRI_mod
# --------------------------
res_tab <- df %>%
  group_split(Species) %>%
  map_dfr(~{
    sp <- unique(.x$Species)
    
    calc_FO_W_raw(.x, cols_items, eps = eps) %>%
      mutate(Species = sp, .before = 1)
  }) %>%
  group_by(Species) %>%
  mutate(
    IRI_mod = FO_pct * W_pct,
    IRI_mod_pct = IRI_mod / sum(IRI_mod, na.rm = TRUE) * 100
  ) %>%
  ungroup()

# Check %W
check_W <- res_tab %>%
  group_by(Species) %>%
  summarise(
    W_total = sum(W_pct, na.rm = TRUE),
    .groups = "drop"
  )

print(check_W)

# Save final table
write_csv(
  res_tab %>% arrange(Species, desc(IRI_mod_pct)),
  "FO_W_IRImod_by_species.csv"
)

# --------------------------
# 4) Prepare Costello data
# --------------------------
res_costello <- res_tab %>%
  mutate(
    FO_prop = FO_pct / 100,
    Pi_pct = W_pct,
    Pi_prop = Pi_pct / 100,
    IRI_mod_pct = ifelse(is.na(IRI_mod_pct), 0, IRI_mod_pct),
    Species_lab = dplyr::recode(
      Species,
      "S. acapulcoensis" = "Stegastes acapulcoensis",
      "S. flavilatus" = "Stegastes flavilatus"
    )
  ) %>%
  filter(!is.na(Pi_pct) & FO_pct > 0)

# --------------------------
# 5) Costello plot
# --------------------------
pal_species <- c(
  "S. acapulcoensis" = "#C16540",
  "S. flavilatus" = "#00C1C8"
)

max_iri <- max(res_costello$IRI_mod_pct, na.rm = TRUE)

ann_df <- data.frame(
  Species_lab = c(
    "Stegastes acapulcoensis",
    "Stegastes flavilatus"
  ),
  x = c(0.93, 0.93),
  y = c(98, 98),
  label = c("bold('A')", "bold('B')")
)

p_costello <- ggplot(
  res_costello,
  aes(
    x = FO_prop,
    y = Pi_pct,
    fill = Species,
    size = IRI_mod_pct,
    label = Item
  )
) +
  geom_point(
    shape = 21,
    colour = "grey30",
    stroke = 0.4,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  geom_text_repel(
    family = "Arial",
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.1,
    segment.alpha = 0.4,
    colour = "black",
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 4.5,
    parse = TRUE
  ) +
  facet_wrap(~ Species_lab, nrow = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    labels = number_format(accuracy = 1)
  ) +
  scale_size_continuous(
    name = "%IRI_mod",
    range = c(2.5, 8),
    limits = c(0, max_iri),
    breaks = c(1, 5, 10, 20, 40),
    labels = function(x) paste0(x, "%"),
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill = NA,
        colour = "black",
        stroke = 0.4
      )
    )
  ) +
  scale_fill_manual(
    values = pal_species,
    limits = c("S. acapulcoensis", "S. flavilatus"),
    drop = FALSE,
    guide = "none"
  ) +
  labs(
    x = "Frequency of occurrence (proportion)",
    y = "Percent weight (%W)",
    size = "%IRI_mod"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(size = 11, family = "Arial"),
    axis.text = element_text(size = 10, colour = "black", family = "Arial"),
    axis.title = element_text(size = 11, colour = "black", family = "Arial"),
    strip.text = element_text(size = 11, family = "Arial", face = "italic"),
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10, family = "Arial"),
    panel.spacing.x = unit(0.8, "lines")
  )

print(p_costello)

# Save Costello figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/costello_prop_facets_singlelegend_IRI.tiff",
  plot = p_costello,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 8,
  height = 6,
  units = "in"
)

# --------------------------
# 6) Calculate Bi by species
# --------------------------
Bi_tab <- res_tab %>%
  group_by(Species) %>%
  summarise(
    n_items = n(),
    sum_p2 = sum((IRI_mod_pct / 100)^2, na.rm = TRUE),
    Bi = (1 / (n_items - 1)) * ((1 / sum_p2) - 1),
    .groups = "drop"
  )

print(Bi_tab)

# Save Bi table
write_csv(Bi_tab, "Bi_by_species.csv")

# --------------------------
# 7) Trophic breadth index plot
# --------------------------
p_Bi <- ggplot(Bi_tab, aes(x = Species, y = Bi, fill = Species)) +
  geom_col(
    width = 0.6,
    color = "black"
  ) +
  scale_fill_manual(
    values = pal_species
  ) +
  scale_x_discrete(
    labels = c(
      "S. acapulcoensis" = expression(italic("S. acapulcoensis")),
      "S. flavilatus" = expression(italic("S. flavilatus"))
    )
  ) +
  theme_classic(base_family = "Arial") +
  labs(
    title = "Trophic breadth index (Bi)",
    x = "Species",
    y = expression(B[i])
  ) +
  theme(
    text = element_text(size = 11, family = "Arial"),
    axis.text = element_text(size = 11, color = "black", family = "Arial"),
    axis.title = element_text(size = 11, color = "black", family = "Arial"),
    legend.position = "none"
  )

print(p_Bi)

# Save Bi figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/trophic_breadth_Bi.tiff",
  plot = p_Bi,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 5,
  height = 5,
  units = "in"
)