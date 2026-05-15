# Load required packages
library(glmmTMB)
library(DHARMa)
library(tidyr)
library(dplyr)
library(ggplot2)
library(car)
library(emmeans)

# Read data
data <- read.csv("Bite_rate.csv", check.names = FALSE)

# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = c(Algal_turf, Coral, Sponge, Sediment, Water_column),
    names_to = "item",
    values_to = "frequency"
  ) %>%
  mutate(
    Species = factor(
      Species,
      levels = c("S. acapulcoensis", "S. flavilatus")
    ),
    Depth = factor(
      Depth,
      levels = c("Shallow", "Intermediate", "Deep")
    ),
    item = factor(
      item,
      levels = c("Algal_turf", "Coral", "Sponge", "Sediment", "Water_column"),
      labels = c("AT", "C", "Sp", "S", "WC")
    ),
    Territory = factor(Territory)
  )

# --------------------------
# Zero check
# --------------------------
cat("\n--- Global zeros ---\n")

data_long %>%
  summarise(
    zeros = sum(frequency == 0),
    n = n(),
    prop_zeros = zeros / n
  ) %>%
  print()

# --------------------------
# Fit model without triple interaction
# --------------------------
model_no_triple <- glmmTMB(
  frequency ~ Species * Depth + Species * item + Depth * item + (1 | Territory),
  data = data_long,
  family = tweedie(link = "log")
)

# Model summary and Type III test
summary(model_no_triple)

car::Anova(model_no_triple, type = 3)

# --------------------------
# Residual diagnostics
# --------------------------
residuals_model <- simulateResiduals(model_no_triple)

plot(residuals_model)

# --------------------------
# Post-hoc comparisons
# --------------------------

# Species comparisons within each item and depth
comparisons_item_depth <- emmeans(
  model_no_triple,
  pairwise ~ Species | item * Depth,
  type = "response"
)

comparisons_item_depth

# Species comparisons within each item
comparisons_item <- emmeans(
  model_no_triple,
  pairwise ~ Species | item,
  type = "response"
)

comparisons_item

# Species comparisons within each depth
comparisons_depth <- emmeans(
  model_no_triple,
  pairwise ~ Species | Depth,
  type = "response"
)

comparisons_depth

# --------------------------
# Model-based marginal estimates
# --------------------------
em_full <- emmeans(
  model_no_triple,
  ~ Species * Depth * item
)

em_df <- summary(
  em_full,
  type = "response"
)

# Ensure factor order
em_df <- em_df %>%
  mutate(
    Species = factor(
      Species,
      levels = c("S. acapulcoensis", "S. flavilatus")
    ),
    Depth = factor(
      Depth,
      levels = c("Shallow", "Intermediate", "Deep")
    ),
    item = factor(
      item,
      levels = c("AT", "C", "Sp", "S", "WC")
    )
  )

# --------------------------
# Point-range plot
# --------------------------
bite_rate_plot <- ggplot(
  em_df,
  aes(
    x = item,
    y = response,
    color = Species
  )
) +
  geom_pointrange(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    position = position_dodge(width = 0.5),
    linewidth = 0.9
  ) +
  facet_wrap(~ Depth) +
  scale_color_manual(
    values = c(
      "S. acapulcoensis" = "#C16540",
      "S. flavilatus" = "#00C1C8"
    ),
    breaks = c(
      "S. acapulcoensis",
      "S. flavilatus"
    ),
    labels = c(
      expression(italic(Stegastes~acapulcoensis)),
      expression(italic(Stegastes~flavilatus))
    )
  ) +
  labs(
    x = "Foraging site",
    y = "Estimated bite rate (± CI)",
    color = NULL,
    shape = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    axis.text = element_text(
      family = "Arial",
      size = 10,
      colour = "black"
    ),
    strip.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    legend.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    legend.position = "top"
  )

print(bite_rate_plot)

# Save point-range plot
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/bite_rate.tiff",
  plot = bite_rate_plot,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 8,
  height = 6,
  units = "in"
)

# --------------------------
# Heatmap
# --------------------------
heatmap_bites <- ggplot(
  em_df,
  aes(
    x = item,
    y = Depth,
    fill = response
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = round(response, 2)),
    family = "Arial",
    size = 3.2,
    color = "black"
  ) +
  facet_wrap(
    ~ Species,
    labeller = labeller(
      Species = c(
        "S. acapulcoensis" = "Stegastes acapulcoensis",
        "S. flavilatus" = "Stegastes flavilatus"
      )
    )
  ) +
  scale_fill_gradient2(
    name = "Estimated\nbite rate",
    low = "white",
    mid = "#D8F4EF",
    high = "#006D61",
    midpoint = median(em_df$response)
  ) +
  labs(
    x = "Foraging site",
    y = "Depth"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    axis.text.x = element_text(
      family = "Arial",
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 10.2,
      colour = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "italic",
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    legend.title = element_text(
      family = "Arial",
      size = 10,
      colour = "black"
    ),
    legend.text = element_text(
      family = "Arial",
      size = 10,
      colour = "black"
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(heatmap_bites)

# Save heatmap
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/heatmap_bite_rate.tiff",
  plot = heatmap_bites,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 8,
  height = 6,
  units = "in"
)
