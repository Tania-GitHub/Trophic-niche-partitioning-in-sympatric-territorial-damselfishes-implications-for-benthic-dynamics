########################################################
# Biomass analysis
# Stegastes acapulcoensis vs. S. flavilatus
########################################################

# glmmTMB version (optional for reproducibility)
packageVersion("glmmTMB")

# Packages
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(ggplot2)
library(dplyr)
library(car)


########################################################
# 1. Read and prepare data
########################################################

data <- read.csv("biomass.csv")

data <- data %>%
  mutate(
    Species = factor(
      Species,
      levels = c(
        "S. acapulcoensis",
        "S. flavilatus"
      )
    ),
    
    Depth = factor(
      Depth,
      levels = c(
        "Shallow",
        "Intermediate",
        "Deep"
      )
    ),
    
    transect = factor(transect)
  )


########################################################
# 2. Preliminary checks
########################################################

# Sample size and total biomass by species and depth
data %>%
  group_by(Species, Depth) %>%
  summarise(
    n = n(),
    sum_biomass = sum(biomass),
    .groups = "drop"
  )

# Proportion of zeros
prop_zeros <- mean(data$biomass == 0)
prop_zeros

# Biomass distribution
ggplot(data, aes(x = biomass)) +
  geom_histogram(
    binwidth = 1,
    fill = "gray70",
    color = "black"
  ) +
  theme_classic()


########################################################
# 3. Tweedie GLMM
########################################################

tweedie_model <- glmmTMB(
  biomass ~ Species * Depth + (1 | transect),
  
  family = tweedie(
    link = "log"
  ),
  
  data = data
)

summary(tweedie_model)

# Type II Wald chi-square tests
car::Anova(
  tweedie_model,
  type = "II"
)


########################################################
# 4. Model diagnostics
########################################################

res <- simulateResiduals(
  tweedie_model,
  n = 1000
)

plot(res)


########################################################
# 5. Estimated marginal means and species contrasts
########################################################

emm <- emmeans(
  tweedie_model,
  ~ Species | Depth,
  type = "response"
)

# Pairwise species comparisons within each depth zone
pairs(emm)


########################################################
# 6. Biomass ratios and 95% confidence intervals
########################################################

ratio_emm <- contrast(
  emm,
  method = "pairwise"
)

ratio_results <- summary(
  ratio_emm,
  infer = c(TRUE, TRUE),
  type = "response"
)

ratio_results

# Convert to data frame for plotting
ratio_df <- as.data.frame(
  ratio_results
)

# Preserve ecological order of reef zones
ratio_df$Depth <- factor(
  ratio_df$Depth,
  levels = c(
    "Shallow",
    "Intermediate",
    "Deep"
  )
)


########################################################
# 7. Estimated biomass plot
########################################################

pred <- ggpredict(
  tweedie_model,
  terms = c(
    "Depth",
    "Species"
  )
)

pd <- position_dodge(
  width = 0.5
)

biomass_plot <- ggplot(
  pred,
  aes(
    x = x,
    y = predicted,
    color = group
  )
) +
  
  geom_point(
    position = pd,
    size = 3
  ) +
  
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    position = pd,
    width = 0.2,
    linewidth = 0.8
  ) +
  
  labs(
    x = "Reef zone",
    y = expression(
      "Estimated biomass (g " * m^{-2} * ")"
    ),
    color = NULL
  ) +
  
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
      expression(
        italic("Stegastes acapulcoensis")
      ),
      expression(
        italic("Stegastes flavilatus")
      )
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 6)
  ) +
  
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    legend.position = "top",
    
    legend.text = element_text(
      family = "Arial",
      size = 11
    ),
    
    text = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 11,
      color = "black"
    )
  )

print(biomass_plot)


########################################################
# 8. Biomass ratio plot
# Magnitude of the difference between species
########################################################

ratio_plot <- ggplot(
  ratio_df,
  aes(
    x = Depth,
    y = ratio,
    group = 1
  )
) +
  
  # Trend across reef zones
  geom_line(
    linewidth = 0.8,
    color = "#006D61"
  ) +
  
  # Estimated biomass ratios
  geom_point(
    size = 3,
    color = "#006D61"
  ) +
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.15,
    linewidth = 0.8,
    color = "#006D61"
  ) +
  
  # Ratio = 1 indicates equal biomass
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  
  # Logarithmic scale because ratios span a wide range
  scale_y_log10() +
  
  # Axis labels
  labs(
    x = "Reef zone",
    y = expression(
      atop(
        "Biomass ratio",
        "(" * italic("S. acapulcoensis") /
          italic("S. flavilatus") * ")"
      )
    )
  ) +
  
  # Theme
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    text = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 11,
      color = "black"
    )
  )

# Display plot
print(ratio_plot)

########################################################
# 9. Combine biomass and biomass ratio plots
########################################################

# Package to combine plots
library(patchwork)

# Add panel labels in the upper-right corner
biomass_plot_A <- biomass_plot +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "A",
    hjust = 1.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold",
    family = "Arial"
  )

ratio_plot_B <- ratio_plot +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "B",
    hjust = 1.3,
    vjust = 1.3,
    size = 5,
    fontface = "bold",
    family = "Arial"
  )

# Combine panels horizontally
biomass_combined <- biomass_plot_A + ratio_plot_B +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

# Display combined figure
print(biomass_combined)

########################################################
# 10. Save combined figure
########################################################

ggsave(
  filename = paste0(
    "F:/Doctorado/TESIS version 3/",
    "capitulo 2.- corregido/figuras/biomass_combined.tiff"
  ),
  
  plot = biomass_combined,
  
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  
  width = 8,
  height = 4,
  units = "in"
)
