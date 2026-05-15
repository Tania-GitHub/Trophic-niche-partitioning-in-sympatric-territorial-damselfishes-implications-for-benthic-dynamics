# Packages
library(tidyverse)
library(ggpubr)

# 1) Read data
data <- read.csv("Morphological metrics.csv", check.names = FALSE)

# Suggested order for the legend
data$species <- factor(
  data$species,
  levels = c("S. acapulcoensis", "S. flavilatus")
)

# 2) Morphological variables (already numeric)
variables <- c("2", "3", "4", "5", "6", "7")
stopifnot(all(variables %in% names(data)))

# 3) Shapiro–Wilk test by variable and species
cat("=== Shapiro–Wilk test by variable and species ===\n")

for (var in variables) {
  
  cat("\n---", var, "---\n")
  
  for (sp in levels(data$species)) {
    
    x <- data %>%
      filter(species == sp) %>%
      pull(all_of(var))
    
    x <- x[is.finite(x)]
    
    if (length(x) < 3L) {
      
      cat("Species:", sp, "- n <", 3, "(Shapiro test not applicable)\n")
      
    } else {
      
      sw <- tryCatch(shapiro.test(x), error = function(e) NULL)
      
      if (is.null(sw)) {
        
        cat("Species:", sp, "- error in Shapiro test\n")
        
      } else {
        
        cat(
          "Species:", sp,
          "- W =", round(unname(sw$statistic), 3),
          "- p =", signif(sw$p.value, 3), "\n"
        )
      }
    }
  }
}

# 4) Histograms and QQ-plots
for (var in variables) {
  
  p_hist <- ggplot(data, aes(x = .data[[var]] / 10)) +  # mm → cm
    geom_histogram(
      bins = 15,
      fill = "#4C72B0",
      color = "white"
    ) +
    facet_wrap(~ species, nrow = 1) +
    labs(
      title = paste("Histogram of metric", var),
      x = paste("Metric", var, "(cm)"),
      y = "Frequency"
    ) +
    theme_minimal(base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_hist)
  
  p_qq <- ggplot(data, aes(sample = .data[[var]] / 10)) +  # mm → cm
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ species, nrow = 1) +
    labs(
      title = paste("QQ plot of metric", var),
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    theme_minimal(base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p_qq)
}

# 5) Long-format data + conversion to cm
data_long <- data %>%
  pivot_longer(
    cols = all_of(variables),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(
    Value = Value / 10,  # mm → cm
    Variable = factor(
      Variable,
      levels = c("2", "3", "4", "5", "6", "7")
    )
  )

# 6) Wilcoxon test by variable
wilcox_df <- data_long %>%
  group_by(Variable) %>%
  summarise(
    p = tryCatch(
      wilcox.test(Value ~ species)$p.value,
      error = function(e) NA_real_
    ),
    
    W = tryCatch(
      unname(wilcox.test(Value ~ species)$statistic),
      error = function(e) NA_real_
    ),
    
    y.position = max(Value, na.rm = TRUE) * 1.05,
    
    group1 = "S. acapulcoensis",
    group2 = "S. flavilatus",
    
    .groups = "drop"
  )
# Wilcoxon results
wilcox_df %>%
  mutate(p = signif(p, 3)) %>%
  select(Variable, W, p) %>%
  print(n = Inf)

# 7) Final plot
pal_species <- c(
  "S. acapulcoensis" = "#C16540",
  "S. flavilatus" = "#00C1C8"
)

labs_species <- c(
  "S. acapulcoensis" = "Stegastes acapulcoensis",
  "S. flavilatus" = "Stegastes flavilatus"
)

p <- ggplot(
  data_long,
  aes(x = Variable, y = Value, fill = species)
) +
  
  geom_boxplot(
    alpha = 0.8,
    color = "black",
    width = 0.7,
    position = position_dodge(width = 0.8),
    outlier.shape = 21
  ) +
  
  stat_summary(
    aes(group = species),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  
  scale_fill_manual(
    values = pal_species,
    breaks = levels(data$species),
    labels = labs_species
  ) +
  
  labs(
    y = "Length (cm)",
    x = "Morphological metric"
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    text = element_text(
      family = "Arial",
      size = 11,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      face = "italic",
      size = 11,
      color = "black"
    ),
    
    legend.position = "top",
    legend.margin = margin(t = 5),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# Save figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/metrics.tiff",
  plot = p,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6,
  height = 6,
  units = "in"
)



