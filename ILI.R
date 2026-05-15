# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(car)

# Load data
data <- read.csv(
  "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/para someter/bases de datos/LI.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

# Calculate the Intestinal Length Index (ILI)
data <- data %>%
  mutate(ILI = IL / SL)

# Evaluate normality
# Histogram
ggplot(data, aes(x = ILI, fill = species)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity", color = "black") +
  theme_classic() +
  facet_wrap(~ species) +
  labs(title = "ILI distribution by species")

# Shapiro-Wilk test by species
cat("Shapiro-Wilk test by species:\n")
by(data$ILI, data$species, shapiro.test)

# Evaluate homogeneity of variances
cat("\nLevene's test:\n")
leveneTest(ILI ~ species, data = data)

# Statistical test
# Student's t-test assuming equal variances
t.test(ILI ~ species, data = data, var.equal = TRUE)

# Visualization with boxplot
# Ensure species is a factor with the correct order
data$species <- factor(
  data$species,
  levels = c("S. acapulcoensis", "S. flavilatus")
)

# Plot
ILI <- ggplot(data, aes(x = species, y = ILI, fill = species)) +
  geom_boxplot(alpha = 0.6, color = "black") +
  
  #stat_compare_means(
  #  method = "t.test",
  #  label = "p",
  #  method.args = list(var.equal = TRUE)
  #) +  # Display p-value without asterisks
  
  scale_fill_manual(values = c(
    "S. acapulcoensis" = "#C16540",
    "S. flavilatus" = "#00C1C8"
  )) +
  
  scale_x_discrete(labels = c(
    "S. acapulcoensis" = expression(italic("S. acapulcoensis")),
    "S. flavilatus" = expression(italic("S. flavilatus"))
  )) +
  
  labs(
    x = NULL,
    y = "Intestinal Length Index (ILI)"
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    text = element_text(size = 12, color = "black", family = "Arial"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

print(ILI)
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 3.- dieta/figures/ILI.tiff",
  plot = ILI,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 18, height = 14, units = "cm"
)
