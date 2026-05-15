# Load packages
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(car)
library(ggplot2)

# Load data
data <- read.csv("Feeding_pressure.csv")

# Ensure variables are factors
data$Species <- factor(data$Species)

data$Depth <- factor(
  data$Depth,
  levels = c("Shallow", "Intermediate", "Deep")
)

# Zero-inflated negative binomial model (ZINB)
zinb_model <- glmmTMB(
  Feeding_pressure ~ Species * Depth,
  
  ziformula = ~1,
  
  family = nbinom2,
  
  data = data
)

summary(zinb_model)

Anova(zinb_model)

# Residual diagnostics
sim <- simulateResiduals(zinb_model)

plot(sim)

# Pairwise comparisons
emmeans(
  zinb_model,
  pairwise ~ Species | Depth,
  type = "response"
)

# Marginal predictions
pred <- emmeans(
  zinb_model,
  ~ Species * Depth,
  type = "response"
)

df_pred <- as.data.frame(pred)

# Publication-ready plot
feeding_pressure_plot <- ggplot(
  df_pred,
  aes(
    x = Depth,
    y = response,
    group = Species,
    color = Species,
    fill = Species
  )
) +
  
  # Error bars (95% CI)
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    
    position = position_dodge(width = 0.25),
    
    width = 0.2,
    
    linewidth = 0.6
  ) +
  
  # Points filled by species
  geom_point(
    position = position_dodge(width = 0.25),
    
    size = 3.5,
    
    shape = 21,
    
    stroke = 0.5
  ) +
  
  labs(
    y = "Feeding pressure (expected mean ± 95% CI)",
    x = "Reef zone",
    color = NULL,
    fill = NULL
  ) +
  
  # Custom colors and fills
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
  
  scale_fill_manual(
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
  
  # Lower visual margin
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    
    legend.position = "top",
    
    legend.direction = "horizontal",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      family = "Arial",
      size = 10
    ),
    
    axis.line = element_line(
      linewidth = 0.4
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4
    ),
    
    plot.margin = margin(10, 10, 10, 10)
  )

print(feeding_pressure_plot)

# Save figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/feeding_pressure.tiff",
  
  plot = feeding_pressure_plot,
  
  device = "tiff",
  
  dpi = 500,
  
  compression = "lzw",
  
  width = 6,
  height = 6,
  units = "in"
)
