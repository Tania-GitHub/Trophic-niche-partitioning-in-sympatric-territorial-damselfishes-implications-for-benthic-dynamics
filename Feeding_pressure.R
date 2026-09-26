###############################################
# FEEDING PRESSURE
# FP = BR × (B / A)
#
# BR = bite rate (bites h-1)
# B  = species biomass per transect (g)
# A  = transect area (80 m2)
###############################################

# Load packages
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(car)
library(ggplot2)
library(dplyr)

###############################################
# Load data
###############################################

data <- read.csv("PA.csv")

# Check data
head(data)
names(data)
str(data)

###############################################
# Factors
###############################################

data$Species <- factor(data$Species)

data$Depth <- factor(
  data$Depth,
  levels = c("Shallow", "Intermediate", "Deep")
)

# Create a unique ID for each of the 36 transects
# because T1-T12 are repeated among reef zones

data$Transect_ID <- interaction(
  data$Depth,
  data$Transect,
  drop = TRUE
)

# Check number of unique transects
nlevels(data$Transect_ID)
# Should return 36


###############################################
# Calculate feeding pressure
###############################################

# Surveyed area per transect (m2)
A <- 80

data <- data %>%
  mutate(
    Biomass_density = Biomass / A,
    Feeding_pressure = Bite_rate * Biomass_density
  )

# Equivalent:
# Feeding_pressure = Bite_rate * (Biomass / 80)


###############################################
# Inspect calculations
###############################################

head(
  data[, c(
    "Transect",
    "Transect_ID",
    "Depth",
    "Species",
    "Bite_rate",
    "Biomass",
    "Biomass_density",
    "Feeding_pressure"
  )]
)

summary(data$Feeding_pressure)

# Percentage of zeros
mean(data$Feeding_pressure == 0) * 100


###############################################
# GLMM
###############################################

# Feeding pressure is continuous, non-negative,
# and includes true zeros.
# Tweedie distribution with log link.

fp_model <- glmmTMB(
  Feeding_pressure ~ Species * Depth +
    (1 | Transect_ID),
  family = tweedie(link = "log"),
  data = data
)

summary(fp_model)


###############################################
# Type II tests
###############################################

Anova(
  fp_model,
  type = 2
)


###############################################
# Residual diagnostics
###############################################

sim <- simulateResiduals(
  fittedModel = fp_model,
  n = 1000
)

plot(sim)

testDispersion(sim)

testZeroInflation(sim)


###############################################
# Estimated marginal means
###############################################

emm <- emmeans(
  fp_model,
  ~ Species | Depth,
  type = "response"
)

emm


###############################################
# Pairwise comparisons between species
# within each reef zone
###############################################

pairs(
  emm,
  adjust = "tukey"
)


###############################################
# Predictions for figure
###############################################

pred <- emmeans(
  fp_model,
  ~ Species * Depth,
  type = "response"
)

df_pred <- as.data.frame(pred)

df_pred


###############################################
# Publication-ready plot
###############################################

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
  
  # 95% confidence intervals
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    position = position_dodge(width = 0.25),
    width = 0.2,
    linewidth = 0.6
  ) +
  
  # Estimated means
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3.5,
    shape = 21,
    stroke = 0.5
  ) +
  
  labs(
    y = expression(
      "Feeding pressure index (bites h"^{-1} ~
        "\u00D7 g m"^{-2} * ")"
    ),
    x = "Reef zone",
    color = NULL,
    fill = NULL
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
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.text = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
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


###############################################
# Save figure
###############################################

ggsave(
  filename =
    "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/feeding_pressure.tiff",
  
  plot = feeding_pressure_plot,
  
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  
  width = 6,
  height = 6,
  units = "in"
)
