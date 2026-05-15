# glmmTMB version (optional for reproducibility)
packageVersion("glmmTMB")

# Essential packages
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(ggplot2)
library(dplyr)
library(car)

# Read dataset
data <- read.csv("biomass.csv")

# View data in RStudio
# View(data)

# Factors and levels
data <- data %>%
  mutate(
    Species   = factor(Species,
                       levels = c("S. acapulcoensis", "S. flavilatus")),
    
    Depth = factor(Depth,
                   levels = c("Shallow", "Intermediate", "Deep")),
    
    transect = factor(transect)
  )

# Quick check
data %>%
  group_by(Species, Depth) %>%
  summarise(
    n = n(),
    sum = sum(biomass),
    .groups = "drop"
  )

# Proportion of zeros
prop_zeros <- mean(data$biomass == 0)

prop_zeros  # ~ 0.37

# General histogram
ggplot(data, aes(x = biomass)) +
  geom_histogram(
    binwidth = 1,
    fill = "gray70",
    color = "black"
  ) +
  theme_classic()

# TWEEDIE MODEL
tweedie_model <- glmmTMB(
  biomass ~ Species * Depth + (1 | transect),
  
  family = tweedie(link = "log"),
  
  data = data
)

summary(tweedie_model)

car::Anova(tweedie_model, type = "II")

# Diagnostics
res <- simulateResiduals(
  tweedie_model,
  n = 1000
)

plot(res)

# Marginal comparisons on response scale
emm <- emmeans(
  tweedie_model,
  ~ Species | Depth,
  type = "response"
)

pairs(emm)

# Predictions for plot
pred <- ggpredict(
  tweedie_model,
  terms = c("Depth", "Species")
)

# Position for dodged points and error bars
pd <- position_dodge(width = 0.5)

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
    y = "Estimated biomass (g/m²)",
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
      expression(italic("Stegastes acapulcoensis")),
      expression(italic("Stegastes flavilatus"))
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 6)
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    legend.text = element_text(
      family = "Arial",
      size = 11
    ),
    
    legend.position = "top",
    
    text = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    
    strip.text = element_text(
      size = 11,
      color = "black"
    )
  )

# Display plot
print(biomass_plot)

# Save figure
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/biomass.tiff",
  
  plot = biomass_plot,
  
  device = "tiff",
  
  dpi = 500,
  
  compression = "lzw",
  
  width = 6,
  height = 6,
  units = "in"
)

