library(vegan)
library(ape)
library(ggplot2)

datos <- read.csv("Compo.csv", header = TRUE)

meta <- datos[, 1:2]
mat  <- datos[, -c(1,2)]

# Quitar Microplasticos
mat <- mat[, colnames(mat) != "Microplasticos"]

# Convertir a numérico
mat <- as.data.frame(lapply(mat, as.numeric))

# Reemplazar NA por 0
mat[is.na(mat)] <- 0

# Distancia Jaccard
dist_jac <- vegdist(mat, method = "jaccard", binary = TRUE)

# PERMANOVA
permanova <- adonis2(dist_jac ~ Especie, data = meta, permutations = 9999)
print(permanova)

# Homogeneidad
#disp <- betadisper(dist_jac, meta$Especie)
#anova(disp)
#permutest(disp)

# =========================
# PCoA
# =========================
pcoa_res <- pcoa(dist_jac)

pcoa_axes <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_axes) <- c("PCoA1", "PCoA2")

pcoa_df <- cbind(meta, pcoa_axes)

eig <- pcoa_res$values$Relative_eig * 100

# =========================
# GRAFICO
# =========================

composicion <-ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Especie, fill = Especie)) +
  stat_ellipse(type = "t", linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.9)+
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("S. acapulcoensis" = "#C16540",
               "S. flavilatus" = "#00C1C8"),
    labels = c("Stegastes acapulcoensis",
               "Stegastes flavilatus")
  ) +
  scale_fill_manual(
    values = c("S. acapulcoensis" = "#C16540",
               "S. flavilatus" = "#00C1C8"),
    labels = c("Stegastes acapulcoensis",
               "Stegastes flavilatus")
  ) +
  labs(x = paste0("PCoA1 (", round(eig[1], 2), "%)"),
       y = paste0("PCoA2 (", round(eig[2], 2), "%)"),
       color = NULL,
       fill = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 11),
    axis.title = element_text(family = "Arial", size = 11),
    axis.text = element_text(family = "Arial", size = 11, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial", size = 11, face = "italic"),
    legend.position = "top"
  )
print(composicion)

# Guardar figura
ggsave(
  filename = "F:/Doctorado/TESIS version 3/capitulo 2.- corregido/figuras/composicion.tiff",
  plot = composicion,
  device = "tiff",
  dpi = 500,
  compression = "lzw",
  width = 6, height = 6, units = "in"
)

# =========================
# SIMPER
# =========================

library(vegan)

# SIMPER por especie
simper_res <- simper(mat, meta$Especie, permutations = 9999)

# Ver resultados
summary(simper_res)

sim <- summary(simper_res)

# Extraer comparación
sim_df <- sim[[1]]

# Filtrar hasta 70%
sim_df$acum <- cumsum(sim_df$average) / sum(sim_df$average)

sim_70 <- sim_df[sim_df$acum <= 0.7, ]

sim_70

#Promedio de disimilitud
simper_res[[1]]$overall


