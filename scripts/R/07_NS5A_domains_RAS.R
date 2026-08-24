############################################################
# 07_NS5A_domains_RAS.R
# Anotación de dominios funcionales de NS5A y resumen de RAS
# Usa la tabla NS5A_RAS_table.csv del script 05
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

############################################################
# 1. Rutas
############################################################

ras_file    <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_table.csv"
output_csv  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_by_domain.csv"
output_plot <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_by_domain.png"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer tabla de RAS
############################################################

ras <- read_csv(ras_file)

############################################################
# 3. Definir dominios de NS5A (H77, GT1a)
# Dominio I:   1–100
# Dominio II:  101–200
# Dominio III: 201–447
# Tus posiciones (28, 30, 31, 32, 58, 92, 93) están todas en Dominio I
############################################################

domain_map <- data.frame(
  Position = c(28, 30, 31, 32, 58, 92, 93),
  Domain   = "Dominio I"
)

############################################################
# 4. Convertir tabla de formato ancho a largo
############################################################

ras_long <- ras %>%
  pivot_longer(
    cols = starts_with("Pos_"),
    names_to = "Pos_col",
    values_to = "AA"
  ) %>%
  mutate(
    Position = as.integer(gsub("Pos_", "", Pos_col))
  ) %>%
  left_join(domain_map, by = "Position")

############################################################
# 5. Añadir columnas TRUE/FALSE de RAS en formato largo
############################################################

ras_flags <- ras %>%
  pivot_longer(
    cols = starts_with("RAS_"),
    names_to = "RAS_col",
    values_to = "Is_RAS"
  ) %>%
  mutate(
    Position = as.integer(gsub("RAS_", "", RAS_col))
  ) %>%
  select(Sequence, Position, Is_RAS)

############################################################
# 6. Unir ambas tablas
############################################################

ras_full <- ras_long %>%
  left_join(ras_flags, by = c("Sequence", "Position"))

############################################################
# 7. Resumen por dominio
############################################################

ras_summary <- ras_full %>%
  group_by(Domain, Position) %>%
  summarise(
    Total_RAS = sum(Is_RAS, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(ras_summary, output_csv)

############################################################
# 8. Figura para el TFM
############################################################

p <- ggplot(ras_summary, aes(x = Position, y = Total_RAS)) +
  geom_col(fill = "#0072B2") +
  labs(
    title = "Frecuencia de RAS por posición (Dominio I, NS5A GT1a)",
    x = "Posición (H77)",
    y = "Número de secuencias con RAS"
  ) +
  theme_minimal(base_size = 14)

ggsave(output_plot, p, width = 10, height = 6)

cat("✔ Tabla de RAS por dominio guardada en:\n", output_csv, "\n")
cat("✔ Figura guardada en:\n", output_plot, "\n")
