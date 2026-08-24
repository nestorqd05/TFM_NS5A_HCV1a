############################################################
# 08_NS5A_RAS_combinations.R
# Resumen de combinaciones de RAS por secuencia
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
output_csv  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_combinations.csv"
output_plot <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_combinations.png"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer tabla de RAS
############################################################

ras <- read_csv(ras_file)

############################################################
# 3. Calcular número de RAS por secuencia
############################################################

ras_counts <- ras %>%
  mutate(
    Total_RAS = RAS_28 + RAS_30 + RAS_31 + RAS_32 + RAS_58 + RAS_92 + RAS_93
  ) %>%
  select(Sequence, Total_RAS)

############################################################
# 4. Identificar combinaciones específicas de RAS
############################################################

ras_combo <- ras %>%
  mutate(
    Combo = paste0(
      ifelse(RAS_28, "28-", ""),
      ifelse(RAS_30, "30-", ""),
      ifelse(RAS_31, "31-", ""),
      ifelse(RAS_32, "32-", ""),
      ifelse(RAS_58, "58-", ""),
      ifelse(RAS_92, "92-", ""),
      ifelse(RAS_93, "93-", "")
    )
  ) %>%
  mutate(
    Combo = ifelse(Combo == "", "Sin RAS", Combo)
  ) %>%
  select(Sequence, Combo)

############################################################
# 5. Resumen de combinaciones
############################################################

combo_summary <- ras_combo %>%
  group_by(Combo) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

write_csv(combo_summary, output_csv)

############################################################
# 6. Figura para el TFM
############################################################

p <- ggplot(combo_summary, aes(x = reorder(Combo, -Count), y = Count)) +
  geom_col(fill = "#D55E00") +
  labs(
    title = "Combinaciones de RAS en NS5A (GT1a)",
    x = "Combinación de RAS",
    y = "Número de secuencias"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_plot, p, width = 10, height = 6)

cat("✔ Tabla de combinaciones guardada en:\n", output_csv, "\n")
cat("✔ Figura guardada en:\n", output_plot, "\n")
