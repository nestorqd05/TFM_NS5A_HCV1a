############################################################
# 09_NS5A_entropy.R
# Cálculo de variabilidad por posición usando entropía de Shannon
# Usa el alineamiento NS5A_1a_aligned.fasta del script 04
############################################################

library(Biostrings)
library(ggplot2)
library(dplyr)

############################################################
# 1. Rutas
############################################################

aligned_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/aligned/NS5A_1a_aligned.fasta"

output_csv  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_entropy.csv"
output_plot <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_entropy.png"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer alineamiento
############################################################

aln <- readAAStringSet(aligned_file)
aln_matrix <- as.matrix(aln)

n_seq <- nrow(aln_matrix)
n_pos <- ncol(aln_matrix)

############################################################
# 3. Función de entropía de Shannon
############################################################

shannon_entropy <- function(column) {
  column <- column[column != "-"]  # eliminar gaps
  freqs <- table(column) / length(column)
  -sum(freqs * log2(freqs))
}

############################################################
# 4. Calcular entropía por posición
############################################################

entropy_values <- sapply(1:n_pos, function(i) shannon_entropy(aln_matrix[, i]))

entropy_df <- data.frame(
  Position = 1:n_pos,
  Entropy = entropy_values
)

write.csv(entropy_df, output_csv, row.names = FALSE)

############################################################
# 5. Figura para el TFM
############################################################

p <- ggplot(entropy_df, aes(x = Position, y = Entropy)) +
  geom_line(color = "#009E73", size = 1) +
  labs(
    title = "Variabilidad por posición en NS5A (Entropía de Shannon)",
    x = "Posición (alineamiento H77)",
    y = "Entropía"
  ) +
  theme_minimal(base_size = 14)

ggsave(output_plot, p, width = 12, height = 6)

cat("✔ Entropía por posición guardada en:\n", output_csv, "\n")
cat("✔ Figura guardada en:\n", output_plot, "\n")
