############################################################
# 05_detect_RAS_NS5A.R
# Detección de RAS en NS5A usando la numeración de H77
# Funciona sobre el alineamiento generado con MUSCLE (script 04)
############################################################

library(Biostrings)
library(stringr)

############################################################
# 1. Rutas
############################################################

aligned_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/aligned/NS5A_1a_aligned.fasta"
output_file  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_RAS_table.csv"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer alineamiento
############################################################

cat("Leyendo alineamiento...\n")
aln <- readAAStringSet(aligned_file)

if (!"H77_reference" %in% names(aln)) {
  stop("❌ ERROR: H77_reference no está en el alineamiento. Revisa el script 04.")
}

cat("✔ Alineamiento cargado.\n")

############################################################
# 3. Extraer la secuencia H77_reference
############################################################

h77 <- aln["H77_reference"]

# Convertir a vector de caracteres
h77_vec <- strsplit(as.character(h77), "")[[1]]

############################################################
# 4. Posiciones de interés (RAS NS5A genotipo 1a)
# Numeración oficial basada en H77
############################################################

ras_positions <- c(28, 30, 31, 32, 58, 92, 93)

# RAS conocidas en NS5A GT1a
ras_catalog <- list(
  "28" = c("M", "V"),
  "30" = c("K", "R"),
  "31" = c("M", "V", "F"),
  "32" = c("L", "M"),
  "58" = c("D", "S"),
  "92" = c("E", "G"),
  "93" = c("H", "N", "C")
)

############################################################
# 5. Función para obtener el aminoácido en la posición H77
#    teniendo en cuenta los gaps del alineamiento
############################################################

get_residue <- function(seq_vec, h77_vec, pos) {
  # Encontrar índice del alineamiento donde H77 tiene ese residuo
  h77_index <- which(h77_vec != "-")[pos]
  return(seq_vec[h77_index])
}

############################################################
# 6. Analizar todas las secuencias
############################################################

results <- data.frame(
  Sequence = names(aln),
  stringsAsFactors = FALSE
)

for (pos in ras_positions) {
  results[[paste0("Pos_", pos)]] <- NA
  results[[paste0("RAS_", pos)]] <- FALSE
}

cat("Analizando secuencias...\n")

for (i in seq_along(aln)) {
  seq_name <- names(aln)[i]
  seq_vec  <- strsplit(as.character(aln[[i]]), "")[[1]]
  
  for (pos in ras_positions) {
    aa <- get_residue(seq_vec, h77_vec, pos)
    results[i, paste0("Pos_", pos)] <- aa
    
    # ¿Es una RAS conocida?
    if (aa %in% ras_catalog[[as.character(pos)]]) {
      results[i, paste0("RAS_", pos)] <- TRUE
    }
  }
}

############################################################
# 7. Guardar tabla de resultados
############################################################

write.csv(results, output_file, row.names = FALSE)

cat("✔ Tabla de RAS generada en:\n", output_file, "\n")
