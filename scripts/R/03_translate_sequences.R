############################################################
# 03_translate_sequences_safe.R
# Traducción segura de secuencias NS5A (nucleótidos → aminoácidos)
# - Detecta el mejor marco de lectura
# - Ignora secuencias corruptas
# - Añade la referencia H77 correctamente
############################################################

library(Biostrings)

############################################################
# 1. Rutas de entrada y salida
############################################################

# FASTA nucleotídico limpio generado por el script 02
input_file  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed/NS5A_1a_clean_nt.fasta"

# FASTA proteico que generará este script
output_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed/NS5A_1a_aa.fasta"

# FASTA con la referencia H77 (extraída del genoma AF009606)
h77_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/ref/H77_NS5A_nt.fasta"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer secuencias nucleotídicas limpias
############################################################

cat("Leyendo secuencias nucleotídicas limpias...\n")
dna <- readDNAStringSet(input_file)
cat("Total de secuencias:", length(dna), "\n")

############################################################
# 3. Función para traducir probando los 3 marcos de lectura
#    - Se prueban frames 0, 1 y 2
#    - Se elige el que produce menos stop codons (*)
#    - Si una secuencia no se puede traducir, devuelve NULL
############################################################

translate_safe <- function(seq) {
  frames <- list()
  
  # Probar los 3 marcos
  for (f in 0:2) {
    try({
      aa <- translate(subseq(seq, start = 1 + f))
      frames[[as.character(f)]] <- aa
    }, silent = TRUE)
  }
  
  # Eliminar traducciones fallidas
  frames <- frames[sapply(frames, length) > 0]
  if (length(frames) == 0) return(NULL)
  
  # Elegir el frame con menos stop codons
  stops <- sapply(frames, function(aa) sum(strsplit(as.character(aa), "")[[1]] == "*"))
  best <- frames[[which.min(stops)]]
  
  return(best)
}

############################################################
# 4. Traducir todas las secuencias
############################################################

cat("Traduciendo secuencias...\n")

aa_list <- list()
valid_names <- c()

for (i in seq_along(dna)) {
  aa <- translate_safe(dna[[i]])
  
  if (!is.null(aa)) {
    aa_list[[length(aa_list) + 1]] <- aa
    valid_names <- c(valid_names, names(dna)[i])
  }
}

aa <- AAStringSet(aa_list)
names(aa) <- valid_names

cat("Secuencias traducidas correctamente:", length(aa), "\n")

############################################################
# 5. Añadir la referencia H77
#    - Se traduce NS5A H77 desde su FASTA nucleotídico
#    - Se añade como última secuencia del FASTA proteico
############################################################

cat("Añadiendo referencia H77...\n")

h77_nt <- readDNAStringSet(h77_file)
h77_aa <- translate(h77_nt)

aa <- c(aa, AAStringSet(h77_aa))
names(aa)[length(aa)] <- "H77_reference"

cat("✔ Referencia H77 añadida correctamente.\n")

############################################################
# 6. Guardar FASTA proteico final
############################################################

writeXStringSet(aa, output_file)

cat("✔ Traducción completada.\n")
cat("Archivo generado:", output_file, "\n")
