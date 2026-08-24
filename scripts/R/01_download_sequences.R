############################################################
# 01_download_sequences.R
# Descarga secuencias NS5A del HCV genotipo 1a (nucleótidos)
# Incluye la secuencia de referencia H77
############################################################

library(rentrez)
library(Biostrings)

# Ruta de salida (nucleótidos)
output_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/raw/NS5A_1a_raw.fasta"

# Crear carpeta si no existe
dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/raw",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 1. Búsqueda general de NS5A genotipo 1a en GenBank
############################################################

query <- "Hepatitis C virus[Organism] AND NS5A[Gene] AND 1a[All Fields]"

cat("Buscando secuencias NS5A 1a en GenBank...\n")
search_results <- entrez_search(db = "nuccore", term = query, retmax = 500)

cat("Secuencias encontradas:", search_results$count, "\n")

seq_ids <- search_results$ids

############################################################
# 2. Añadir explícitamente la referencia H77
# (ejemplo: accesión AF009606, ajusta si usas otra)
############################################################

ref_h77_id <- "AF009606"  # sustituye por el ID exacto que decidas usar

seq_ids <- unique(c(seq_ids, ref_h77_id))

cat("Total de IDs a descargar (incluyendo H77):", length(seq_ids), "\n")

############################################################
# 3. Descargar todas las secuencias en formato FASTA
############################################################

cat("Descargando secuencias...\n")
fasta_list <- lapply(seq_ids, function(id){
  entrez_fetch(db = "nuccore", id = id, rettype = "fasta")
})

############################################################
# 4. Guardar archivo FASTA nucleotídico
############################################################

writeLines(unlist(fasta_list), output_file)

cat("✔ Archivo FASTA nucleotídico generado:\n", output_file, "\n")
