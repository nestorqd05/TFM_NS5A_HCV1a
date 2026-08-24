############################################################
# 02_clean_nucleotides_strict.R
# Limpieza estricta de secuencias NS5A antes de traducir
############################################################

library(Biostrings)

input_file  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/raw/NS5A_1a_raw.fasta"
output_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed/NS5A_1a_clean_nt.fasta"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed",
           recursive = TRUE, showWarnings = FALSE)

cat("Leyendo secuencias...\n")
dna <- readDNAStringSet(input_file)
cat("Secuencias originales:", length(dna), "\n")

############################################################
# 1. Convertir todo a mayúsculas y eliminar caracteres invisibles
############################################################

dna_clean <- lapply(dna, function(seq){
  s <- toupper(as.character(seq))
  s <- gsub("[^ACGTN]", "N", s)   # reemplaza TODO lo que no sea A/C/G/T/N
  DNAString(s)
})

dna_clean <- DNAStringSet(dna_clean)

############################################################
# 2. Eliminar secuencias con demasiadas Ns
############################################################

max_Ns <- 150
dna_clean <- dna_clean[vcountPattern("N", dna_clean) < max_Ns]

cat("Secuencias tras limpieza estricta:", length(dna_clean), "\n")

############################################################
# 3. Guardar FASTA limpio
############################################################

writeXStringSet(dna_clean, output_file)
cat("✔ FASTA nucleotídico limpio guardado en:\n", output_file, "\n")

