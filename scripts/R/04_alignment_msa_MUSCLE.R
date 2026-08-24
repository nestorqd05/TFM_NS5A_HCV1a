############################################################
# 04_alignment_msa_MUSCLE.R
# Alineamiento múltiple de NS5A (proteína) usando MUSCLE
# Sin DECIPHER, sin MAFFT, sin msaConvert
############################################################

library(Biostrings)
library(msa)

############################################################
# 1. Rutas
############################################################

input_file  <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/processed/NS5A_1a_aa.fasta"
output_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/aligned/NS5A_1a_aligned.fasta"

dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/aligned",
           recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Leer secuencias proteicas
############################################################

cat("Leyendo secuencias proteicas...\n")
aa <- readAAStringSet(input_file)
cat("Total de secuencias:", length(aa), "\n")

if (!"H77_reference" %in% names(aa)) {
  stop("❌ ERROR: H77_reference no está en el FASTA proteico. Revisa el script 03.")
}

############################################################
# 3. Eliminar stop codons (*) para evitar avisos de MUSCLE
############################################################

aa_clean <- AAStringSet(gsub("\\*", "", as.character(aa)))
names(aa_clean) <- names(aa)

############################################################
# 4. Alineamiento con MUSCLE vía msa
############################################################

cat("Iniciando alineamiento con MUSCLE...\n")

aln <- msa(aa_clean,
           method = "Muscle",
           order = "input",
           type = "protein")

cat("✔ Alineamiento completado.\n")

############################################################
# 5. EXTRAER EL ALINEAMIENTO DIRECTAMENTE DEL OBJETO msa
#    aln@unmasked es un AAStringSet con el alineamiento final
############################################################

aligned <- aln@unmasked

############################################################
# 6. Guardar FASTA alineado
############################################################

writeXStringSet(aligned, output_file)

cat("✔ Archivo alineado guardado en:\n", output_file, "\n")
