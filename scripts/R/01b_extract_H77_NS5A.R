############################################################
# 01b_extract_H77_NS5A.R
# Extrae la región NS5A de la cepa H77 (AF009606 / NC_004102)
############################################################

library(rentrez)
library(Biostrings)

# Descargar genoma completo H77
cat("Descargando genoma H77...\n")
h77_genome <- entrez_fetch(db="nuccore", id="AF009606", rettype="fasta")

# Guardar temporalmente
temp_file <- tempfile(fileext=".fasta")
writeLines(h77_genome, temp_file)

# Leer secuencia
genome <- readDNAStringSet(temp_file)

# Coordenadas NS5A en H77
start <- 6258
end   <- 7603

ns5a_h77 <- subseq(genome[[1]], start=start, end=end)

# Guardar FASTA
output_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/ref/H77_NS5A_nt.fasta"
dir.create("C:/Users/nesto/Documents/TFM_NS5A_HCV1a/ref", recursive=TRUE, showWarnings=FALSE)

writeXStringSet(DNAStringSet(ns5a_h77), output_file, format="fasta")

cat("✔ NS5A H77 extraída y guardada en:\n", output_file, "\n")
