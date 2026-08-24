############################################################
# 10_NS5A_structure_RAS_mapping.R
# Mapeo estructural de RAS en NS5A GT1a usando modelo AlphaFold BFVD
############################################################

library(bio3d)
library(dplyr)

############################################################
# 1. Rutas
############################################################

pdb_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/structure/NS5A_GT1a_AF.pdb"
annot_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/structure/NS5A_RAS_annotations.txt"

############################################################
# 2. Leer estructura
############################################################

pdb <- read.pdb(pdb_file)

############################################################
# 3. Posiciones de RAS (H77 numbering)
############################################################

ras_positions <- c(28, 30, 31, 32, 58, 92, 93)

############################################################
# 4. Extraer coordenadas de las posiciones RAS
############################################################

coords <- pdb$atom %>%
  filter(resno %in% ras_positions, elety == "CA") %>%
  select(resno, x, y, z)

write.table(coords, annot_file, row.names = FALSE, quote = FALSE)

cat("✔ Coordenadas de RAS guardadas en:\n", annot_file, "\n")
cat("✔ Estructura lista para visualizar en PyMOL/ChimeraX.\n")
