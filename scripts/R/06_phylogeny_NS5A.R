############################################################
# 06_phylogeny_NS5A.R
# Filogenia de NS5A usando alineamiento MUSCLE (script 04)
# - Árbol NJ con bootstrap
# - Exporta árbol en Newick y PNG
############################################################

library(Biostrings)
library(phangorn)
library(ape)

############################################################
# 1. Rutas
############################################################

aligned_file <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/data/aligned/NS5A_1a_aligned.fasta"

output_tree_newick <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_tree.newick"
output_tree_png     <- "C:/Users/nesto/Documents/TFM_NS5A_HCV1a/results/NS5A_tree.png"

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
# 3. Convertir alineamiento a formato phangorn
############################################################

cat("Convirtiendo alineamiento a formato phyDat...\n")

aln_phydat <- phyDat(as.matrix(aln), type = "AA")

############################################################
# 4. Distancias + Árbol NJ
############################################################

cat("Calculando matriz de distancias...\n")
dist_matrix <- dist.ml(aln_phydat)   # distancia máxima verosimilitud

cat("Construyendo árbol NJ...\n")
tree_nj <- NJ(dist_matrix)

############################################################
# 5. Bootstrap (100 réplicas)
############################################################

cat("Realizando bootstrap (100 réplicas)...\n")
bs <- bootstrap.phyDat(aln_phydat, FUN = function(x) NJ(dist.ml(x)), bs = 100)

# Añadir valores de bootstrap al árbol
tree_nj$node.label <- bs

############################################################
# 6. Guardar árbol en formato Newick
############################################################

write.tree(tree_nj, file = output_tree_newick)
cat("✔ Árbol guardado en formato Newick:\n", output_tree_newick, "\n")

############################################################
# 7. Guardar imagen PNG del árbol
############################################################

png(output_tree_png, width = 1200, height = 900)
plot(tree_nj, main = "Árbol filogenético NS5A (NJ + bootstrap)")
nodelabels(tree_nj$node.label, adj = c(1.2), frame = "none", cex = 0.7)
dev.off()

cat("✔ Imagen PNG del árbol guardada en:\n", output_tree_png, "\n")
