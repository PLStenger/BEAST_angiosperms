library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)

setwd("~XXXXX/05_tree_affter_beast")

# 1. READING AND CORRECTION OF DENSITY TREES

# Run this first in your terminal
# sed 's/\[&[^]]*\]//g' XXXXXX/Murat_parameters_cleaned-SpeciesTreeAlignment_filtered.trees > XXXXXX/cleaned_trees.trees
# sed -i '' '/Translate/,/;/d' XXXXXX/cleaned_trees.trees

# Then:
trees_multi <- read.tree("cleaned_trees.trees", tree.names = TRUE)
tip_counts <- sapply(trees_multi, function(tree) length(tree$tip.label))
trees_filtered <- trees_multi[ sapply(trees_multi, function(tr) length(tr$tip.label) > 1) ]

file_path <- "Murat_parameters_cleaned_calibrated-SpeciesTreeAlignment_filtered.trees"
lines <- readLines(file_path)
start_translate <- grep("Translate", lines)
end_translate <- grep(";", lines[start_translate:length(lines)])[1] + start_translate - 1
translate_lines <- lines[(start_translate + 1):(end_translate - 1)]
translate_lines <- gsub(",", "", translate_lines)
translate_lines <- trimws(translate_lines) 
mapping <- sapply(strsplit(translate_lines, " "), function(x) x[2])
names(mapping) <- sapply(strsplit(translate_lines, " "), function(x) x[1])
replace_tips <- function(tree, map) {
  tree$tip.label <- unname(map[tree$tip.label])
  return(tree)
}
trees_multi_corrected <- lapply(trees_multi, replace_tips, map = mapping)
class(trees_multi_corrected) <- "multiPhylo"
trees_filtered <- trees_multi_corrected[sapply(trees_multi_corrected, function(tr) length(tr$tip.label) > 1)]

# 2. IMPORT TREE CONSENSUS MCC AND ALIGNMENT LABELS
library(treeio)
mcc_tree <- read.beast("Murat_parameters_cleaned_calibrated-SpeciesTreeAlignment_filtered_treeAnnotator.tree")

# 3. SUPERIMPOSED GRAPH OF THE TWO TREES
# Create the density tree
p_density <- ggdensitree(trees_filtered, alpha = 0.3, colour = 'steelblue')

# Create the MCC tree (with HPD data)
p_mcc_base <- ggtree(mcc_tree, colour = "black", size = 0.7)
p_mcc_base <- revts(p_mcc_base)  # inversion éventuelle du temps

# Retrieve data from MCC
mcc_data <- p_mcc_base$data

# Extract MCC + HPD data
mcc_data <- p_mcc_base$data
hpd_bounds <- do.call(rbind, mcc_data$height_0.95_HPD)
mcc_data$xmin <- -hpd_bounds[, 1]
mcc_data$xmax <- -hpd_bounds[, 2]

# Graph
p <- p_density +
  geom_tree(data = mcc_data, aes(x = x, y = y), color = "black", size = 0.7) +
  geom_segment(data = mcc_data,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               color = "blue", size = 2, alpha = 0.5) +
  geom_tiplab(data = mcc_data, size = 3) +
  theme_tree2() 

print(p)

