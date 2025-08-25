# ================== PACKAGES ==================
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)

# ================== PARAMS ====================
setwd("~/Documents/INRAE_PaleoLab/03_angiosperms_origins/05_tree_affter_beast")
target <- "Nymphaea_colorata"

# sed 's/\[&[^]]*\]//g' 03_angiosperms_origins/05_tree_affter_beast/Murat_parameters_cleaned_calibrated_all_prior-SpeciesTreeAlignment_filtered.trees > 03_angiosperms_origins/05_tree_affter_beast/cleaned_trees_lognormal.trees
# sed -i '' '/Translate/,/;/d' 03_angiosperms_origins/05_tree_affter_beast/cleaned_trees_lognormal.trees

# ================== 1) READING & TRANSLATING TREES (POST-BEAST) ==================
# Trees from the posterior distribution (densitree)
trees_multi <- read.tree("cleaned_trees_lognormal.trees", tree.names = TRUE)

# Translation table (numeric -> species name) extracted from the Translate block of the .trees
file_path <- "Murat_parameters_cleaned_calibrated_all_prior-SpeciesTreeAlignment_filtered.trees"
lines <- readLines(file_path)
start_translate <- grep("\\bTranslate\\b", lines)
end_translate <- grep(";", lines[start_translate:length(lines)])[1] + start_translate - 1
translate_lines <- lines[(start_translate + 1):(end_translate - 1)]
translate_lines <- gsub(",", "", translate_lines)
translate_lines <- trimws(translate_lines)

parts <- strsplit(translate_lines, "\\s+")
ids   <- sapply(parts, `[[`, 1)
labs  <- sapply(parts, `[[`, 2)
mapping <- setNames(labs, ids)

replace_tips <- function(tree, map) {
  # Replace numeric IDs with species labels
  # IDs not found remain unchanged
  tree$tip.label <- unname(ifelse(tree$tip.label %in% names(map), map[tree$tip.label], tree$tip.label))
  return(tree)
}

trees_multi_corrected <- lapply(trees_multi, replace_tips, map = mapping)
class(trees_multi_corrected) <- "multiPhylo"

# ================== 2) PURGE ACROSS ALL TREES + HARMONIZATION OF TIPS ==================
# Drop Nymphaea_colorata in each tree
trees_purged <- lapply(trees_multi_corrected, function(tr) {
  if (!is.null(tr$tip.label) && target %in% tr$tip.label) {
    tr <- ape::drop.tip(tr, target)
  }
  tr
})

# Option: drop degenerate trees (< 2 tips)
trees_purged <- trees_purged[sapply(trees_purged, function(tr) length(tr$tip.label) > 1)]

# IMPORTANT: to avoid "ghost lineages", force ALL trees to share EXACTLY the same set of tips
# -> keep only the intersection of labels present in all remaining trees
if (length(trees_purged) < 2) {
  stop("After purge, fewer than 2 trees remain for densitree.")
}
common_tips <- Reduce(intersect, lapply(trees_purged, function(tr) tr$tip.label))

# Also remove the target species from the intersection if it is there by accident (safety)
common_tips <- setdiff(common_tips, target)

if (length(common_tips) < 2) {
  stop("Common tip set insufficient (<2) after harmonization. Check source data.")
}

# Harmonize each tree to the same tip set (order free, but identical across trees)
trees_aligned <- lapply(trees_purged, function(tr) ape::keep.tip(tr, common_tips))

# Set back as multiPhylo
class(trees_aligned) <- "multiPhylo"

# Controls
cat("Nymphaea present after purge (should be 0): ",
    sum(sapply(trees_aligned, function(tr) target %in% tr$tip.label)), "\n")
cat("Number of trees in densitree: ", length(trees_aligned), "\n")
cat("Number of common tips: ", length(common_tips), "\n")

# ================== 3) MCC TREE : READING, PURGE, ALIGNMENT TO COMMON TIPS ==================
mcc_tree <- read.beast("Murat_parameters_cleaned_calibrated_all_prior-SpeciesTreeAlignment_filtered_TreeAnnotator.trees")

# Drop target species in treedata (keeps HPD annotations)
if (target %in% mcc_tree@phylo$tip.label) {
  mcc_tree <- treeio::drop.tip(mcc_tree, target)
}

# Align MCC to same tip set as densitree (avoid y-shifts and “holes”)
extra_in_mcc <- setdiff(mcc_tree@phylo$tip.label, common_tips)
if (length(extra_in_mcc) > 0) {
  mcc_tree <- treeio::drop.tip(mcc_tree, extra_in_mcc)
}

# If MCC had tips absent from densitree, they are removed. Ensures visual consistency.

# ================== 4) GGTREE MCC & HPD ==================
p_mcc_base <- ggtree(mcc_tree, colour = "black", size = 0.7)
p_mcc_base <- revts(p_mcc_base)
mcc_data <- p_mcc_base$data

# Robust extraction of HPDs (handles NAs)
hpd_bounds <- do.call(
  rbind,
  lapply(mcc_data$height_0.95_HPD, function(x) {
    if (is.null(x) || any(is.na(x))) c(NA_real_, NA_real_) else as.numeric(x)
  })
)
mcc_data$xmin <- -hpd_bounds[, 1]
mcc_data$xmax <- -hpd_bounds[, 2]

phylo_mcc <- as.phylo(mcc_tree)

# ================== 5) GROUPS/FAMILIES & BRANCH COLORING ==================
# (... unchanged lists of species)

mcc_data <- mcc_data %>%
  mutate(
    family = case_when(
      label %in% brassicaceae_species ~ "Brassicaceae",
      label %in% rosaceae_species ~ "Rosaceae",
      label %in% fabaceae_species ~ "Fabaceae",
      label %in% curcubitaceae_species ~ "Curcubitaceae",
      label %in% solanaceae_species ~ "Solanaceae",
      label %in% poaceae_species ~ "Poaceae",
      TRUE ~ "Other"
    ),
    label_italic = paste0("italic('", label, "')"),
    hpd_color = case_when(
      label %in% brassicaceae_species ~ "goldenrod",
      label %in% rosaceae_species ~ "deeppink",
      label %in% fabaceae_species ~ "forestgreen",
      label %in% curcubitaceae_species ~ "blueviolet",
      label %in% solanaceae_species ~ "brown3",
      label %in% poaceae_species ~ "cornflowerblue",
      TRUE ~ "blue"
    )
  )

# nodes of colored clades
# (... unchanged clade_nodes code)

mcc_data <- mcc_data %>%
  mutate(
    branch_col = case_when(
      node %in% clade_nodes_brassicaceae ~ "goldenrod",
      node %in% clade_nodes_rosaceae ~ "deeppink",
      node %in% clade_nodes_fabaceae ~ "forestgreen",
      node %in% clade_nodes_curcubitaceae ~ "blueviolet",
      node %in% clade_nodes_solanaceae ~ "brown3",
      node %in% clade_nodes_poaceae ~ "cornflowerblue",
      TRUE ~ "black"
    )
  )

# ================== 6) FINAL PLOT (DENSITREE + MCC + HPD + FOSSILS) ==================
x_breaks <- seq(-200, 0, by = 10)
x_labels <- seq(-200, 0, by = 10)

# Keep only MCC lines that have valid HPD
mcc_hpd <- mcc_data %>% filter(!is.na(xmin), !is.na(xmax))

# Remove underscores in labels
mcc_data <- mcc_data %>%
  mutate(
    label_clean = gsub("_", " ", label),  # replace _ with space
    label_italic_clean = paste0("italic('", label_clean, "')")  # keep italics with new label
  )

p <- ggdensitree(trees_aligned, alpha = 0.3, colour = "steelblue") +
  geom_tree(...) + # unchanged graphics code
  geom_segment(...) +
  geom_tiplab(
    data = mcc_data,
    aes(label = label_italic_clean, color = family),
    size = 3, parse = TRUE, inherit.aes = FALSE
  ) +
  scale_color_manual(...) +
  theme_tree2() +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, minor_breaks = NULL)

# Fossils (vertical lines + labels)
fossil_dates <- data.frame(
  date  = c(-65, -125, -73.5, -101.6, -34, -52, -60, -47, -82),
  label = c("Poaceae fossils", "CrownEudicots fossils", "Fabaceae fossils",
            "Rosaceae fossils", "Brassicaceae fossils", "Solanaceae fossils",
            "Vitaceae fossils", "Gentianales fossils", "Fagaceae fossils"),
  color = c("cornflowerblue", "grey50", "forestgreen", "deeppink",
            "goldenrod", "brown3", "black", "black", "black")
)
fossil_dates$ypos <- max(mcc_data$y, na.rm = TRUE) * 0.99

p <- p +
  geom_vline(...) +
  geom_text(...) +
  scale_x_continuous(...)

print(p)

################################################################################################
# To identify the nodes and retrieve HPD intervals
################################################################################################

# Add node numbers to the existing plot
p <- p + geom_text2(aes(label = node), hjust = -0.3, vjust = -0.5, size = 3, color = "darkred")
print(p)

# Function to extract HPD (min, max) and median height for a node
get_node_HPD <- function(treedata, node_id) {
  ...
}

# Example:
get_node_HPD(mcc_tree, 135)

################################################################################################
# Function to extract HPD intervals for the MRCA of two given species
################################################################################################

extract_HPD_MRCA <- function(tree, mcc_data, sp1, sp2) {
  # Check if species are in tips
  if (!(sp1 %in% tree$tip.label)) {
    warning(paste("Species", sp1, "not found in tree"))
    return(NULL)
  }
  ...
  
  # Return HPD bounds (corrected for timescale)
  return(data.frame(
    species_pair = paste(sp1, "-", sp2),
    lower_HPD = -hpd_row$xmin,
    upper_HPD = -hpd_row$xmax
  ))
}

# species comparisons list
comp_esp <- list(
  ...
)

# Apply and aggregate HPDs
hpd_results <- do.call(rbind, lapply(comp_esp, function(pair) {
  extract_HPD_MRCA(phylo_mcc, mcc_data, pair[1], pair[2])
}))
print(hpd_results)


# ================== 7) ESTIMATING mu (mutation rate) ==================
library(coda)

# Load the .log files into R
logfile <- read.table("Murat_parameters_cleaned_calibrated.log", header=TRUE, comment.char = "#")
# (other log files loaded similarly...)

# Check column names
colnames(logfile)

# Select mutation rate column (example: 'clock.rate')
# Strict clock
mu <- logfile$rate.mean

# Relaxed clock (UCLN)
mu <- logfile$ucldMean

# Since these are amino acids:
# UCRelaxedClockModel with uncorrelated lognormal distribution
# BEAST rates ("mu") are substitutions/site/MY.
# To obtain per year rate, divide by 1e6.
mu_nucleo <- mu/1000000

# Descriptive statistics
summary(mu)
# (comments with stats remain unchanged)
