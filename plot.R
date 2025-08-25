# (english translation at the end of the script)

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

# sed 's/\[&[^]]*\]//g' /Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/05_tree_affter_beast/Murat_parameters_cleaned_calibrated_all_prior-SpeciesTreeAlignment_filtered.trees > /Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/05_tree_affter_beast/cleaned_trees_lognormal.trees
# sed -i '' '/Translate/,/;/d' /Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/05_tree_affter_beast/cleaned_trees_lognormal.trees

# ================== 1) LECTURE & TRADUCTION DES ARBRES (POST-BEAST) ==================
# Arbres de la distribution (densitree)
trees_multi <- read.tree("cleaned_trees_lognormal.trees", tree.names = TRUE)

# Table de traduction (numérique -> nom d'espèce) extraite du bloc Translate du .trees
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
  # Remplace les identifiants numériques par les labels d'espèces
  # Les ids absents restent inchangés
  tree$tip.label <- unname(ifelse(tree$tip.label %in% names(map), map[tree$tip.label], tree$tip.label))
  return(tree)
}

trees_multi_corrected <- lapply(trees_multi, replace_tips, map = mapping)
class(trees_multi_corrected) <- "multiPhylo"

# ================== 2) PURGE DANS TOUS LES ARBRES + HARMONISATION DES TIPS ==================
# Drop Nymphaea_colorata dans chaque arbre
trees_purged <- lapply(trees_multi_corrected, function(tr) {
  if (!is.null(tr$tip.label) && target %in% tr$tip.label) {
    tr <- ape::drop.tip(tr, target)
  }
  tr
})

# Option: on retire les arbres dégénérés (< 2 tips)
trees_purged <- trees_purged[sapply(trees_purged, function(tr) length(tr$tip.label) > 1)]

# IMPORTANT : pour éviter les "traits fantômes", on force TOUS les arbres à partager EXACTEMENT le même jeu de tips
# -> on garde l'intersection des labels présents dans tous les arbres restants
if (length(trees_purged) < 2) {
  stop("Après purge, il reste moins de 2 arbres pour le densitree.")
}
common_tips <- Reduce(intersect, lapply(trees_purged, function(tr) tr$tip.label))

# Enlève aussi l'espèce cible si elle restait dans l'intersection par accident (sécurité)
common_tips <- setdiff(common_tips, target)

if (length(common_tips) < 2) {
  stop("Jeu de tips commun insuffisant (<2) après harmonisation. Vérifie les données sources.")
}

# Harmonise chaque arbre au même set de tips (ordre libre, mais identique pour tous)
trees_aligned <- lapply(trees_purged, function(tr) ape::keep.tip(tr, common_tips))

# Reclasser en multiPhylo
class(trees_aligned) <- "multiPhylo"

# Contrôles
cat("Nymphaea présente après purge (doit être 0) : ",
    sum(sapply(trees_aligned, function(tr) target %in% tr$tip.label)), "\n")
cat("Nombre d'arbres dans le densitree : ", length(trees_aligned), "\n")
cat("Nombre de tips communs : ", length(common_tips), "\n")

# ================== 3) MCC : LECTURE, PURGE, ALIGNEMENT AUX TIPS COMMUNS ==================
mcc_tree <- read.beast("Murat_parameters_cleaned_calibrated_all_prior-SpeciesTreeAlignment_filtered_TreeAnnotator.trees")

# Drop l'espèce cible dans le treedata (préserve les annotations HPD)
if (target %in% mcc_tree@phylo$tip.label) {
  mcc_tree <- treeio::drop.tip(mcc_tree, target)
}

# Aligne le MCC sur le même set de tips que le densitree (évite décalages y et “trous”)
extra_in_mcc <- setdiff(mcc_tree@phylo$tip.label, common_tips)
if (length(extra_in_mcc) > 0) {
  mcc_tree <- treeio::drop.tip(mcc_tree, extra_in_mcc)
}

# Si le MCC a des tips absents du densitree, ils ont été supprimés. On garantit la cohérence visuelle.

# ================== 4) GGTREE DU MCC & HPD ==================
p_mcc_base <- ggtree(mcc_tree, colour = "black", size = 0.7)
p_mcc_base <- revts(p_mcc_base)
mcc_data <- p_mcc_base$data

# Extraction robuste des HPD (gère les NAs)
hpd_bounds <- do.call(
  rbind,
  lapply(mcc_data$height_0.95_HPD, function(x) {
    if (is.null(x) || any(is.na(x))) c(NA_real_, NA_real_) else as.numeric(x)
  })
)
mcc_data$xmin <- -hpd_bounds[, 1]
mcc_data$xmax <- -hpd_bounds[, 2]

phylo_mcc <- as.phylo(mcc_tree)

# ================== 5) GROUPES/FAMILLES & COLORATION DES BRANCHES ==================
brassicaceae_species <- c("Arabidopsis_thaliana", "Arabidopsis_lyrata", "Arabidopsis_halleri",
                          "Capsella_rubella", "Thellungiella_parvul", "Brassica_oleracea",
                          "Brassica_napus_C", "Brassica_rapa", "Brassica_juncea", "Brassica_nigra")
rosaceae_species <- c("Rosa_chinensis", "Malus_domestica", "Fragaria_vesca",
                      "Pyrus_bretschneideri", "Prunus_persica", "Prunus_mume")

fabaceae_species <- c("Vigna_radiata", "Vigna_angularis", "Phaseolus_vulgaris", "Cajanus_cajan", "Lupinus_albus", "Glycine_max", "Medicago_truncatula", "Lotus_japonicus", "Lupinus_angustifoliu", "Arachis_duranensis", "Pisum_sativum", "Cicer_arietinum")

curcubitaceae_species <- c("Cucumis_sativus", "Cucumis_melo", "Lagenaria_siceraria", "Citrullus_lanatus", "Cucurbita_moschata", "Cucurbita_maxima", "Cucurbita_pepo_subsp")

solanaceae_species <- c("Solanum_melongena", "Solanum_lycopersicum", "Capsicum_annuum", "Nicotiana_Tabacum")

poaceae_species <- c("TdzA", "TdsA", "W6xA", "Triticum_urartu", "Hordeum_vulgare", "Aegilops_tauschii", "W6xD", "W6xB", "TdsB", "TdzB", "Sorghum_bicolor", "Miscanthus_sinensis", "Zea_mays_AGPv3", "Setaria_italica"   ,   "Panicum_virgatum"    , "Oropetium_thomaeum"  , "Oryza_sativa"     ,    "Brachypodium_distach")

mcc_data <- mcc_data %>%
  mutate(
    famille = case_when(
      label %in% brassicaceae_species ~ "Brassicaceae",
      label %in% rosaceae_species ~ "Rosaceae",
      label %in% fabaceae_species ~ "Fabaceae",
      label %in% curcubitaceae_species ~ "Curcubitaceae",
      label %in% solanaceae_species ~ "Solanaceae",
      label %in% poaceae_species ~ "Poaceae",
      TRUE ~ "Other"
    ),
    label_italic = paste0("italic('", label, "')"),
    couleur_hpd = case_when(
      label %in% brassicaceae_species ~ "goldenrod",
      label %in% rosaceae_species ~ "deeppink",
      label %in% fabaceae_species ~ "forestgreen",
      label %in% curcubitaceae_species ~ "blueviolet",
      label %in% solanaceae_species ~ "brown3",
      label %in% poaceae_species ~ "cornflowerblue",
      TRUE ~ "blue"
    )
  )

# noeuds des clades colorés
br_present <- brassicaceae_species[brassicaceae_species %in% phylo_mcc$tip.label]
ro_present <- rosaceae_species[rosaceae_species %in% phylo_mcc$tip.label]
fa_present <- fabaceae_species[fabaceae_species %in% phylo_mcc$tip.label]
cu_present <- curcubitaceae_species[curcubitaceae_species %in% phylo_mcc$tip.label]
so_present <- solanaceae_species[solanaceae_species %in% phylo_mcc$tip.label]
po_present <- poaceae_species[poaceae_species %in% phylo_mcc$tip.label]

clade_nodes_brassicaceae <- if (length(br_present) > 1) {
  mrca <- getMRCA(phylo_mcc, br_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

clade_nodes_rosaceae <- if (length(ro_present) > 1) {
  mrca <- getMRCA(phylo_mcc, ro_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

clade_nodes_fabaceae <- if (length(fa_present) > 1) {
  mrca <- getMRCA(phylo_mcc, fa_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

clade_nodes_curcubitaceae <- if (length(cu_present) > 1) {
  mrca <- getMRCA(phylo_mcc, cu_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

clade_nodes_solanaceae <- if (length(so_present) > 1) {
  mrca <- getMRCA(phylo_mcc, so_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

clade_nodes_poaceae <- if (length(po_present) > 1) {
  mrca <- getMRCA(phylo_mcc, po_present)
  c(mrca, unlist(phangorn::Descendants(phylo_mcc, mrca, type = "all")))
} else integer(0)

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

# ================== 6) GRAPHIQUE FINAL (DENSITREE + MCC + HPD + FOSSILES) ==================
x_breaks <- seq(-200, 0, by = 10)
x_labels <- seq(-200, 0, by = 10)

# Filtrer les lignes MCC qui ont de vraies HPD
mcc_hpd <- mcc_data %>% filter(!is.na(xmin), !is.na(xmax))

# Enlever les underscore des labels
mcc_data <- mcc_data %>%
  mutate(
    label_clean = gsub("_", " ", label),  # remplace _ par un espace dans les labels
    label_italic_clean = paste0("italic('", label_clean, "')")  # pour garder l’italique avec le nouveau label
  )

p <- ggdensitree(trees_aligned, alpha = 0.3, colour = "steelblue") +
  geom_tree(
    data = mcc_data,
    aes(x = x, y = y, color = branch_col),
    size = 0.7,
    inherit.aes = FALSE,
    lineend = "round"
  ) +
  geom_segment(
    data = mcc_hpd,
    aes(x = xmin, xend = xmax, y = y, yend = y, color = couleur_hpd),
    size = 2, alpha = 0.5, inherit.aes = FALSE
  ) +
  geom_tiplab(
    data = mcc_data,
    aes(label = label_italic_clean, color = famille),
    size = 3, parse = TRUE, inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Brassicaceae" = "goldenrod",
      "Rosaceae" = "deeppink",
      "Fabaceae" = "forestgreen",
      "Curcubitaceae" = "blueviolet",
      "Solanaceae" = "brown3",
      "Poaceae" = "cornflowerblue",
      "Other" = "black",
      "goldenrod" = "goldenrod",
      "deeppink" = "deeppink",
      "forestgreen" = "forestgreen",
      "blueviolet" = "blueviolet",
      "brown3" = "brown3",
      "cornflowerblue" = "cornflowerblue",
      "black" = "black",
      "blue" = "blue"
    ),
    guide = "none"
  ) +
  theme_tree2() +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, minor_breaks = NULL)

# Vitaceae -->  https://www.researchgate.net/profile/Steven-Manchester/publication/381882047_Cenozoic_seeds_of_Vitaceae_reveal_a_deep_history_of_extinction_and_dispersal_in_the_Neotropics/links/6745360bb5bd9d17d608579d/Cenozoic-seeds-of-Vitaceae-reveal-a-deep-history-of-extinction-and-dispersal-in-the-Neotropics.pdf
# Gentianales (contient coffee) --> https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0126690&type=printable
# Fagaceae -> http://www.fridgeirgrimsson.com/kerfi/wp-content/uploads/2017/02/2016_FG_GWG_RZ_TD.pdf

# Fossiles (lignes verticales + labels)
fossil_dates <- data.frame(
  date  = c(-65, -125, -73.5, -101.6, -34, -52, -60, -47, -82),
  label = c("Poaceae fossils", "CrownEudicots fossils", "Fabaceae fossils", "Rosaceae fossils", "Brassicaceae fossils", "Solanaceae fossils", "Vitaceae fossils", "Gentianales fossils", "Fagaceae fossils"),
  color = c("cornflowerblue", "grey50", "forestgreen", "deeppink", "goldenrod", "brown3", "black", "black", "black")
)
#fossil_dates$ypos <- max(mcc_data$y, na.rm = TRUE) * seq(0.90, 0.98, length.out = nrow(fossil_dates))
fossil_dates$ypos <- max(mcc_data$y, na.rm = TRUE) * 0.99


p <- p +
  geom_vline(data = fossil_dates, aes(xintercept = date, color = color),
             linetype = "dashed", size = 0.8) +
  geom_text(data = fossil_dates, aes(x = date, y = ypos, label = label, color = color),
            angle = 90, vjust = 1, hjust = 1, size = 4, inherit.aes = FALSE) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0.15)) # ajoute 10% d’espace à droite
  )


print(p)



################################################################################################
# Pour conniatre les noeufs et retrouver les bornes HPD
################################################################################################


# Ajouter les numéros de noeuds sur le graphique existant
p <- p +
  geom_text2(aes(label = node), hjust = -0.3, vjust = -0.5, size = 3, color = "darkred")


print(p)


# Fonction pour extraire les HPD (min, max) et la hauteur médiane pour un noeud
get_node_HPD <- function(treedata, node_id) {
  d <- treedata@data
  if (!node_id %in% d$node) {
    stop("Node ID not found in the tree data.")
  }
  
  node_info <- d[d$node == node_id, ]
  hpd <- unlist(node_info$height_0.95_HPD)
  
  return(list(
    node = node_id,
    median = node_info$height,
    hpd_min = hpd[1],
    hpd_max = hpd[2]
  ))
}

# Exemple d'utilisation :
get_node_HPD(mcc_tree, 135)


################################################################################################
# Fonction pour extraire les bornes HPD pour la MRCA de deux espèces données
################################################################################################

extract_HPD_MRCA <- function(tree, mcc_data, sp1, sp2) {
  # Vérifie si les espèces sont dans les tips
  if (!(sp1 %in% tree$tip.label)) {
    warning(paste("Espèce", sp1, "non trouvée dans l'arbre"))
    return(NULL)
  }
  if (!(sp2 %in% tree$tip.label)) {
    warning(paste("Espèce", sp2, "non trouvée dans l'arbre"))
    return(NULL)
  }
  
  # Trouve le noeud MRCA des deux espèces
  node_mrca <- ape::getMRCA(tree, c(sp1, sp2))
  
  # Extrait la ligne correspondante dans mcc_data (node = node_mrca)
  hpd_row <- mcc_data %>% filter(node == node_mrca)
  
  if (nrow(hpd_row) == 0) {
    warning(paste("MRCA non trouvé dans mcc_data pour", sp1, "et", sp2))
    return(NULL)
  }
  
  # Retourne les bornes HPD (corrigées pour l’échelle de temps)
  return(data.frame(
    species_pair = paste(sp1, "-", sp2),
    lower_HPD = -hpd_row$xmin,
    upper_HPD = -hpd_row$xmax
  ))
}

# Liste correcte des comparaisons
comp_esp <- list(
  c("Manihot_esculenta", "Populus_trichocarpa"),
  c("Gossypium_raimondii", "Theobroma_cacao"),
  c("Populus_trichocarpa", "Theobroma_cacao"),
  c("Prunus_persica", "Quercus_robur")
)

comp_esp <- list(
  c("Manihot_esculenta", "Populus_trichocarpa"),
  c("Gossypium_raimondii", "Theobroma_cacao"),
  c("Populus_trichocarpa", "Theobroma_cacao"),
  c("Brassica_nigra", "Brassica_rapa"),
  c("Brassica_rapa", "Thellungiella_parvula"),
  c("Arabidopsis_halleri", "Arabidopsis_lyrata"),
  c("Arabidopsis_halleri", "Arabidopsis_thaliana"),
  c("Arabidopsis_thaliana", "Capsella_rubella"),
  c("Arabidopsis_thaliana", "Brassica_rapa"),
  c("Arabidopsis_thaliana", "Carica_papaya"),
  c("Arabidopsis_thaliana", "Populus_trichocarpa"),
  c("Arabidopsis_thaliana", "Citrus_clementina"),
  c("Prunus_mume", "Prunus_persica"),
  c("Malus_domestica", "Pyrus_bretschneideri"),
  c("Malus_domestica", "Prunus_persica"),
  c("Fragaria_vesca", "Rosa_chinensis"),
  c("Fragaria_vesca", "Prunus_persica"),
  c("Prunus_persica", "Quercus_robur"),
  c("Arabidopsis_thaliana", "Prunus_persica"),
  c("Arabidopsis_thaliana", "Eucalyptus_grandis"),
  c("Vigna_angularis", "Vigna_radiata"),
  c("Phaseolus_vulgaris", "Vigna_radiata"),
  c("Glycine_max", "Phaseolus_vulgaris"),
  c("Cajanus_cajan", "Phaseolus_vulgaris"),
  c("Medicago_truncatula", "Pisum_sativum"),
  c("Cicer_arietinum", "Medicago_truncatula"),
  c("Lotus_japonicus", "Medicago_truncatula"),
  c("Medicago_truncatula", "Phaseolus_vulgaris"),
  c("Lupinus_albus", "Lupinus_angustifolius"),
  c("Lupinus_albus", "Phaseolus_vulgaris"),
  c("Arachis_duranensis", "Phaseolus_vulgaris"),
  c("Citrullus_lanatus", "Lagenaria_siceraria"),
  c("Cucumis_melo", "Cucumis_sativus"),
  c("Citrullus_lanatus", "Cucumis_melo"),
  c("Cucurbita_moschata", "Cucurbita_pepo_subsp_pepo"),
  c("Cucurbita_maxima", "Cucurbita_pepo_subsp_pepo"),
  c("Cucumis_melo", "Cucurbita_pepo_subsp_pepo"),
  c("Cucumis_melo", "Phaseolus_vulgaris"),
  c("Arabidopsis_thaliana", "Cucumis_melo"),
  c("Arabidopsis_thaliana", "Vitis_vinifera"),
  c("Solanum_lycopersicum", "Solanum_melongena"),
  c("Capsicum_annuum", "Solanum_lycopersicum"),
  c("Nicotiana_tabacum", "Solanum_lycopersicum"),
  c("Coffea_canephora", "Solanum_lycopersicum"),
  c("Sesamum_indicum", "Solanum_lycopersicum"),
  c("Daucus_carota", "Solanum_lycopersicum"),
  c("Helianthus_annuus", "Lactuca_sativa"),
  c("Lactuca_sativa", "Solanum_lycopersicum"),
  c("Arabidopsis_thaliana", "Solanum_lycopersicum"),
  c("Trochodendron_aralioides", "Tetracentron_sinense"),
  c("Buxus_sinica", "Tetracentron_sinense"),
  c("Buxus_sinica", "Nelumbo_nucifera"),
  c("Arabidopsis_thaliana", "Nelumbo_nucifera"),
  c("Arabidopsis_thaliana", "Ceratophyllum_demersum"),
  c("Cinnamomum_micranthum", "Liriodendron_chinense"),
  c("Arabidopsis_thaliana", "Liriodendron_chinense"),
  c("Aegilops_tauschii", "Triticum_urartu"),
  c("Aegilops_tauschii", "Hordeum_vulgare"),
  c("Aegilops_tauschii", "Brachypodium_distachyon"),
  c("Brachypodium_distachyon", "Oryza_sativa"),
  c("Miscanthus_sinensis", "Sorghum_bicolor"),
  c("Sorghum_bicolor", "Zea_mays"),
  c("Panicum_virgatum", "Setaria_italica"),
  c("Sorghum_bicolor", "Setaria_italica"),
  c("Oropetium_thomaeum", "Sorghum_bicolor"),
  c("Oryza_sativa", "Sorghum_bicolor"),
  c("Ananas_comosus", "Oryza_sativa"),
  c("Elaeis_guineensis", "Musa_acuminata"),
  c("Dioscorea_alata", "Elaeis_guineensis"),
  c("Asparagus_officinalis", "Elaeis_guineensis"),
  c("Elaeis_guineensis", "Oryza_sativa"),
  c("Colocasia_esculenta", "Spirodela_polyrhiza"),
  c("Oryza_sativa", "Spirodela_polyrhiza"),
  c("Acorus_tatarinowii", "Oryza_sativa"),
  c("Arabidopsis_thaliana", "Oryza_sativa"),
  c("Amborella_trichopoda", "Oryza_sativa"),
  c("Amborella_trichopoda", "Arabidopsis_thaliana"),
  c("Arabidopsis_thaliana", "Nymphaea_colorata"),
  c("Nymphaea_colorata", "Oryza_sativa"),
  c("Arabidopsis_thaliana", "Ginkgo_biloba"),
  c("Ginkgo_biloba", "Oryza_sativa"),
  c("Arabidopsis_thaliana", "Arabidopsis_thaliana"),
  c("Vitis_vinifera", "Vitis_vinifera"),
  c("Brachypodium_distachyon", "Brachypodium_distachyon"),
  c("Oryza_sativa", "Oryza_sativa"),
  c("Sorghum_bicolor", "Sorghum_bicolor"),
  c("Oryza_sativa", "Oryza_sativa"),
  c("Oryza_sativa", "Oryza_sativa")
)




# Appliquer et agréger les HPD
hpd_results <- do.call(rbind, lapply(comp_esp, function(pair) {
  extract_HPD_MRCA(phylo_mcc, mcc_data, pair[1], pair[2])
}))

print(hpd_results)




# Récupération du mu (taux de mutation)

# Installer et charger le package coda si besoin
#install.packages("coda")
library(coda)

# Charger le fichier .log dans R 
logfile <- read.table("Murat_parameters_cleaned_calibrated.log", header=TRUE, comment.char = "#")
logfile <- read.table("Murat_parameters_cleaned_calibrated_all_prior.log", header=TRUE, comment.char = "#")
logfile <- read.table("Murat_parameters_cleaned_calibrated_all_prior_all_normal.log", header=TRUE, comment.char = "#")
logfile <- read.table("brassicaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("rosaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("fabaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("curcubitaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("solanaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("poaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("curcubitaceae_fossils.log", header=TRUE, comment.char = "#")



# Vérifie les noms de colonnes
colnames(logfile)

# Sélectionne la colonne du taux de mutation (exemple: 'clock.rate')
# Pour horloge stricte
mu <- logfile$rate.mean

# Pour horloge détendue (UCLN) (horloge relaxée)
mu <- logfile$ucldMean

# Comme ce sont des Acides aminés
# UCRelaxedClockModel (un modèle à taux variables selon une distribution log-normale uncorrelated)
# Les taux BEAST ("mu") pour ces données d'acides aminés calibrées en millions d'années sont en substitutions/site/million d'années.
# Pour obtenir un taux par année, il faut diviser par un million.
mu_nucleo <- mu/1000000

# Statistiques descriptives
summary(mu)

# brassicaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001553 0.005451 0.006434 0.007637 0.007720 1.000000 

# rosaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001592 0.005805 0.011536 0.249213 0.091985 2.211061 

# fabaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002716 0.005414 0.006422 0.113997 0.061663 1.000000 

# curcubitaceae
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0160  0.5838  1.3496  2.8502  3.0764 62.2695 

# Solanaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001136 0.003211 0.003543 0.006687 0.003893 1.000000 

# poaceae
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05377 0.07567 0.08424 0.11692 0.09780 1.00000 


# All normal
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06410 0.06763 0.07118 0.07939 0.07429 1.00000 

# All lognormal
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1630  0.1919  0.2014  0.1999  0.2069  1.0000 



# Curcubitacee fossils
# > summary(mu)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0009531 0.0018509 0.0020697 0.0082197 0.0024814 1.0000000 



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
