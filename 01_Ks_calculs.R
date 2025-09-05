library(tidyverse)
library(boot)

# Charger le fichier
ks <- read.delim("/Users/stengerpierre-louis/Documents/INRAE_PaleoLab/03_angiosperms_origins/01_data/Ks_results/gene_pairs_022023.Ks_results.tsv", header = TRUE)

# Restreindre l'intervalle de Ks
ks.sp <- ks %>% filter(Ks_PAML > 0.01 & Ks_PAML < 3)

# Créer la colonne species_pair
ks.sp <- ks.sp %>% mutate(species_pair = paste(Sp1, Sp2, sep = "_"))

# === Code pour une ou plusieurs comparaisons ciblées ===

# Spécifier les comparaisons à tester, par exemple Os_Os
species_pairs_to_test <- c("Vv_Vv") # Ajoute d'autres comparaisons dans ce vecteur si besoin

# Fonction pour estimer le pic de la densité
find_peak <- function(data) {
  dens <- density(data$Ks_PAML, na.rm = TRUE)
  dens$x[which.max(dens$y)]
}

# Fonction de statistique à passer au bootstrap
statistic_function <- function(data, indices) {
  boot_data <- data[indices, ]
  find_peak(boot_data)
}

# Résultat final
results <- data.frame(species_pair = character(),
                      peak = numeric(),
                      lower_bound = numeric(),
                      upper_bound = numeric(),
                      stringsAsFactors = FALSE)

for (species in species_pairs_to_test) {
  filtered_data <- ks.sp %>% filter(species_pair == species)
  if (nrow(filtered_data) > 0) {
    b1 <- boot(data = filtered_data, statistic = statistic_function, R = 10000)
    boot_ci <- boot.ci(b1, type = "bca")
    results <- rbind(results, data.frame(
      species_pair = species,
      peak = find_peak(filtered_data),
      lower_bound = boot_ci$bca[4],
      upper_bound = boot_ci$bca[5]
    ))
  }
}

print(results)

