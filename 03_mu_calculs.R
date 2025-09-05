# ================== PACKAGES ==================
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)

# Récupération du mu (taux de mutation)

# Installer et charger le package coda si besoin
#install.packages("coda")
library(coda)

# Charger le fichier .log dans R 
#logfile <- read.table("Murat_parameters_cleaned_calibrated.log", header=TRUE, comment.char = "#")
#logfile <- read.table("Murat_parameters_cleaned_calibrated_all_prior.log", header=TRUE, comment.char = "#")
#logfile <- read.table("Murat_parameters_cleaned_calibrated_all_prior_all_normal.log", header=TRUE, comment.char = "#")
logfile <- read.table("brassicaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("rosaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("fabaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("curcubitaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("curcubitaceae_fossils.log", header=TRUE, comment.char = "#")
logfile <- read.table("solanaceae.log", header=TRUE, comment.char = "#")
logfile <- read.table("poaceae.log", header=TRUE, comment.char = "#")

logfile <- read.table("Murat_parameters_cleaned_calibrated_all_prior.log", header=TRUE, comment.char = "#")
logfile <- read.table("Gymno.log", header=TRUE, comment.char = "#")

logfile <- read.table("Grand_Teta.log", header=TRUE, comment.char = "#")
logfile <- read.table("rho.log", header=TRUE, comment.char = "#")
logfile <- read.table("gamma.log", header=TRUE, comment.char = "#")


mu <- logfile$ucldMean
summary(mu)


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
# NEW:
# > summary(mu)                                                               
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001540 0.005416 0.006410 0.007464 0.007687 1.000000

# rosaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001592 0.005805 0.011536 0.249213 0.091985 2.211061 
# NEW:
#> summary(mu)                                                           
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001295 0.002542 0.003005 0.121085 0.010661 2.211061 

# fabaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002716 0.005414 0.006422 0.113997 0.061663 1.000000 
# NEW:
#> summary(mu)                                                           
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001541 0.004596 0.005051 0.022082 0.005708 1.000000 

# curcubitaceae
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0160  0.5838  1.3496  2.8502  3.0764 62.2695 

# Solanaceae
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001136 0.003211 0.003543 0.006687 0.003893 1.000000 

# poaceae
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05377 0.07567 0.08424 0.11692 0.09780 1.00000 

# new poacee:
# > summary(mu)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02645 0.04252 0.04714 0.05210 0.05196 1.00000 


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
# NEW:
# > summary(mu)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008717 0.0017960 0.0019898 0.0042977 0.0022478 1.0000000 



# Murat_parameters_cleaned_calibrated_all_prior.log
# > summary(mu)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008717 0.0017960 0.0019898 0.0042977 0.0022478 1.0000000 







# Grand_Teta.log

# rho.log

# gamma.log
# > summary(mu)                                                        
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1525  0.1770  0.2581  0.3139  0.4344  1.0000 




