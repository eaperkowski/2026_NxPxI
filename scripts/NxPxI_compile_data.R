# NxPxI_compile_data.R: Script that compiles primary datasheet.
# Main functions: clean raw LI-COR files, determine focal and
# chlorophyll leaf area, add focal/chlorophyll/Pfraction leaves
# into total leaf area calculation, fits A/Ci curves, and 
# extracts snapshot photosynthetic measurements. 
# Note: # all paths assume the folder containing this script 
# is the root directory

#####################################################################
# Libraries, data files, custom functions
#####################################################################
# Libraries
library(tidyverse)
library(LeafArea)
library(plantecophys)

# Load custom functions
R.utils::sourceDirectory("../functions/")

# Read relevant .csv files
tla <- read.csv("../data_sheets/NxPxI_tla.csv") %>%
  separate(id, sep = "_", into = c("inoc", "n_trt", "p_trt", "block", "ind"),
           remove = F) %>%
  mutate(n_trt = gsub("N", "", n_trt),
         p_trt = gsub("P", "", p_trt))

#####################################################################
# Determine chlorophyll and focal leaf areas
#####################################################################
ij_path <- "/Applications/ImageJ.app/"
ij_photos <- "../leaf_scans/actual_leaves/"

# Estimate leaf area
chl_focal_area <- run.ij(path.imagej = ij_path,
                         set.directory = ij_photos,
                         set.memory = 20,
                         distance.pixel = 118.017,
                         known.distance = 1, 
                         low.circ = 0.05,
                         trim.pixel = 0,
                         low.size = 0.1)

# Some light cleaning, add estimated Pfraction leaf area
chl_focal_area2 <- chl_focal_area %>%
  separate(sample, into = c("type", "id")) %>%
  pivot_wider(names_from = type, values_from = total.leaf.area) %>%
  mutate(pfraction_area_cm2 = (chl + focal) / 2) %>%
  dplyr::select(ind = id, chl_area_cm2 = chl,
                focal_area_cm2 = focal,
                pfraction_area_cm2)

# Merge individual leaves with rest of TLA, update TLA calculation
tla <- tla %>%
  full_join(chl_focal_area2, by = c("ind")) %>%
  mutate(tla_cm2 = tla_cm2 + chl_area_cm2 + focal_area_cm2 + pfraction_area_cm2)
# write.csv(tla, "../data_sheets/NxPxI_tla_full.csv", row.names = F)

########car#####################################################################
# Clean LI-COR files
#####################################################################

## Clean licor files
# clean_licor_files(directory_path = "../licor_data/licor_raw/",
#                   write_directory = "../licor_data/licor_cleaned/")

# Read cleaned LI-COR files
licor_files <- list.files(path = "../licor_data/licor_cleaned/",
                          recursive = T,
                          pattern = "\\.csv$",
                          full.names = T)
licor_files <- setNames(licor_files, 
                  stringr::str_extract(basename(licor_files),
                                       ".*(?=\\.csv)"))

# Read all files and merge into central data frame
licor_merged <- plyr::rbind.fill(lapply(licor_files, read.csv)) %>%
  mutate(date = lubridate::ymd_hms(date),
         date_only = stringr::word(date, 1)) %>%
  mutate(date = str_c(date_only, " ", hhmmss)) %>%
  mutate(id = ifelse(id == "3", "34", id),
         id = ifelse(id == "3_real", "3", id),
         Qin_cuvette = 2000,
         keep_row = "yes") %>%
  dplyr::select(obs, time, elapsed, date, date_only, hhmmss:Qin,
                Qin_cuvette, Qabs:SS_r, keep_row) %>%
  arrange(machine, date, id)

unique(licor_merged$id)

#####################################################################
# Fit LI-COR files
#####################################################################

# -----
# 001 
# -----
# Fit curve
aci_001 <- licor_merged %>% filter(id == "1") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_001)

# Create data frame
aci_coefs <- data.frame(id = "001", t(coef(aci_001)))

# -----
# 002
# -----
# Fit curve
aci_002 <- licor_merged %>% filter(id == "2") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_002)

# Add to data frame
aci_coefs[2,] <- c(id = "002", t(coef(aci_002)))

# -----
# 003
# -----
# Fit curve
aci_003 <- licor_merged %>% filter(id == "3") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_003)

# Create data frame
aci_coefs[3,] <- c(id = "003", t(coef(aci_003)))

# -----
# 004
# -----
# Fit curve
aci_004 <- licor_merged %>% filter(id == "4") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_004)

# Create data frame
aci_coefs[4,] <- c(id = "004", t(coef(aci_004)))

# -----
# 005
# -----
# Fit curve
aci_005 <- licor_merged %>% filter(id == "5") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_005)

# Create data frame
aci_coefs[5,] <- c(id = "005", t(coef(aci_005)))

# -----
# 006
# -----
# Fit curve
aci_006 <- licor_merged %>% filter(id == "6") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_006)

# Create data frame
aci_coefs[6,] <- c(id = "006", t(coef(aci_006)))

# -----
# 007
# -----
# Fit curve
aci_007 <- licor_merged %>% filter(id == "7") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_007)

# Create data frame
aci_coefs[7,] <- c(id = "007", t(coef(aci_007)))

# -----
# 008
# -----
# Fit curve
aci_008 <- licor_merged %>% filter(id == "8") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_008)

# Create data frame
aci_coefs[8,] <- c(id = "008", t(coef(aci_008)))

# -----
# 009
# -----
# Fit curve
aci_009 <- licor_merged %>% filter(id == "9") %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                         PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_009)

# Create data frame
aci_coefs[9,] <- c(id = "009", t(coef(aci_009)))

# -----
# 010
# -----
# Fit curve
licor_merged$keep_row[4769] <- "no"

aci_010 <- licor_merged %>% filter(id == "10" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_010)
summary(aci_010)

# Create data frame
aci_coefs[10,] <- c(id = "010", t(coef(aci_010)))

# -----
# 011
# -----
# Fit curve
aci_011 <- licor_merged %>% filter(id == "11" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_011)

# Create data frame
aci_coefs[11,] <- c(id = "011", t(coef(aci_011)))

# -----
# 012
# -----
# Fit curve
aci_012 <- licor_merged %>% filter(id == "12" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_012)

# Create data frame
aci_coefs[12,] <- c(id = "012", t(coef(aci_012)))

# -----
# 013
# -----
# Fit curve
aci_013 <- licor_merged %>% filter(id == "13" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_013)

# Create data frame
aci_coefs[13,] <- c(id = "013", t(coef(aci_013)))

# -----
# 014
# -----
# Fit curve
aci_014 <- licor_merged %>% filter(id == "14" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_014)

# Create data frame
aci_coefs[14,] <- c(id = "014", t(coef(aci_014)))

# -----
# 015
# -----
# Fit curve
aci_015 <- licor_merged %>% filter(id == "15" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_015)

# Create data frame
aci_coefs[15,] <- c(id = "015", t(coef(aci_015)))

# -----
# 016
# -----
# Fit curve
aci_016 <- licor_merged %>% filter(id == "16" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_016)

# Create data frame
aci_coefs[16,] <- c(id = "016", t(coef(aci_016)))

# -----
# 017
# -----
# Fit curve
aci_017 <- licor_merged %>% filter(id == "17" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_017)

# Create data frame
aci_coefs[17,] <- c(id = "017", t(coef(aci_017)))

# -----
# 018
# -----
# Fit curve
aci_018 <- licor_merged %>% filter(id == "18" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_018)

# Create data frame
aci_coefs[18,] <- c(id = "018", t(coef(aci_018)))

# -----
# 019
# -----
# Fit curve
aci_019 <- licor_merged %>% filter(id == "19" & A > 3)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_019)

# Create data frame
aci_coefs[19,] <- c(id = "019", t(coef(aci_019)))

# -----
# 020
# -----
# Fit curve
aci_020 <- licor_merged %>% filter(id == "20" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_020)

# Create data frame
aci_coefs[20,] <- c(id = "020", t(coef(aci_020)))

# -----
# 021
# -----
# Fit curve
aci_021 <- licor_merged %>% filter(id == "21" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_021)

# Create data frame
aci_coefs[21,] <- c(id = "021", t(coef(aci_021)))

# -----
# 022
# -----
# Fit curve
aci_022 <- licor_merged %>% filter(id == "22" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_022)

# Create data frame
aci_coefs[22,] <- c(id = "022", t(coef(aci_022)))

# -----
# 023
# -----
# Fit curve
aci_023 <- licor_merged %>% filter(id == "23" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_023)

# Create data frame
aci_coefs[23,] <- c(id = "023", NA, NA, NA, NA)

# -----
# 024
# -----
# Fit curve
aci_024 <- licor_merged %>% filter(id == "24" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_024)

# Create data frame
aci_coefs[24,] <- c(id = "024", t(coef(aci_024)))

# -----
# 025
# -----
# Fit curve
aci_025 <- licor_merged %>% filter(id == "25" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_025)

# Create data frame
aci_coefs[25,] <- c(id = "025", t(coef(aci_025)))

# -----
# 026
# -----
# Fit curve
aci_026 <- licor_merged %>% filter(id == "26" & A > 5)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 800)
plot(aci_026)

# Create data frame
aci_coefs[26,] <- c(id = "026", t(coef(aci_026)))

# -----
# 027
# -----
# Fit curve
aci_027 <- licor_merged %>% filter(id == "27" & A > 4)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 300)
plot(aci_027)

# Create data frame
aci_coefs[27,] <- c(id = "027", t(coef(aci_027)))

# -----
# 028
# -----
# Fit curve
aci_028 <- licor_merged %>% filter(id == "28" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_028)

# Create data frame
aci_coefs[28,] <- c(id = "028", t(coef(aci_028)))

# -----
# 029
# -----
# Fit curve
aci_029 <- licor_merged %>% filter(id == "29" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_029)

# Create data frame
aci_coefs[29,] <- c(id = "029", t(coef(aci_029)))

# -----
# 030
# -----
# Fit curve
aci_030 <- licor_merged %>% filter(id == "30" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_030)

# Create data frame
aci_coefs[30,] <- c(id = "030", t(coef(aci_030)))

## write.csv(aci_coefs, "../data_sheets/NxPxI_phys.csv", row.names = F)










