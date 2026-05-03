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

# Read .csv files
## Load file with only total leaf area (to be merged with focal,
## chlorophyll, and Pfraction leaves below) 
# tla <- read.csv("../data_sheets/NxPxI_tla.csv") %>%
#   separate(id, sep = "_", into = c("inoc", "n_trt", "p_trt", "block", "ind"),
#            remove = F) %>%
#   mutate(n_trt = gsub("N", "", n_trt),
#          p_trt = gsub("P", "", p_trt))

# Load file with total leaf area, includnig focal, chlorophyll, and
# P fraction leaves
tla_full <- read.csv("../data_sheets/NxPxI_tla_full.csv") %>%
  mutate(ind = str_pad(ind, width = 3, pad = "0"))

#####################################################################
# Determine chlorophyll and focal leaf areas
#####################################################################
## Image J paths and scan directory
# ij_path <- "/Applications/ImageJ.app/"
# ij_photos <- "../leaf_scans/actual_leaves/"

## Estimate leaf area
# chl_focal_area <- run.ij(path.imagej = ij_path,
#                          set.directory = ij_photos,
#                          set.memory = 20,
#                          distance.pixel = 118.017,
#                          known.distance = 1, 
#                          low.circ = 0.05,
#                          trim.pixel = 0,
#                          low.size = 0.1)

## Some light cleaning, add estimated Pfraction leaf area as average of 
## the focal and chlorophyll leaves (can do this given nonindependence of
## trifoliate leaf set)
# chl_focal_area2 <- chl_focal_area %>%
#   separate(sample, into = c("type", "id")) %>%
#   pivot_wider(names_from = type, values_from = total.leaf.area) %>%
#   mutate(pfraction_area_cm2 = (chl + focal) / 2) %>%
#   dplyr::select(ind = id, chl_area_cm2 = chl,
#                 focal_area_cm2 = focal,
#                 pfraction_area_cm2)

## Merge individual leaves with rest of TLA, update TLA calculation
# tla <- tla %>%
#   full_join(chl_focal_area2, by = c("ind")) %>%
#   mutate(tla_cm2 = tla_cm2 + chl_area_cm2 + focal_area_cm2 + pfraction_area_cm2)
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
aci_coefs[22,] <- c(id = "022", t(coef(aci_022)))

# -----
# 023
# -----

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
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

# Add to data frame
aci_coefs[30,] <- c(id = "030", t(coef(aci_030)))

# -----
# 031
# -----
# Fit curve
aci_031 <- licor_merged %>% filter(id == "31" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_031)

# Add to data frame
aci_coefs[31,] <- c(id = "031", t(coef(aci_031)))


# -----
# 032
# -----
# Fit curve
aci_032 <- licor_merged %>% filter(id == "32" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_032)

# Add to data frame
aci_coefs[32,] <- c(id = "032", t(coef(aci_032)))

# -----
# 033
# -----
# Fit curve
aci_033 <- licor_merged %>% filter(id == "33" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_033)

# Add to data frame
aci_coefs[33,] <- c(id = "033", t(coef(aci_033)))

# -----
# 034
# -----
# Fit curve
aci_034 <- licor_merged %>% filter(id == "34" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_034)

# Add to data frame
aci_coefs[34,] <- c(id = "034", t(coef(aci_034)))

# -----
# 035
# -----
# Fit curve
aci_035 <- licor_merged %>% filter(id == "35" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_035)

# Add to data frame
aci_coefs[35,] <- c(id = "035", t(coef(aci_035)))

# -----
# 036
# -----
# Fit curve
aci_036 <- licor_merged %>% filter(id == "36" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_036)

# Add to data frame
aci_coefs[36,] <- c(id = "036", t(coef(aci_036)))

# -----
# 037
# -----
# Fit curve
aci_037 <- licor_merged %>% filter(id == "37" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_037)

# Add to data frame
aci_coefs[37,] <- c(id = "037", t(coef(aci_037)))

# -----
# 038
# -----
# Fit curve
aci_038 <- licor_merged %>% filter(id == "38" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_038)

# Add to data frame
aci_coefs[38,] <- c(id = "038", t(coef(aci_038)))

# -----
# 039
# -----
# Fit curve
aci_039 <- licor_merged %>% filter(id == "39" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_039)

# Add to data frame
aci_coefs[39,] <- c(id = "039", t(coef(aci_039)))

# -----
# 040
# -----
# Fit curve
aci_040 <- licor_merged %>% filter(id == "40" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_040)

# Add to data frame
aci_coefs[40,] <- c(id = "040", t(coef(aci_040)))

# -----
# 041
# -----
# Fit curve
## No curve

# Add to data frame
aci_coefs[41,] <- c(id = "041", NA, NA, NA, NA)

# -----
# 042
# -----
# Fit curve
aci_042 <- licor_merged %>% filter(id == "42" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_042)

# Add to data frame
aci_coefs[42,] <- c(id = "042", t(coef(aci_042)))

# -----
# 043
# -----
# Fit curve
aci_043 <- licor_merged %>% filter(id == "43" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_043)

# Add to data frame
aci_coefs[43,] <- c(id = "043", t(coef(aci_043)))

# -----
# 044
# -----
# Fit curve
aci_044 <- licor_merged %>% filter(id == "44" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_044)

# Add to data frame
aci_coefs[44,] <- c(id = "044", t(coef(aci_044)))

# -----
# 045
# -----
# Fit curve
aci_045 <- licor_merged %>% filter(id == "45" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_045)

# Add to data frame
aci_coefs[45,] <- c(id = "045", t(coef(aci_045)))

# -----
# 046
# -----
# Fit curve
aci_046 <- licor_merged %>% filter(id == "46" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_046)

# Add to data frame
aci_coefs[46,] <- c(id = "046", t(coef(aci_046)))

# -----
# 047
# -----
# Fit curve
aci_047 <- licor_merged %>% filter(id == "47" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_047)

# Add to data frame
aci_coefs[47,] <- c(id = "047", t(coef(aci_047)))

# -----
# 048
# -----
# Fit curve
aci_048 <- licor_merged %>% filter(id == "48" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_048)

# Add to data frame
aci_coefs[48,] <- c(id = "048", t(coef(aci_048)))

# -----
# 049
# -----
# Fit curve
aci_049 <- licor_merged %>% filter(id == "49" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_049)

# Add to data frame
aci_coefs[49,] <- c(id = "049", t(coef(aci_049)))

# -----
# 050
# -----
# Fit curve
aci_050 <- licor_merged %>% filter(id == "50" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_050)

# Add to data frame
aci_coefs[50,] <- c(id = "050", t(coef(aci_050)))

# -----
# 051
# -----

# Add to data frame
aci_coefs[51,] <- c(id = "051", NA, NA, NA, NA)

# -----
# 052
# -----
# Fit curve
aci_052 <- licor_merged %>% filter(id == "52" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_052)

# Add to data frame
aci_coefs[52,] <- c(id = "052", t(coef(aci_052)))

# -----
# 053
# -----
# Fit curve
aci_053 <- licor_merged %>% filter(id == "53" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_053)

# Add to data frame
aci_coefs[53,] <- c(id = "053", t(coef(aci_053)))

# -----
# 054
# -----
# Fit curve
aci_054 <- licor_merged %>% filter(id == "54" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_054)

# Add to data frame
aci_coefs[54,] <- c(id = "054", t(coef(aci_054)))

# -----
# 055
# -----
# Fit curve
aci_055 <- licor_merged %>% filter(id == "55" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_055)

# Add to data frame
aci_coefs[55,] <- c(id = "055", t(coef(aci_055)))

# -----
# 056
# -----
# Fit curve
aci_056 <- licor_merged %>% filter(id == "56" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_056)

# Add to data frame
aci_coefs[56,] <- c(id = "056", t(coef(aci_056)))

# -----
# 057
# -----
# Fit curve
aci_057 <- licor_merged %>% filter(id == "57" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_057)

# Add to data frame
aci_coefs[57,] <- c(id = "057", t(coef(aci_057)))

# -----
# 058
# -----
# Fit curve
aci_058 <- licor_merged %>% filter(id == "58" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_058)

# Add to data frame
aci_coefs[58,] <- c(id = "058", t(coef(aci_058)))

# -----
# 059
# -----
# Fit curve
aci_059 <- licor_merged %>% filter(id == "59" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_059)

# Add to data frame
aci_coefs[59,] <- c(id = "059", t(coef(aci_059)))

# -----
# 060
# -----
# Fit curve
aci_060 <- licor_merged %>% filter(id == "60" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_060)

# Add to data frame
aci_coefs[60,] <- c(id = "060", t(coef(aci_060)))

# -----
# 061
# -----
# Fit curve
aci_061 <- licor_merged %>% filter(id == "61" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_061)

# Add to data frame
aci_coefs[61,] <- c(id = "061", t(coef(aci_061)))

# -----
# 062
# -----
# Fit curve
aci_062 <- licor_merged %>% filter(id == "62" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_062)

# Add to data frame
aci_coefs[62,] <- c(id = "062", t(coef(aci_062)))

# -----
# 063
# -----
# Fit curve
aci_063 <- licor_merged %>% filter(id == "63" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_063)

# Add to data frame
aci_coefs[63,] <- c(id = "063", t(coef(aci_063)))

# -----
# 064
# -----
# Fit curve
aci_064 <- licor_merged %>% filter(id == "64" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 500)
plot(aci_064)

# Add to data frame
aci_coefs[64,] <- c(id = "064", t(coef(aci_064)))

# -----
# 065
# -----
# Fit curve
aci_065 <- licor_merged %>% filter(id == "65_b" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_065)

# Add to data frame
aci_coefs[65,] <- c(id = "065", t(coef(aci_065)))

# -----
# 066
# -----
# Fit curve
aci_066 <- licor_merged %>% filter(id == "66" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_066)

# Add to data frame
aci_coefs[66,] <- c(id = "066", t(coef(aci_066)))

# -----
# 067
# -----
# Fit curve
aci_067 <- licor_merged %>% filter(id == "67" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_067)

# Add to data frame
aci_coefs[67,] <- c(id = "067", t(coef(aci_067)))

# -----
# 068
# -----
# Fit curve
aci_068 <- licor_merged %>% filter(id == "68" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_068)

# Add to data frame
aci_coefs[68,] <- c(id = "068", t(coef(aci_068)))

# -----
# 069
# -----
# Fit curve
aci_069 <- licor_merged %>% filter(id == "69" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_069)

# Add to data frame
aci_coefs[69,] <- c(id = "069", t(coef(aci_069)))

# -----
# 070
# -----
# Fit curve
aci_070 <- licor_merged %>% filter(id == "70" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_070)

# Add to data frame
aci_coefs[70,] <- c(id = "070", t(coef(aci_070)))

# -----
# 071
# -----
# Fit curve
aci_071 <- licor_merged %>% filter(id == "71" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_071)

# Add to data frame
aci_coefs[71,] <- c(id = "071", t(coef(aci_071)))

# -----
# 072
# -----
# Fit curve
aci_072 <- licor_merged %>% filter(id == "72" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_072)

# Add to data frame
aci_coefs[72,] <- c(id = "072", t(coef(aci_072)))

# -----
# 073
# -----
# Fit curve
aci_073 <- licor_merged %>% filter(id == "73" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_073)

# Add to data frame
aci_coefs[73,] <- c(id = "073", t(coef(aci_073)))

# -----
# 074
# -----
# Fit curve
aci_074 <- licor_merged %>% filter(id == "74" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_074)

# Add to data frame
aci_coefs[74,] <- c(id = "074", t(coef(aci_074)))

# -----
# 075
# -----
# Fit curve
aci_075 <- licor_merged %>% filter(id == "75" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_075)

# Add to data frame
aci_coefs[75,] <- c(id = "075", t(coef(aci_075)))

# -----
# 076
# -----
# Fit curve
aci_076 <- licor_merged %>% filter(id == "76" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_076)

# Add to data frame
aci_coefs[76,] <- c(id = "076", t(coef(aci_076)))

# -----
# 077
# -----
# Fit curve
aci_077 <- licor_merged %>% filter(id == "77" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_077)

# Add to data frame
aci_coefs[77,] <- c(id = "077", t(coef(aci_077)))

# -----
# 078
# -----
# Fit curve
aci_078 <- licor_merged %>% filter(id == "78" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_078)

# Add to data frame
aci_coefs[78,] <- c(id = "078", t(coef(aci_078)))

# -----
# 079
# -----
# Fit curve
aci_079 <- licor_merged %>% filter(id == "79" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_079)

# Add to data frame
aci_coefs[79,] <- c(id = "079", t(coef(aci_079)))

# -----
# 080
# -----
# Fit curve
aci_080 <- licor_merged %>% filter(id == "80" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_080)

# Add to data frame
aci_coefs[80,] <- c(id = "080", t(coef(aci_080)))

# -----
# 081
# -----
# Fit curve
aci_081 <- licor_merged %>% filter(id == "81" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_081)

# Add to data frame
aci_coefs[81,] <- c(id = "081", t(coef(aci_081)))

# -----
# 082
# -----
# Fit curve
aci_082 <- licor_merged %>% filter(id == "82" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_082)

# Add to data frame
aci_coefs[82,] <- c(id = "082", t(coef(aci_082)))

# -----
# 083
# -----
# Fit curve
aci_083 <- licor_merged %>% filter(id == "83" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_083)

# Add to data frame
aci_coefs[83,] <- c(id = "083", t(coef(aci_083)))

# -----
# 084
# -----

# Add to data frame
aci_coefs[84,] <- c(id = "084", NA, NA, NA, NA)

# -----
# 085
# -----
# Fit curve
aci_085 <- licor_merged %>% filter(id == "85" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_085)

# Add to data frame
aci_coefs[85,] <- c(id = "085", t(coef(aci_085)))

# -----
# 086
# -----
# Fit curve
aci_086 <- licor_merged %>% filter(id == "86" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_086)

# Add to data frame
aci_coefs[86,] <- c(id = "086", t(coef(aci_086)))

# -----
# 087
# -----
# Fit curve
aci_087 <- licor_merged %>% filter(id == "87" & keep_row == "yes" & A > 0)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_087)

# Add to data frame
aci_coefs[87,] <- c(id = "087", t(coef(aci_087)))

# -----
# 088
# -----
# Fit curve
aci_088 <- licor_merged %>% filter(id == "88" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_088)

# Add to data frame
aci_coefs[88,] <- c(id = "088", t(coef(aci_088)))

# -----
# 089
# -----
# Fit curve
aci_089 <- licor_merged %>% filter(id == "89" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_089)

# Add to data frame
aci_coefs[89,] <- c(id = "089", t(coef(aci_089)))

# -----
# 090
# -----
# Fit curve
aci_090 <- licor_merged %>% filter(id == "90" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(aci_090)

# Add to data frame
aci_coefs[90,] <- c(id = "090", t(coef(aci_090)))

# -----
# 091
# -----
# Add to data frame
aci_coefs[91,] <- c(id = "091", NA, NA, NA, NA)

# -----
# 092
# -----
# Fit curve
aci_092 <- licor_merged %>% filter(id == "92" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_092)

# Add to data frame
aci_coefs[92,] <- c(id = "092", t(coef(aci_092)))

# -----
# 093
# -----
# Fit curve
aci_093 <- licor_merged %>% filter(id == "93" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_093)

# Add to data frame
aci_coefs[93,] <- c(id = "093", t(coef(aci_093)))

# -----
# 094
# -----
# Add to data frame
aci_coefs[94,] <- c(id = "094", NA, NA, NA, NA)

# -----
# 095
# -----
# Fit curve
aci_095 <- licor_merged %>% filter(id == "95" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_095)

# Add to data frame
aci_coefs[95,] <- c(id = "095", t(coef(aci_095)))

# -----
# 096
# -----
# Fit curve
aci_096 <- licor_merged %>% filter(id == "96" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_096)

# Add to data frame
aci_coefs[96,] <- c(id = "096", t(coef(aci_096)))

# -----
# 097
# -----
# Fit curve
aci_097 <- licor_merged %>% filter(id == "97" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_097)

# Add to data frame
aci_coefs[97,] <- c(id = "097", t(coef(aci_097)))

# -----
# 098
# -----
# Fit curve
aci_098 <- licor_merged %>% filter(id == "98" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_098)

# Add to data frame
aci_coefs[98,] <- c(id = "098", t(coef(aci_098)))

# -----
# 099
# -----
# Fit curve
aci_099 <- licor_merged %>% filter(id == "99" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_099)

# Add to data frame
aci_coefs[99,] <- c(id = "099", t(coef(aci_099)))

# -----
# 100
# -----
# Fit curve
aci_100 <- licor_merged %>% filter(id == "100" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_100)

# Add to data frame
aci_coefs[100,] <- c(id = "100", t(coef(aci_100)))

# -----
# 101
# -----
# Fit curve
aci_101 <- licor_merged %>% filter(id == "101" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_101)

# Add to data frame
aci_coefs[101,] <- c(id = "101", t(coef(aci_101)))

# -----
# 102
# -----
# Fit curve
aci_102 <- licor_merged %>% filter(id == "102" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_102)

# Add to data frame
aci_coefs[102,] <- c(id = "102", t(coef(aci_102)))

# -----
# 103
# -----
# Fit curve
aci_103 <- licor_merged %>% filter(id == "103" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_103)

# Add to data frame
aci_coefs[103,] <- c(id = "103", t(coef(aci_103)))

# -----
# 104
# -----
# Fit curve
aci_104 <- licor_merged %>% filter(id == "104" & keep_row == "yes" & A < 30)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_104)

# Add to data frame
aci_coefs[104,] <- c(id = "104", t(coef(aci_104)))

# -----
# 105
# -----
# Fit curve
aci_105 <- licor_merged %>% filter(id == "105" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_105)

# Add to data frame
aci_coefs[105,] <- c(id = "105", t(coef(aci_105)))

# -----
# 106
# -----
# Fit curve
aci_106 <- licor_merged %>% filter(id == "106" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_106)

# Add to data frame
aci_coefs[106,] <- c(id = "106", t(coef(aci_106)))

# -----
# 107
# -----
# Fit curve
aci_107 <- licor_merged %>% filter(id == "107" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_107)

# Add to data frame
aci_coefs[107,] <- c(id = "107", t(coef(aci_107)))

# -----
# 108
# -----
# Fit curve
aci_108 <- licor_merged %>% filter(id == "108" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_108)

# Add to data frame
aci_coefs[108,] <- c(id = "108", t(coef(aci_108)))

# -----
# 109
# -----
# Fit curve
aci_109 <- licor_merged %>% filter(id == "109" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_109)

# Add to data frame
aci_coefs[109,] <- c(id = "109", t(coef(aci_109)))

# -----
# 110
# -----
# Fit curve
aci_110 <- licor_merged %>% filter(id == "110" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_110)

# Add to data frame
aci_coefs[110,] <- c(id = "110", t(coef(aci_110)))

# -----
# 111
# -----
# Fit curve
aci_111 <- licor_merged %>% filter(id == "111" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_111)

# Add to data frame
aci_coefs[111,] <- c(id = "111", t(coef(aci_111)))

# -----
# 112
# -----
# Fit curve
aci_112 <- licor_merged %>% filter(id == "112" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_112)

# Add to data frame
aci_coefs[112,] <- c(id = "112", t(coef(aci_112)))

# -----
# 113
# -----
# Fit curve
aci_113 <- licor_merged %>% filter(id == "113" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_113)

# Add to data frame
aci_coefs[113,] <- c(id = "113", t(coef(aci_113)))

# -----
# 114
# -----
# Fit curve
aci_114 <- licor_merged %>% filter(id == "114" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_114)

# Add to data frame
aci_coefs[114,] <- c(id = "114", t(coef(aci_114)))

# -----
# 115
# -----
# Fit curve
aci_115 <- licor_merged %>% filter(id == "115_b" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_115)

# Add to data frame
aci_coefs[115,] <- c(id = "115", t(coef(aci_115)))

# -----
# 116
# -----
# Fit curve
aci_116 <- licor_merged %>% filter(id == "116" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_116)

# Add to data frame
aci_coefs[116,] <- c(id = "116", t(coef(aci_116)))

# -----
# 117
# -----
# Fit curve
aci_117 <- licor_merged %>% filter(id == "117" & keep_row == "yes" & A < 18)  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_117)

# Add to data frame
aci_coefs[117,] <- c(id = "117", t(coef(aci_117)))

# -----
# 118
# -----
# Fit curve
aci_118 <- licor_merged %>% filter(id == "118" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_118)

# Add to data frame
aci_coefs[118,] <- c(id = "118", t(coef(aci_118)))

# -----
# 119
# -----
# Fit curve
aci_119 <- licor_merged %>% filter(id == "119" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_119)

# Add to data frame
aci_coefs[119,] <- c(id = "119", t(coef(aci_119)))

# -----
# 120
# -----
# Fit curve
aci_120 <- licor_merged %>% filter(id == "120" & keep_row == "yes")  %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin"),
         Tcorrect = FALSE, fitTPU = TRUE)
plot(aci_120)

# Add to data frame
aci_coefs[120,] <- c(id = "120", t(coef(aci_120)))


#####################################################################
# Snapshot physiology measurements
#####################################################################
phys_data <- licor_merged %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  filter(id != "65" & id != "115") %>%
  mutate(id = ifelse(id == "115_b", 
                     "115",
                     ifelse(id == "65_b", "65", id)),
         id = str_pad(id, 3, pad = "0")) %>%
  dplyr::select(id, machine, date_only, A, Ci, Ca, gsw, Tleaf) %>%
  full_join(aci_coefs, by = "id") %>%
  filter(!is.na(id)) %>%
  mutate(Ci_Ca = Ci / Ca,
         iwue = A / gsw,
         l = stomatal_limitation(A_net = A, 
                                 Vcmax = as.numeric(Vcmax), 
                                 leaf.temp = Tleaf,
                                 temp = "C")$l,
         ind = id)

head(phys_data)

#####################################################################
# Merge physiology measurements with leaf area
#####################################################################
# Merge phys data with leaf area
compiled_df <- phys_data %>%
  full_join(tla_full, by = "ind") %>%
  dplyr::select(id = id.y, inoc:p_trt, ind, block, machine:l,
                tla_cm2:pfraction_area_cm2) %>%
  arrange(ind)

## Write compiled dataset
# write.csv(compiled_df, "../data_sheets/NxPxI_compiled_data.csv", 
#           row.names = F)
