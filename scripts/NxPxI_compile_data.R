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
  separate(id, sep = "_", 
           into = c("inoc", "n_trt", "p_trt", "block", "ind"), 
           remove = F) %>%
  mutate(n_trt = gsub("N", "", n_trt),
         p_trt = gsub("P", "", p_trt))
head(tla)

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
                         low.size = 0.1,
                         save.image = T)

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
  
write.csv(tla, "../data_sheets/NxPxI_tla_full.csv", row.names = F)

hist(tla$tla_cm2)

ggplot(data = tla, 
       aes(x = p_trt, y = tla_cm2, fill = n_trt)) +
  geom_jitter(size = 3, shape = 21) +
  geom_smooth(method = 'lm', aes(color = n_trt)) +
  facet_grid(~inoc) +
  theme_bw(base_size = 20)

library(lme4)
library(car)
library(emmeans)

tla$tla_cm2 <- as.numeric(tla$tla_cm2)
tla$n_trt <- factor(tla$n_trt, levels = c("0", "210"))
tla$p_trt <- as.numeric(tla$p_trt)


tla_lmer <- lmer(sqrt(tla_cm2) ~ inoc * p_trt * n_trt + (1 | block),
            data = tla)

shapiro.test(residuals(tla_lmer))

Anova(tla_lmer)

test(emtrends(tla_lmer, ~1, "p_trt"))
emmeans(tla_lmer, pairwise~n_trt)

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
  mutate(id = ifelse(id == "3", id == ""),
         Qin_cuvette = 2000) %>%
  dplyr::select(obs, time, elapsed, date, date_only, hhmmss:Qin,
                Qin_cuvette, Qabs:SS_r) %>%
  arrange(machine, date, id)

#####################################################################
# Fit LI-COR files
#####################################################################

# -----
# 001 
# -----

# Fit curve
aci_001 <- 









