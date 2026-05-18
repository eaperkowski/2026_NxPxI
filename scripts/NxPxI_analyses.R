# NxPxI_analyses.R: Script investigates treatment effects
# on leaf nutrient content, leaf physiology, biomass production,
# and biomass partitioning

######################################################################
# Libraries, data files, custom functions
#####################################################################
# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggpubr)
library(multcomp)
library(MuMIn)

# Load compiled data sheet
df <- read.csv("../data_sheets/NxPxI_compiled_data.csv") %>%
  mutate(p_trt = as.numeric(p_trt),
         n_trt = factor(n_trt, levels = c("0", "210")),
         inoc.ntrt = factor(str_c(inoc, ".", n_trt),
                            levels = c("NI.0", "YI.0", "NI.210", "YI.210"))) %>%
  filter(id != "YI_0N_7.5P_B3_066" & id != "YI_0N_0P_B4_063")
head(df)

full.cols <- c("#ffdb9e", "#b2182b", "#92c5d3", "#1e4a8c")
int.cols <- c("#d87964", "#5887af")

###############
# Anet
###############
# Model
anet_model <- lmer(A ~ inoc * n_trt * p_trt + (1 | block) + (1 | machine),
                   data = subset(df, p_trt < 33))

# Check normality assumptions
plot(anet_model)
qqnorm(residuals(anet_model))
qqline(residuals(anet_model))
hist(residuals(anet_model))
shapiro.test(residuals(anet_model))
outlierTest(anet_model)

# Model output
summary(anet_model)
Anova(anet_model)

# Plot
anet_plot <- ggplot(data = subset(df, p_trt < 32), 
                    aes(x = p_trt, y = A, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  labs(x = "P fertilization (ppm)",
     y = expression(bolditalic("A")[bold("net")]*bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
     color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 24, 21, 24),
                                                 alpha = 1))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
anet_plot

###############
# gsw
###############
# Model
gsw_model <- lmer(gsw ~ inoc * n_trt * p_trt + (1 | block) + (1 | machine),
            data = subset(df, p_trt < 33))

# Check normality assumptions
plot(gsw_model)
qqnorm(residuals(gsw_model))
qqline(residuals(gsw_model))
hist(residuals(gsw_model))
shapiro.test(residuals(gsw_model))
outlierTest(gsw_model)

# Model output
summary(gsw_model)
Anova(gsw_model)

# Post-hoc comparisons: Inoculation-by-N fertilization interaction
cld(emmeans(gsw_model, pairwise~inoc*n_trt))
## Inoculated plants have greater stomtal conductance, but only under low N

# Post-hoc comparisons: N fertilization-by-P fertilization interaction
test(emtrends(gsw_model, ~n_trt, "p_trt"))
## P fertilization decreases stomatal conductance under 0 ppm N, does not
## change under 210 ppm N

# Plot regression line prep
gsw_regline <- emmeans(gsw_model, ~inoc*n_trt*p_trt, 
                         type = "response", 
                         at = list(p_trt = seq(0, 31, 1))) %>%
  data.frame() %>%
  mutate(inoc.ntrt = factor(str_c(inoc, ".", n_trt),
                            levels = c("NI.0", "YI.0", "NI.210", "YI.210")),
         linetype = ifelse(n_trt == "0", "dashed", "solid"))


# Plot
gsw_plot <- ggplot(data = subset(df, p_trt < 32), 
                    aes(x = p_trt, y = gsw, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4,  width = 0.5) +
  geom_ribbon(data = gsw_regline,
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = inoc.ntrt),
              alpha = 0.3) +
  geom_smooth(data = gsw_regline, 
              aes(y = emmean, color = inoc.ntrt, linetype = linetype), 
              method = "loess", size = 1.5) +
  scale_color_manual(values = full.cols,
                     labels = c("0 ppm N, uninoculated",
                                "0 ppm N, inoculated",
                                "210 ppm N, uninoculated",
                                "210 ppm N, inoculated")) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  labs(x = "P fertilization (ppm)",
       y = expression(bolditalic("g")[bold("sw")]*bold(" (mol m"^"-2"*" s"^"-1"*")")),
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
gsw_plot

###############
# Vcmax
###############
# Model
vcmax_model <- lmer(Vcmax ~ inoc * n_trt * p_trt + (1 | block),
                    data = subset(df, p_trt < 33))

# Check normality assumptions
plot(vcmax_model)
qqnorm(residuals(vcmax_model))
qqline(residuals(vcmax_model))
hist(residuals(vcmax_model))
shapiro.test(residuals(vcmax_model))
outlierTest(vcmax_model)

# Model output
summary(vcmax_model)
Anova(vcmax_model)

# Post-hoc comparisons: inoculation-by-N fertilization interaction
cld(emmeans(vcmax_model, pairwise~inoc*n_trt))
## Vcmax25 is greater in inoculated plants, but only under 0 ppm N

# Plot
vcmax_plot <- ggplot(data = subset(df, p_trt < 32), 
                    aes(x = p_trt, y = Vcmax, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  labs(x = "P fertilization (ppm)",
       y = expression(bolditalic("V")[bold("cmax25")]*bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
vcmax_plot

###############
# Jmax
###############
# Model
jmax_model <- lmer(sqrt(Jmax) ~ inoc * n_trt * p_trt + (1 | block),
              data = subset(df, p_trt < 33))

# Check normality assumptions
plot(jmax_model)
qqnorm(residuals(jmax_model))
qqline(residuals(jmax_model))
hist(residuals(jmax_model))
shapiro.test(residuals(jmax_model))
outlierTest(jmax_model)

# Model output
summary(jmax_model)
Anova(jmax_model)

# Post-hoc comparisons: P fertilization-by-inoculation-by-N fertilization
cld(emmeans(jmax_model, pairwise~inoc*n_trt))
## Jmax25 is greater in inoculated plants, but only under 0 ppm N

# Plot
jmax_plot <- ggplot(data = subset(df, p_trt < 32), 
                     aes(x = p_trt, y = Jmax, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  labs(x = "P fertilization (ppm)",
       y = expression(bolditalic("J")[bold("max25")]*bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
jmax_plot

###############
# Jmax : Vcmax
###############
# Model
jvmax_model <- lmer(log(Jmax_Vcmax) ~ inoc * n_trt * p_trt + (1 | block),
                    data = subset(df, p_trt < 33))

# Check normality assumptions
plot(jvmax_model)
qqnorm(residuals(jvmax_model))
qqline(residuals(jvmax_model))
hist(residuals(jvmax_model))
shapiro.test(residuals(jvmax_model))
outlierTest(jvmax_model)

# Model output
summary(jvmax_model)
Anova(jvmax_model)

# Post-hoc comparisons: N fertilization-by-P fertilization interaction
test(emtrends(jvmax_model, pairwise~n_trt, "p_trt"))
## P fertilization significantly increases Jmax:Vcmax, but only under 0 ppm N

test(emtrends(jvmax_model, pairwise~n_trt*inoc, "p_trt"))

# Plot regression line prep
jvmax_regline <- emmeans(jvmax_model, ~inoc*n_trt*p_trt, 
                         type = "response", 
                         at = list(p_trt = seq(0, 31, 1))) %>%
  data.frame() %>%
  mutate(inoc.ntrt = factor(str_c(inoc, ".", n_trt),
                             levels = c("NI.0", "YI.0", "NI.210", "YI.210")),
         linetype = ifelse(n_trt == "0", "solid", "dashed"))

# Plot
jvmax_plot <- ggplot(data = subset(df, p_trt < 32), 
                     aes(x = p_trt, y = Jmax_Vcmax, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  geom_ribbon(data = jvmax_regline,
              aes(y = response, ymin = lower.CL, ymax = upper.CL, fill = inoc.ntrt),
              alpha = 0.3) +
  geom_smooth(data = jvmax_regline, 
              aes(y = response, color = inoc.ntrt, linetype = linetype), 
              method = "loess", size = 1.5) +
  scale_color_manual(values = full.cols,
                     labels = c("0 ppm N, uninoculated",
                                "0 ppm N, inoculated",
                                "210 ppm N, uninoculated",
                                "210 ppm N, inoculated")) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "P fertilization (ppm)",
       y = expression(bolditalic("J")[bold("max25")]*bold(":")*bolditalic("V")[bold("cmax25")]*bold(" (unitless)")),
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
jvmax_plot

###############
# iWUE
###############
# Remove outliers
df$iwue[19] <- NA

# Model
iwue_model <- lmer(log(iwue) ~ inoc * n_trt * p_trt + (1 | block),
                    data = subset(df, p_trt < 33))

# Check normality assumptions
plot(iwue_model)
qqnorm(residuals(iwue_model))
qqline(residuals(iwue_model))
hist(residuals(iwue_model))
shapiro.test(residuals(iwue_model))
outlierTest(iwue_model)

# Model output
summary(iwue_model)
Anova(iwue_model)

# Pairwise comparisons: N fertilization-by-P fertilization interaction
test(emtrends(iwue_model, pairwise~n_trt, "p_trt"))
## P fertilization increases iWUE, but only under 0 ppm N

# Pairwise comparisons: N fertilization effect
emmeans(iwue_model, pairwise~n_trt, type = "response")
## iWUE decreases with increasing N fertilization

# Plot regression line prep
iwue_regline <- emmeans(iwue_model, ~inoc*n_trt*p_trt, 
                         type = "response", 
                         at = list(p_trt = seq(0, 31, 1))) %>%
  data.frame() %>%
  mutate(inoc.ntrt = factor(str_c(inoc, ".", n_trt),
                            levels = c("NI.0", "YI.0", "NI.210", "YI.210")),
         linetype = ifelse(n_trt == "0", "solid", "dashed"))

# Plot
iwue_plot <- ggplot(data = subset(df, p_trt < 32), 
                     aes(x = p_trt, y = iwue, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  geom_ribbon(data = iwue_regline,
              aes(y = response, ymin = lower.CL, ymax = upper.CL, fill = inoc.ntrt),
              alpha = 0.3) +
  geom_smooth(data = iwue_regline, 
              aes(y = response, color = inoc.ntrt, linetype = linetype), 
              method = "loess", size = 1.5) +
  scale_color_manual(values = full.cols,
                     labels = c("0 ppm N, uninoculated",
                                "0 ppm N, inoculated",
                                "210 ppm N, uninoculated",
                                "210 ppm N, inoculated")) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "P fertilization (ppm)",
       y = expression(bolditalic("i")*bold("WUE ("*mu*"mol mol"^"-1"*")")),
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
iwue_plot

###############
# Stom-lim
###############
# Remove outliers
df$l[20] <- NA

# Model
l_model <- lmer(l ~ inoc * n_trt * p_trt + (1 | block),
                   data = subset(df, p_trt < 33))

# Check normality assumptions
plot(l_model)
qqnorm(residuals(l_model))
qqline(residuals(l_model))
hist(residuals(l_model))
shapiro.test(residuals(l_model))
outlierTest(l_model)

# Model output
summary(l_model)
Anova(l_model)

# Post-hoc comparisons: marginal inoculation-by-P fertilization interaction
test(emtrends(l_model, pairwise~inoc, "p_trt"))
## Inoculated plants experience a marginal increase in stomatal limitation with
## increasing P fertilization, no response in uninoculated plants

# Post-hoc comparisons: N fertilization effect
emmeans(l_model, pairwise~n_trt)
## N fertilization decreases stomatal limitation

# Post-hoc comparisons: Inoculation effect
emmeans(l_model, pairwise~inoc)
## Inoculation decreases stomatal limitation

test(emtrends(l_model, pairwise~inoc*n_trt, "p_trt"))

# Plot regression line prep
l_regline <- emmeans(l_model, ~inoc*n_trt*p_trt, 
                        type = "response", 
                        at = list(p_trt = seq(0, 31, 1))) %>%
  data.frame() %>%
  mutate(inoc.ntrt = factor(str_c(inoc, ".", n_trt),
                            levels = c("NI.0", "YI.0", "NI.210", "YI.210")),
         linetype = ifelse(n_trt == "0" & inoc == "YI", "solid", "dashed"))

# Plot
l_plot <- ggplot(data = subset(df, p_trt < 32), 
                    aes(x = p_trt, y = l, fill = inoc.ntrt)) +
  geom_jitter(aes(shape = inoc), size = 4, width = 0.5) +
  geom_ribbon(data = l_regline,
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = inoc.ntrt),
              alpha = 0.3) +
  geom_smooth(data = l_regline, 
              aes(y = emmean, color = inoc.ntrt, linetype = linetype), 
              method = "loess", size = 1.5) +
  scale_color_manual(values = full.cols,
                     labels = c("0 ppm N, uninoculated",
                                "0 ppm N, inoculated",
                                "210 ppm N, uninoculated",
                                "210 ppm N, inoculated")) +
  scale_fill_manual(values = full.cols,
                    labels = c("0 ppm N, uninoculated",
                               "0 ppm N, inoculated",
                               "210 ppm N, uninoculated",
                               "210 ppm N, inoculated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "P fertilization (ppm)",
       y = "Stomatal limitation (unitless)",
       color = "Treatment", fill = "Treatment") +
  guides(shape = "none",
         linetype = "none",
         fill = guide_legend(override.aes = list(shape = c(21, 21, 24, 24)))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(size = 1.25),
        legend.text.align = 0)
l_plot


png("../drafts/figs/NxPxI_snapshot_phys.png", 
    height = 10, width = 14, units = "in", res = 600)
ggarrange(anet_plot, gsw_plot, iwue_plot, l_plot,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "right",
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)"))
dev.off()

png("../drafts/figs/NxPxI_photoCapacity.png", 
    height = 10, width = 14, units = "in", res = 600)
ggarrange(vcmax_plot, jmax_plot, jvmax_plot,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "right",
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)"))
dev.off()
