
# Code author : Ihsan Fadilah
# Project     : Levels of 5,6-orthoquinone derivative of primaquine in urine
#               correlate with CYP2D6 genotype-predicted metabolic activity
#               score
# About       : Graphical representation of data,
#               graphical summaries of data,
#               modelling (graphical depictions, numerical estimates)

# Setup -------------------------------------------------------------------
library(tidyverse)
library(haven)
library(here)
library(rms)
library(ggbeeswarm)
library(lme4)
library(pracma)
library(zoo)

# Plot customisation
# Set plots to some format
primary <- c('#272636', '#E3932B', '#FECD5B', '#536E85', '#FFFFFF')
special <- c('#AF9259', '#918580', '#8F2736', '#C82F46')
theme_set(theme_bw())

## Further refinement
theme_update(
  text = element_text(size = 10, family = "Foundry Sterling"), # Font
  plot.title = element_text(hjust = 0),      # Centre-align title
  plot.subtitle = element_text(hjust = 0),   # Centre-align subtitle
  legend.title = element_text(colour = 'black',
                              size = 8, 
                              face = 'bold'),# Legend title
  legend.position = 'right',                 # Move legend
  legend.background = element_blank(),       # Remove legend background
  legend.box.background = element_blank(),   # Remove lengend-box background
  legend.spacing.y = unit(0.01, 'mm'),       # Make legend closer
  legend.key.height = unit(0.4, "cm"),       # Make legend closer
  # panel.grid.minor = element_blank(),      # Remove minor lines
  panel.grid.minor.x = element_blank(),      # Remove minor lines on the x axis
  axis.title.x = element_text(hjust = 1),    # Move title for x-axis
  axis.ticks = element_blank(),              # Remove axis ticks
  aspect.ratio = 1,
  axis.title.y = element_text(hjust = 0.5)   # Move title for y-axis
)

# Data --------------------------------------------------------------------
prima0to2 <- read_rds(file = here('data', 'prima0to2.rds')) |> 
  group_by(timepoint_cont) |> 
  mutate(methb_max = max(methb, na.rm = T),
         methb_q50 = median(methb, na.rm = T),
         methb_mean = mean(methb, na.rm = T)) |> 
  ungroup() |> 
  group_by(day, patid) |> 
  mutate(methb_max_dayid = max(methb, na.rm = T),
         methb_max_dayid = if_else(methb_max_dayid == -Inf, NA, methb_max_dayid),
         methb_q50_dayid = median(methb, na.rm = T)) |> 
  ungroup() |> 
  mutate(cyp_cat1 = factor(cyp_cat1, levels = c('Intermediate', 'Normal', 'Ultrarapid')))

# Linear interpolation for AUC calculation
prima0to2_imp <- prima0to2 |> 
  mutate(imp_methb = if_else(is.na(methb), 'Imputed', 'Observed')) |> 
  group_by(day, patid) |> 
  mutate(methb = na.approx(methb, timepoint, na.rm = FALSE),
         methb = na.locf(methb, na.rm = FALSE),
         methb = na.locf(methb, fromLast = TRUE, na.rm = FALSE)) |> 
  ungroup()

prima0to2_imp |> 
  drop_na(methb) %>% # Do better here!
  ggplot() +
  geom_line(aes(x = timepoint,
                y = log2(methb),
                group = patid,
                colour = imp_methb),
            alpha = 0.7, linewidth = 0.8) +
  facet_grid(~day) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       colour = '',
       caption = 'Daily pre-dose measurements at timepoint 0') +
  scale_colour_manual(
    values = c('Imputed' = '#C82F46',
               'Observed' = '#536E85')
  )
  
n_patid <- prima0to2$patid |> unique() |> length()
patid <- prima0to2$patid |> unique()
timepoint <- prima0to2$timepoint |> unique()

# Day 0
auc_day0 <- rep(NA, n_patid)
for (i in 1:n_patid) {
  id <- patid[i]
  x <- filter(prima0to2_imp, day == 'Day 0', patid == id) |> select(methb) |> as.vector()
  methb <- x$methb
  auc <- trapz(timepoint, methb)
  auc_day0[i] <- auc
}

day0_auc <- tibble(patid = patid, methb_auc_dayid = auc_day0, day = 'Day 0')

# Day 1
auc_day1 <- rep(NA, n_patid)
for (i in 1:n_patid) {
  id <- patid[i]
  x <- filter(prima0to2_imp, day == 'Day 1', patid == id) |> select(methb) |> as.vector()
  methb <- x$methb
  auc <- trapz(timepoint, methb)
  auc_day1[i] <- auc
}

day1_auc <- tibble(patid = patid, methb_auc_dayid = auc_day1, day = 'Day 1')

# Day 2
auc_day2 <- rep(NA, n_patid)
for (i in 1:n_patid) {
  id <- patid[i]
  x <- filter(prima0to2_imp, day == 'Day 2', patid == id) |> select(methb) |> as.vector()
  methb <- x$methb
  auc <- trapz(timepoint, methb)
  auc_day2[i] <- auc
}

day2_auc <- tibble(patid = patid, methb_auc_dayid = auc_day2, day = 'Day 2')

methb_auc <- bind_rows(day0_auc, day1_auc, day2_auc)

prima0to2 <- left_join(prima0to2, methb_auc, by = c('day', 'patid'))

# At timepoint == 0 only
prima_base <- prima0to2 |>
  filter(timepoint == 0) |> 
  select(-c(clinid))

# Graphical, non-summary --------------------------------------------------
# Primaquine daily dose
prima_base |> 
  ggplot() +
  geom_histogram(aes(pqmgkgday, fill = group), bins = 23) +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  scale_x_continuous(limits = c(0, 1.5),
                     breaks = seq(0, 1.5, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(aspect.ratio = 0.75) +
  labs(x = '\nPrimaquine daily dose (mg/kg)',
       y = 'Count\n',
       fill = 'Treatment allocation') +
  scale_fill_manual(
    values = c('Intervention' = '#C82F46',
               'Control' = '#536E85')
  )

# Methaemoglobin over timepoints, by day
prima0to2 %>%
  drop_na(methb) %>% # Do better here!
  ggplot() +
  geom_line(aes(x = timepoint,
                y = log2(methb),
                group = patid),
            alpha = 0.3, linewidth = 0.8) +
  facet_grid(~day) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoint 0')

# Methaemoglobin over timepoints, connected
prima0to2 %>%
  drop_na(methb) %>% # Do better here!
  ggplot() +
  annotate("rect", xmin = c(8, 17), xmax = c(9, 18), ymin = -Inf, ymax = Inf,
           fill = '#918580', alpha = 0.3) +
  geom_vline(xintercept = c(0, 9, 18),
             linetype = 'dotted') +
  geom_line(aes(x = timepoint_cont,
                y = log2(methb),
                group = patid),
            alpha = 0.3, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 26),
                     breaks = seq(0, 26, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 0.4) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoints 0, 9, and 18')

# Methaemoglobin over timepoints, connected (by CYP2D6)
prima0to2 %>%
  drop_na(methb) %>% # Do better here!
  ggplot() +
  annotate("rect", xmin = c(8, 17), xmax = c(9, 18), ymin = -Inf, ymax = Inf,
           fill = '#918580', alpha = 0.3) +
  geom_vline(xintercept = c(0, 9, 18),
             linetype = 'dotted') +
  geom_line(aes(x = timepoint_cont,
                y = log2(methb),
                group = patid,
                colour = cyp_cat1),
            alpha = 0.7, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 26),
                     breaks = seq(0, 26, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  theme(aspect.ratio = 0.4,
        legend.key.width = unit(0.3, 'cm')) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       colour = 'CYP2D6 activity',
       caption = 'Daily pre-dose measurements at timepoints 0, 9, and 18') +
  scale_colour_manual(
    values = c('Intermediate' = '#C82F46',
               'Normal' = '#536E85')
  )

# Last timepoint but base timepoint?
# Important for overall time-series
# Likely not completely reliable after scrutinising the raw dataset

# Graphical, summary ------------------------------------------------------
# Lines, by day (median)
prima0to2_q50 <- prima0to2 |>
  select(patid, day, timepoint, methb_q50) |> 
  mutate(type = 'Median')

prima0to2_max <- prima0to2 |>
  select(patid, day, timepoint, methb_max) |> 
  mutate(type = 'Maximum')

prima0to2_mean <- prima0to2 |>
  select(patid, day, timepoint, methb_mean) |> 
  mutate(type = 'Mean')

prima0to2 %>%
  ggplot() +
  geom_line(aes(x = timepoint, y = log2(methb_q50)),
            alpha = 0.3, linewidth = 0.8) +
  geom_line(aes(x = timepoint, y = log2(methb_mean)),
            alpha = 0.3, linewidth = 0.8) +
  geom_line(aes(x = timepoint, y = log2(methb_max)),
            alpha = 0.3, linewidth = 0.8) +
  geom_point(aes(x = timepoint, y = log2(methb_q50), colour = type),
             size = 0.9, data = prima0to2_q50) +
  geom_point(aes(x = timepoint, y = log2(methb_mean), colour = type),
             size = 0.9, data = prima0to2_mean) +
  geom_point(aes(x = timepoint, y = log2(methb_max), colour = type),
             size = 0.9, data = prima0to2_max) +
  facet_grid(~day) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.5, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoint 0') + 
  theme(legend.key.width = unit(0.1, 'cm')) +
  scale_colour_manual(
    name = "Summary metric",  # Legend title
    values = c('Median' = '#E3932B',
               'Mean' = '#536E85',
               'Maximum' = '#C82F46')
  )

# Boxplots, by day
prima0to2 |> 
  mutate(timepoint = factor(timepoint)) |> 
  ggplot() +
  geom_boxplot(aes(x = timepoint, y = log2(methb)),
               varwidth = TRUE, fill = 'transparent',
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.7, linewidth = 0.3) +
  facet_wrap(~day) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoint 0')

# Boxplots over timepoints, connected
prima0to2 |> 
  ggplot() +
  annotate("rect", xmin = c(8, 17), xmax = c(9, 18),
           ymin = -Inf, ymax = Inf,
           fill = '#918580', alpha = 0.3) +
  geom_vline(xintercept = c(0, 9, 18),
             linetype = 'dotted') +
  geom_boxplot(aes(x = timepoint_cont, y = log2(methb), group = timepoint_cont),
               varwidth = TRUE, fill = 'transparent',
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.55, linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 26.6),
                     breaks = seq(0, 26, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.25, 16)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  theme(panel.grid.major.x = element_blank(),
        aspect.ratio = 0.425) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoints 0, 9, and 18')

# Boxplot, by day
## Maximum
prima_base |> 
  ggplot() +
  geom_beeswarm(aes(y = log2(methb_max_dayid), x = day),
                groupOnX = F, alpha = 0.3, size = 1.25) +
  geom_boxplot(aes(y = log2(methb_max_dayid), x = day),
               varwidth = T, # !!!
               fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.5, linewidth = 0.5) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Maximum methaemoglobin (%)\n')

## AUC
upper_half <- 2^((log2(64) + log2(128))/2)
lower_half <- 2^((log2(2) + log2(4))/2)
prima_base |> 
  ggplot() +
  geom_beeswarm(aes(y = log2(methb_auc_dayid), x = day),
                groupOnX = F, alpha = 0.3, size = 1.25) +
  geom_boxplot(aes(y = log2(methb_auc_dayid), x = day),
               varwidth = T, # !!!
               fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.5, linewidth = 0.5) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Methaemoglobin AUC\n')

# Lines, by day and CYP2D6 group
## Maximum
prima_base |>
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_point(aes(y = log2(methb_max_dayid), x = day, colour = cyp_cat1),
             alpha = 1, size = 1.2) +
  geom_line(aes(y = log2(methb_max_dayid), x = day, group = patid, colour = cyp_cat1),
            alpha = 0.3, size = 0.7) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Maximum methaemoglobin (%)\n',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Intermediate' = '#C82F46',
               'Normal' = '#536E85')
  )

## AUC
prima_base |>
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_point(aes(y = log2(methb_auc_dayid), x = day, colour = cyp_cat1),
             alpha = 1, size = 1.2) +
  geom_line(aes(y = log2(methb_auc_dayid), x = day, group = patid, colour = cyp_cat1),
            alpha = 0.3, size = 0.7) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Methaemoglobin AUC\n',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Intermediate' = '#C82F46',
               'Normal' = '#536E85')
  )

# Maximum methaemoglobin, by day and some factor
## Primaquine dose
### Maximum
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  geom_point(aes(y = log2(methb_max_dayid), x = pqmgkgday),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = pqmgkgday, y = log2(methb_max_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0.9, 1.3),
                     breaks = seq(0.8, 1.4, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = 'Maximum methaemoglobin (%)\n',
       x = '\nPrimaquine daily dose (mg/kg)',
       colour = '') +
  theme(legend.key.width = unit(0.5, 'cm'),
        panel.spacing = unit(1.3, "lines"),) +
  facet_wrap(~day)

### AUC
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  geom_point(aes(y = log2(methb_auc_dayid), x = pqmgkgday),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = pqmgkgday, y = log2(methb_auc_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0.9, 1.3),
                     breaks = seq(0.8, 1.4, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = 'Methaemoglobin AUC\n',
       x = '\nPrimaquine daily dose (mg/kg)',
       colour = '') +
  theme(legend.key.width = unit(0.5, 'cm'),
        panel.spacing = unit(1.3, "lines")) +
  facet_wrap(~day)

## CYP2D6
### Maximum
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_jitter(aes(y = log2(methb_max_dayid), x = day, colour = cyp_cat1),
                  alpha = 0.3, size = 1.25, width = 0.05) +
  geom_boxplot(aes(y = log2(methb_max_dayid), x = day, colour = cyp_cat1),
               varwidth = F, # !!!
               fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.65, linewidth = 0.5) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Maximum methaemoglobin (%)\n',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Intermediate' = '#C82F46',
               'Normal' = '#536E85')
  )

### AUC
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_jitter(aes(y = log2(methb_auc_dayid), x = day, colour = cyp_cat1),
                alpha = 0.3, size = 1.25, width = 0.05) +
  geom_boxplot(aes(y = log2(methb_auc_dayid), x = day, colour = cyp_cat1),
               varwidth = F, # !!!
               fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.65, linewidth = 0.5) +
  theme(panel.grid.major.x = element_blank(),
        legend.key.width = unit(0.4, 'cm')) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '',
       y = 'Methaemoglobin AUC\n',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Intermediate' = '#C82F46',
               'Normal' = '#536E85')
  )

## 5,6-orthoquinone
### Maximum
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_point(aes(x = log2(orthoq), y = log2(methb_max_dayid)),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = log2(orthoq), y = log2(methb_max_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_x_continuous(limits = log2(c(64, 3000)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\n5,6-orthoquinone (ng/ml)',
       y = 'Maximum methaemoglobin (%)\n',
       colour = 'CYP2D6 activity') +
  theme(legend.key.width = unit(0.5, 'cm'),
        legend.position = 'top',
        panel.spacing = unit(1.1, "lines")) +
  facet_wrap(~day)

### AUC
prima_base |> 
  ggplot() +
  geom_point(aes(x = log2(orthoq), y = log2(methb_auc_dayid)),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = log2(orthoq), y = log2(methb_auc_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_x_continuous(limits = log2(c(64, 3000)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\n5,6-orthoquinone (ng/ml)',
       y = 'Methaemoglobin AUC\n',
       colour = '') +
  theme(legend.key.width = unit(0.5, 'cm'),
        panel.spacing = unit(1.1, "lines")) +
  facet_wrap(~day)

## G6PD
### Maximum
prima_base |> 
  ggplot() +
  geom_point(aes(x = g6pd1, y = log2(methb_max_dayid)),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = g6pd1, y = log2(methb_max_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_x_continuous(limits = c(6, 12.5),
                     breaks = seq(0, 20, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nG6PD activity (IU/gHb)',
       y = 'Maximum methaemoglobin (%)\n',
       colour = '') +
  theme(legend.key.width = unit(0.5, 'cm')) +
  facet_wrap(~day)

### AUC
prima_base |> 
  ggplot() +
  geom_point(aes(x = g6pd1, y = log2(methb_auc_dayid)),
             size = 1.5, alpha = 1) +
  geom_smooth(aes(x = g6pd1, y = log2(methb_auc_dayid)),
              method = "lm", se = F, linewidth = 0.8, colour = '#C82F46') +
  scale_x_continuous(limits = c(6, 12.5),
                     breaks = seq(0, 20, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nG6PD activity (IU/gHb)',
       y = 'Methaemoglobin AUC\n',
       colour = '') +
  theme(legend.key.width = unit(0.5, 'cm')) +
  facet_wrap(~day)

## Methaemoglobin with daily dose, by CYP2D6
### Maximum
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  geom_point(aes(y = log2(methb_max_dayid), x = pqmgkgday, colour = cyp_cat1),
             size = 1.5, alpha = 0.35) +
  geom_smooth(aes(x = pqmgkgday, y = log2(methb_max_dayid), colour = cyp_cat1),
              method = "lm", se = F, linewidth = 0.8) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0.9, 1.3),
                     breaks = seq(0.8, 1.4, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = 'Maximum methaemoglobin (%)\n',
       x = '\nPrimaquine daily dose (mg/kg)',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Normal' = '#536E85',
               'Intermediate' = '#C82F46')
  ) +
  theme(legend.key.width = unit(0.5, 'cm'),
        legend.position = 'top',
        panel.spacing = unit(1.3, "lines"),
        plot.margin = margin(10, 10, 10, 10)) +
  facet_wrap(~day)

### AUC
prima_base |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  geom_point(aes(y = log2(methb_auc_dayid), x = pqmgkgday, colour = cyp_cat1),
             size = 1.5, alpha = 0.35) +
  geom_smooth(aes(x = pqmgkgday, y = log2(methb_auc_dayid), colour = cyp_cat1),
              method = "lm", se = F, linewidth = 0.8) +
  scale_y_continuous(limits = log2(c(lower_half, upper_half)),
                     breaks = log2(c(2, 4, 8, 16, 32, 64, 128)),
                     labels = c(2, 4, 8, 16, 32, 64, 128),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0.9, 1.3),
                     breaks = seq(0.8, 1.4, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = 'Methaemoglobin AUC\n',
       x = '\nPrimaquine daily dose (mg/kg)',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Normal' = '#536E85',
               'Intermediate' = '#C82F46')
  ) +
  theme(legend.key.width = unit(0.5, 'cm'),
        legend.position = 'top',
        panel.spacing = unit(1.3, "lines"),
        plot.margin = margin(10, 10, 10, 10)) +
  facet_wrap(~day)

## 5,6-orthoquinone with ...
cyp_score <- prima_base$cyp_score1 |> unique()
cyp_score <- cyp_score[2:length(cyp_score)]

### Maximum
prima_base |> 
  drop_na(cyp_cat1) |> 
  filter(day == 'Day 2') |> 
  ggplot() +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  geom_point(aes(y = log2(orthoq), x = pqmgkgday, colour = cyp_cat1),
             size = 1.5, alpha = 0.35) +
  geom_smooth(aes(x = pqmgkgday, y = log2(orthoq), colour = cyp_cat1),
              method = "lm", se = F, linewidth = 0.8) +
  scale_y_continuous(limits = log2(c(64, 3000)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0.9, 1.3),
                     breaks = seq(0.8, 1.4, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = '5,6-orthoquinone (ng/ml)\n',
       x = '\nPrimaquine daily dose (mg/kg)',
       colour = 'CYP2D6 activity') +
  scale_colour_manual(
    values = c('Normal' = '#E3932B',
               'Intermediate' = '#C82F46',
               'Ultrarapid' = '#536E85')
  ) +
  theme(legend.key.width = unit(0.5, 'cm'),
        # legend.position = 'top',
        panel.spacing = unit(1.3, "lines"),
        plot.margin = margin(10, 10, 10, 10))

# ### Misc.
# prima_base |> 
#   drop_na(cyp_cat1) |> 
#   ggplot() +
#   geom_vline(xintercept = 1, linetype = 'dotted') +
#   geom_point(aes(y = log2(orthoq), x = cyp_score1),
#              size = 1.5, alpha = 0.35) +
#   geom_smooth(aes(x = cyp_score1, y = log2(orthoq)),
#               method = "lm", se = F, linewidth = 0.8) +
#   scale_y_continuous(limits = log2(c(64, 3000)),
#                      breaks = log2(2^seq(1, 15, by = 1)),
#                      labels = 2^seq(1, 15, by = 1),
#                      expand = expansion(mult = c(0, 0))) +
#   scale_x_continuous(breaks = seq(0.25, 2, by = 0.25)) +
#   labs(y = '5,6-orthoquinone (ng/ml)\n',
#        x = '\nCYP2D6 activity score')

# Models ------------------------------------------------------------------
prima_base_d2 <- prima_base |>
  filter(day == 'Day 2') |>
  drop_na(methb_max_dayid) |> 
  mutate(log2_orthoq = log2(orthoq),
         log2_methb = log2(methb_max_dayid),
         log2_methb_auc = log2(methb_auc_dayid),
         methb = methb_max_dayid,
         methb_auc = methb_auc_dayid,
         orthoq1000 = orthoq/1000)
dd <- datadist(prima_base_d2)
options(datadist = 'dd')

# # AUC
# q2a_auc <- ols(log2_methb_auc ~ cyp_cat1, data = prima_base_d2); anova(q2a_auc)
# ttest_auc <- t.test(log2_methb_auc ~ cyp_cat1, data = prima_base_d2, var.equal = T); ttest_auc
# ttest_mean_auc <- ttest_auc$estimate
# ttest_df <- tibble(
#   cyp_cat1 = c('Intermediate', 'Normal'),
#   methb = ttest_mean_auc,
#   methb2 = 2^methb
# )
# q2a2_auc <- ols(log2_methb_auc ~ cyp_score1, data = prima_base_d2); anova(q2a2_auc)
# q2a3_auc <- ols(log2_methb_auc ~ cyp_score2, data = prima_base_d2); anova(q2a3_auc)

# Maximum
# Question 2a
# On day 2, methaemoglobin with CYP2D6 activity
q2a <- lm(log2_methb ~ cyp_cat1, data = prima_base_d2); anova(q2a)
performance::check_model(q2a)
performance::check_heteroscedasticity(q2a)
residualsq2a <- residuals(q2a)
qqnorm(residualsq2a); qqline(residualsq2a, col = '#C82F46'); shapiro.test(residualsq2a)

ttest <- t.test(log2_methb ~ cyp_cat1, data = prima_base_d2, var.equal = T); ttest
ttest_mean <- ttest$estimate
ttest_df <- tibble(
  cyp_cat1 = c('Intermediate', 'Normal'),
  methb = ttest_mean,
  methb2 = 2^methb
)

q2a2_lm <- lm(log2_methb ~ cyp_score1, data = prima_base_d2); anova(q2a2_lm)
q2a2 <- ols(log2_methb ~ cyp_score1, data = prima_base_d2); anova(q2a2)
performance::check_model(q2a2_lm)
performance::check_heteroscedasticity(q2a2_lm)
residualsq2a2 <- residuals(q2a2_lm)
qqnorm(residualsq2a2); qqline(residualsq2a2, col = '#C82F46'); shapiro.test(residualsq2a2)

q2a3_lm <- lm(log2_methb ~ cyp_score2, data = prima_base_d2); anova(q2a3_lm)
q2a3 <- ols(log2_methb ~ cyp_score2, data = prima_base_d2); anova(q2a3)
performance::check_model(q2a3_lm)
performance::check_heteroscedasticity(q2a3_lm)
residualsq2a3 <- residuals(q2a3_lm)
qqnorm(residualsq2a3); qqline(residualsq2a3, col = '#C82F46'); shapiro.test(residualsq2a3)

AIC(q2a); AIC(q2a2); AIC(q2a3)
BIC(q2a); BIC(q2a2); BIC(q2a3)

q2a2_pred <- Predict(q2a2,
                     cyp_score1 = seq(0, 2.5, by = 0.01),
                     ref.zero = F) |>
  as_tibble() |> 
  mutate(yhat = 2^yhat,
         lower = 2^lower,
         upper = 2^upper); q2a2_pred
q2a3_pred <- Predict(q2a3,
                     cyp_score2 = seq(0, 2.5, by = 0.01),
                     ref.zero = F) |>
  as_tibble() |> 
  mutate(yhat = 2^yhat,
         lower = 2^lower,
         upper = 2^upper); q2a2_pred

prima_base_d2 |> 
  ggplot() +
  geom_jitter(aes(cyp_cat1, methb), width = 0.05,
             colour = 'black', alpha = 0.3, size = 2) +
  geom_point(aes(x = cyp_cat1, y = 2^methb),
             data = ttest_df, shape = 5, size = 2) +
  geom_boxplot(aes(cyp_cat1, methb),
               varwidth = TRUE, fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.7, linewidth = 0.3, alpha = 0.3) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity',
       y = 'Day-2 methaemoglobin (%)\n')
prima_base_d2 |> 
  ggplot() +
  geom_point(aes(cyp_score1, methb),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score1, ymin = lower, ymax = upper),
              data = q2a2_pred, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score1, y = yhat),
            data = q2a2_pred, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 2, by = 0.25),
                     limits = c(0.2, 2.05),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = 'Day-2 methaemoglobin (%)\n')
prima_base_d2 |> 
  ggplot() +
  geom_point(aes(cyp_score2, methb),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score2, ymin = lower, ymax = upper),
              data = q2a3_pred, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score2, y = yhat),
            data = q2a3_pred, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 2.5, by = 0.25),
                     limits = c(0.2, 2.3),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = 'Day-2 methaemoglobin (%)\n')

# Question 2b
# On day 2, methaemoglobin with 5,6-orthoquinone
q2b0 <- ols(methb ~ orthoq, data = prima_base_d2); anova(q2b0)
q2b0_lm <- lm(methb ~ orthoq, data = prima_base_d2); anova(q2b0_lm)
performance::check_model(q2b0_lm)
performance::check_heteroscedasticity(q2b0_lm)
performance::check_normality(q2b0_lm)
residualsq2b0_lm <- residuals(q2b0_lm)
qqnorm(residualsq2b0_lm); qqline(residualsq2b0_lm, col = "red"); shapiro.test(residualsq2b0_lm)

q2b <- ols(log2_methb ~ log2_orthoq, data = prima_base_d2); anova(q2b)
q2b_lm <- lm(log2_methb ~ log2_orthoq, data = prima_base_d2); anova(q2b_lm)
performance::check_model(q2b_lm)
performance::check_heteroscedasticity(q2b_lm)
performance::check_normality(q2b_lm)
residualsq2b_lm <- residuals(q2b_lm)
qqnorm(residualsq2b_lm); qqline(residualsq2b_lm, col = "red"); shapiro.test(residualsq2b_lm)

2^q2b$coefficients
2^confint(q2b)

q2b_ari <- ols(log2_methb ~ orthoq1000, data = prima_base_d2); anova(q2b_ari)
q2b_ari_lm <- lm(log2_methb ~ orthoq1000, data = prima_base_d2); anova(q2b_ari_lm)
performance::check_model(q2b_ari_lm)
performance::check_heteroscedasticity(q2b_ari_lm)
performance::check_normality(q2b_ari_lm)
residualsq2b_ari_lm <- residuals(q2b_ari_lm)
qqnorm(residualsq2b_ari_lm); qqline(residualsq2b_ari_lm, col = "red"); shapiro.test(residualsq2b_ari_lm)

2^q2b_ari$coefficients
2^confint(q2b_ari)

AIC(q2b0_lm, q2b_lm, q2b_ari_lm)

# 2^q2b_ari2$coefficients
# 2^confint(q2b_ari2)

q2b_pred <- Predict(q2b,
                    log2_orthoq = log2(seq(50, 2500, by = 1)),
                    ref.zero = F) |>
  as_tibble(); q2b_pred
q2b_ari_pred <- Predict(q2b_ari,
                        orthoq1000 = seq(50, 2500, by = 1)/1000,
                        ref.zero = F) |>
  as_tibble() |> 
  mutate(yhat = 2^yhat,
         lower = 2^lower,
         upper = 2^upper); q2b_ari_pred
# q2b_ari_pred2 <- Predict(q2b_ari2,
#                         orthoq1000 = seq(50, 2500, by = 1)/1000,
#                         ref.zero = F) |>
#   as_tibble() |> 
#   mutate(yhat = 2^yhat,
#          lower = 2^lower,
#          upper = 2^upper); q2b_ari_pred2

# q2b_pred2 <- Predict(q2b,
#                      conf.type = 'individual',
#                      orthoq = seq(50, 2500, by = 1),
#                      ref.zero = F) |>
#   as_tibble() |> 
#   mutate(yhat = 2^yhat,
#          lower = 2^lower,
#          upper = 2^upper)
# q2b_pred2

(q2b_ori_plot <- q2b_pred |> 
  mutate(yhat = 2^yhat,
         lower = 2^lower,
         upper = 2^upper,
         log2_orthoq = 2^log2_orthoq) |> 
  ggplot() +
  geom_point(data = prima_base_d2,
             colour = 'black', alpha = 0.3,
             aes(orthoq, y = methb_max_dayid)) +
  # geom_ribbon(aes(x = orthoq, ymin = lower, ymax = upper),
  #             alpha = 0.2, data = q2b_pred2) +
  geom_ribbon(data = q2b_ari_pred,
              aes(x = 1000 * orthoq1000, ymin = lower, ymax = upper),
              alpha = 0.25, fill = 'black') +
    geom_line(aes(x = 1000 * orthoq1000, y = yhat),
              data = q2b_ari_pred, colour = 'black') +
  # geom_ribbon(aes(x = log2_orthoq, ymin = lower, ymax = upper),
  #             alpha = 0.25) +
  # geom_line(aes(x = log2_orthoq, y = yhat), colour = 'black') +
  scale_x_continuous(limits = c(50, 2000),
                     breaks = seq(0, 2000, by = 250),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\n5,6-orthoquinone (ng/ml)',
       y = 'Day-2 methaemoglobin (%)\n'))

q2b_ari_pred |>
  mutate(yhat = log2(yhat),
         lower = log2(lower),
         upper = log2(upper)) |> 
    ggplot() +
    geom_point(data = prima_base_d2,
               colour = 'black', alpha = 0.3,
               aes(orthoq, y = log2(methb_max_dayid))) +
    geom_ribbon(aes(x = 1000 * orthoq1000, ymin = lower, ymax = upper),
                alpha = 0.25, fill = 'black') +
    geom_line(aes(x = 1000 * orthoq1000, y = yhat),
              colour = 'black') +
    scale_x_continuous(limits = c(50, 2000),
                       breaks = seq(0, 2000, by = 250),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = log2(c(1, 16)),
                       breaks = log2(2^seq(1, 15, by = 1)),
                       labels = 2^seq(1, 15, by = 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = '\n5,6-orthoquinone (ng/ml)',
         y = 'Day-2 methaemoglobin (%)\n')

# Non-methaemoglobin ------------------------------------------------------
prima_base_d2_nonmet <- prima_base |>
  filter(day == 'Day 2') |>
  drop_na(orthoq, cyp2d6) |> 
  mutate(log2_orthoq = log2(orthoq),
         log2_methb = log2(methb_max_dayid),
         log2_methb_auc = log2(methb_auc_dayid),
         methb = methb_max_dayid,
         methb_auc = methb_auc_dayid,
         mr_recip = 1/metab_ratio,
         log2_mr_recip = log2(mr_recip),
         orthoq1000 = orthoq/1000,
         cyp_cat1_bin = if_else(cyp_score1 >= 1.25, 'Normal or ultrarapid', 'Intermediate'))
dd_nonmet <- datadist(prima_base_d2_nonmet)
options(datadist = 'dd_nonmet')

# T-test
cypcat_oqt_lm <- lm(log2_orthoq ~ cyp_cat1_bin, data = prima_base_d2_nonmet); anova(cypcat_oqt_lm)
cypcat_oqt <- ols(log2_orthoq ~ cyp_cat1_bin, data = prima_base_d2_nonmet); anova(cypcat_oqt)
performance::check_model(cypcat_oqt_lm)
performance::check_heteroscedasticity(cypcat_oqt_lm)
performance::check_normality(cypcat_oqt_lm)
residualscypcat_oqt_lm <- residuals(cypcat_oqt_lm)
qqnorm(residualscypcat_oqt_lm); qqline(residualscypcat_oqt_lm, col = "red"); shapiro.test(residualscypcat_oqt_lm)

ttest_cypcat_oqt <- t.test(log2_orthoq ~ cyp_cat1_bin, data = prima_base_d2_nonmet); 
performance::check_heteroscedasticity(cypcat_oqt)
residualst <- residuals(cypcat_oqt)
qqnorm(residualst); qqline(residualst, col = "red"); shapiro.test(residualst)

ttest_mean_oqt <- ttest_cypcat_oqt$estimate
(ttest_df_oq <- tibble(
  cyp_cat1_bin = c('Intermediate', 'Normal or ultrarapid'),
  oq = ttest_mean_oq,
  oq2 = 2^oq
))

prima_base_d2_nonmet |>
  ggplot() +
  geom_jitter(aes(cyp_cat1_bin, log2_orthoq), width = 0.05,
                colour = 'black', alpha = 0.3, size = 2) +
  geom_point(aes(x = cyp_cat1_bin, y = oq),
             data = ttest_df_oq, shape = 5, size = 2) +
  geom_boxplot(aes(cyp_cat1_bin, log2_orthoq),
               varwidth = TRUE, fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.7, linewidth = 0.3, alpha = 0.3) +
  scale_y_continuous(limits = log2(c(64, 8192)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity',
       y = '5,6-orthoquinone (ng/ml)\n')

# Heteroscedastic 1
cypcat_oq_cont <- ols(log2_orthoq ~ cyp_score1, data = prima_base_d2_nonmet); anova(cypcat_oq_cont)
cypcat_oq_cont_lm <- lm(log2_orthoq ~ cyp_score1, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_lm)
performance::check_model(cypcat_oq_cont_lm)
performance::check_heteroscedasticity(cypcat_oq_cont_lm)
performance::check_normality(cypcat_oq_cont_lm)
residualscypcat_oq_cont_lm <- residuals(cypcat_oq_cont_lm)
qqnorm(residualscypcat_oq_cont_lm); qqline(residualscypcat_oq_cont_lm, col = "red"); shapiro.test(residualscypcat_oq_cont_lm)

cypcat_oq_cont_pred <- Predict(cypcat_oq_cont,
                     cyp_score1 = seq(0, 3.5, by = 0.01),
                     ref.zero = F) |>
  as_tibble(); cypcat_oq_cont_pred

prima_base_d2_nonmet |> 
  ggplot() +
  geom_point(aes(cyp_score1, log2(orthoq)),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score1, ymin = lower, ymax = upper),
              data = cypcat_oq_cont_pred, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score1, y = yhat),
            data = cypcat_oq_cont_pred, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 3, by = 0.25),
                     limits = c(0.15, 3.15),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(64, 8192)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = '5,6-orthoquinone (ng/ml)\n')

# Heteroscedastic 2
cypcat_oq_cont2 <- ols(log2_orthoq ~ cyp_score2, data = prima_base_d2_nonmet); anova(cypcat_oq_cont2)
cypcat_oq_cont_lm2 <- lm(log2_orthoq ~ cyp_score2, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_lm2)
performance::check_model(cypcat_oq_cont_lm2)
performance::check_heteroscedasticity(cypcat_oq_cont_lm2)
performance::check_normality(cypcat_oq_cont_lm2)
residualscypcat_oq_cont_lm2 <- residuals(cypcat_oq_cont_lm2)
qqnorm(residualscypcat_oq_cont_lm2); qqline(residualscypcat_oq_cont_lm2, col = "red"); shapiro.test(residualscypcat_oq_cont_lm2)

cypcat_oq_cont_pred2 <- Predict(cypcat_oq_cont2,
                               cyp_score2 = seq(0, 3.5, by = 0.01),
                               ref.zero = F) |>
  as_tibble(); cypcat_oq_cont_pred2

prima_base_d2_nonmet |> 
  ggplot() +
  geom_point(aes(cyp_score2, log2(orthoq)),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score2, ymin = lower, ymax = upper),
              data = cypcat_oq_cont_pred2, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score2, y = yhat),
            data = cypcat_oq_cont_pred2, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 3, by = 0.25),
                     limits = c(0.15, 3.15),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(64, 8192)),
                     breaks = log2(2^seq(1, 15, by = 1)),
                     labels = 2^seq(1, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = '5,6-orthoquinone (ng/ml)\n')

# CYP2D6 non-normality of residuals
cypcat_oq_cont_rev <- ols(cyp_score1 ~ orthoq, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_rev)
cypcat_oq_cont_rev_lm <- lm(cyp_score1 ~ orthoq, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_rev_lm)
performance::check_model(cypcat_oq_cont_rev_lm)
performance::check_heteroscedasticity(cypcat_oq_cont_rev_lm)
performance::check_normality(cypcat_oq_cont_rev_lm)
residualscypcat_oq_cont_rev_lm <- residuals(cypcat_oq_cont_rev_lm)
qqnorm(residualscypcat_oq_cont_rev_lm); qqline(residualscypcat_oq_cont_rev_lm, col = "red"); shapiro.test(residualscypcat_oq_cont_rev_lm)

cypcat_oq_cont_pred_rev <- Predict(cypcat_oq_cont_rev,
                                   orthoq = seq(20, 8000, by = 1),
                                   ref.zero = F) |>
  as_tibble(); cypcat_oq_cont_pred_rev

prima_base_d2_nonmet |> 
  ggplot() +
  geom_point(aes(y = cyp_score1, x = orthoq),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(x = orthoq, ymin = lower, ymax = upper),
              data = cypcat_oq_cont_pred_rev, alpha = 0.25, fill = 'black') +
  geom_line(aes(x = orthoq, y = yhat),
            data = cypcat_oq_cont_pred_rev) +
  # scale_x_continuous(breaks = seq(0.25, 3, by = 0.25),
  #                    limits = c(0.15, 3.15),
  #                    expand = expansion(mult = c(0, 0))) +
  # scale_y_continuous(limits = log2(c(64, 8192)),
  #                    breaks = log2(2^seq(1, 15, by = 1)),
  #                    labels = 2^seq(1, 15, by = 1),
  #                    expand = expansion(mult = c(0, 0))) +
  labs(y = 'CYP2D6 activity score\n',
       x = '\n5,6-orthoquinone (ng/ml)')

cypcat_oq_cont_rev2 <- ols(cyp_score2 ~ orthoq, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_rev2)
cypcat_oq_cont_rev_lm2 <- lm(cyp_score2 ~ orthoq, data = prima_base_d2_nonmet); anova(cypcat_oq_cont_rev_lm2)
performance::check_model(cypcat_oq_cont_rev_lm2)
performance::check_heteroscedasticity(cypcat_oq_cont_rev_lm2)
performance::check_normality(cypcat_oq_cont_rev_lm2)
residualscypcat_oq_cont_rev_lm2 <- residuals(cypcat_oq_cont_rev_lm2)
qqnorm(residualscypcat_oq_cont_rev_lm); qqline(residualscypcat_oq_cont_rev_lm, col = "red"); shapiro.test(residualscypcat_oq_cont_rev_lm)

# Reciprocal of metabolic ratio i.e., 1/(urine pq/urin oq)
# T-test
cypcat_mr_lm <- lm(log2_mr_recip ~ cyp_cat1_bin, data = prima_base_d2_nonmet); anova(cypcat_mr_lm)
cypcat_mr <- ols(log2_mr_recip ~ cyp_cat1_bin, data = prima_base_d2_nonmet); anova(cypcat_mr)
performance::check_model(cypcat_mr_lm)
performance::check_heteroscedasticity(cypcat_mr_lm)
performance::check_normality(cypcat_mr_lm)
residualscypcat_mr_lm <- residuals(cypcat_mr_lm)
qqnorm(residualscypcat_mr_lm); qqline(residualscypcat_mr_lm, col = "red"); shapiro.test(residualscypcat_mr_lm)

ttest_cypcat_mr <- t.test(log2_mr_recip ~ cyp_cat1_bin, data = prima_base_d2_nonmet)

ttest_mean_mr <- ttest_cypcat_mr$estimate
(ttest_df_mr <- tibble(
  cyp_cat1_bin = c('Intermediate', 'Normal or ultrarapid'),
  mr = ttest_mean_mr,
  mr2 = 2^mr
))

prima_base_d2_nonmet |>
  ggplot() +
  geom_jitter(aes(cyp_cat1_bin, log2_mr_recip), width = 0.05,
              colour = 'black', alpha = 0.3, size = 2) +
  geom_point(aes(x = cyp_cat1_bin, y = mr),
             data = ttest_df_mr, shape = 5, size = 2) +
  geom_boxplot(aes(cyp_cat1_bin, log2_mr_recip),
               varwidth = TRUE, fill = 'transparent', outlier.shape = NA,
               outlier.alpha = 0.5, outlier.size = 0.7,
               width = 0.7, linewidth = 0.3, alpha = 0.3) +
  scale_y_continuous(limits = log2(c(0.07, 10)),
                     breaks = log2(2^seq(-5, 15, by = 1)),
                     labels = 2^seq(-5, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity',
       y = 'Reciprocal of metabolic ratio\n')

# Homoscedastic + normal errors!
cypcat_mr_cont <- ols(log2_mr_recip ~ cyp_score1, data = prima_base_d2_nonmet); anova(cypcat_mr_cont)
cypcat_mr_cont_lm <- lm(log2_mr_recip ~ cyp_score1, data = prima_base_d2_nonmet); anova(cypcat_mr_cont_lm)
performance::check_model(cypcat_mr_cont_lm)
performance::check_heteroscedasticity(cypcat_mr_cont_lm)
performance::check_normality(cypcat_mr_cont_lm)
residualscypcat_mr_cont_lm <- residuals(cypcat_mr_cont_lm)
qqnorm(residualscypcat_mr_cont_lm); qqline(residualscypcat_mr_cont_lm, col = "red"); shapiro.test(residualscypcat_mr_cont_lm)

cypcat_mr_cont_pred <- Predict(cypcat_mr_cont,
                               cyp_score1 = seq(0, 3.5, by = 0.01),
                               ref.zero = F) |>
  as_tibble(); cypcat_mr_cont_pred

prima_base_d2_nonmet |> 
  ggplot() +
  geom_point(aes(cyp_score1, log2_mr_recip),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score1, ymin = lower, ymax = upper),
              data = cypcat_mr_cont_pred, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score1, y = yhat),
            data = cypcat_mr_cont_pred, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 3, by = 0.25),
                     limits = c(0.15, 3.15),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.07, 10)),
                     breaks = log2(2^seq(-5, 15, by = 1)),
                     labels = 2^seq(-5, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = 'Reciprocal of metabolic ratio\n')

# Homoscedastic + normal error! (2)
cypcat_mr_cont2 <- ols(log2_mr_recip ~ cyp_score2, data = prima_base_d2_nonmet); anova(cypcat_mr_cont2)
cypcat_mr_cont_lm2 <- lm(log2_mr_recip ~ cyp_score2, data = prima_base_d2_nonmet); anova(cypcat_mr_cont_lm2)
performance::check_model(cypcat_mr_cont_lm2)
performance::check_heteroscedasticity(cypcat_mr_cont_lm2)
performance::check_normality(cypcat_mr_cont_lm2)
residualscypcat_mr_cont_lm2 <- residuals(cypcat_mr_cont_lm2)
qqnorm(residualscypcat_mr_cont_lm2); qqline(residualscypcat_mr_cont_lm2, col = "red"); shapiro.test(residualscypcat_mr_cont_lm2)

cypcat_mr_cont_pred2 <- Predict(cypcat_mr_cont2,
                                cyp_score2 = seq(0, 3.5, by = 0.01),
                               ref.zero = F) |>
  as_tibble(); cypcat_mr_cont_pred2

prima_base_d2_nonmet |> 
  ggplot() +
  geom_point(aes(cyp_score2, log2_mr_recip),
             colour = 'black', alpha = 0.3) +
  geom_ribbon(aes(cyp_score2, ymin = lower, ymax = upper),
              data = cypcat_mr_cont_pred2, alpha = 0.25, fill = 'black') +
  geom_line(aes(cyp_score2, y = yhat),
            data = cypcat_mr_cont_pred2, fill = 'black') +
  scale_x_continuous(breaks = seq(0.25, 3, by = 0.25),
                     limits = c(0.15, 3.15),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.07, 10)),
                     breaks = log2(2^seq(-5, 15, by = 1)),
                     labels = 2^seq(-5, 15, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nCYP2D6 activity score',
       y = 'Reciprocal of metabolic ratio\n')











