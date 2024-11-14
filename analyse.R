
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
  mutate(cyp_cat1 = factor(cyp_cat1, levels = c('Intermediate', 'Normal')))

prima_base <- prima0to2 |>
  filter(timepoint == 0) |> 
  select(-c(clinid, weight, pqmg_per_tab, pqmgday, pqmgkgday))

# Graphical, non-summary --------------------------------------------------
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

# Last timepoint but base timepoint?
# Important for overall time-series

# Graphical, summary ------------------------------------------------------
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

# Boxplot, by day
prima0to2 |> 
  filter(timepoint == 0) |> 
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

# Lines, by day and CYP2D6 group
prima0to2 |>
  drop_na(cyp_cat1) |> 
  filter(timepoint == 0) |> 
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

# Maximum methaemoglobin, by day and some factor
## CYP2D6
prima0to2 |> 
  filter(timepoint == 0) |> 
  drop_na(cyp_cat1) |> 
  ggplot() +
  geom_beeswarm(aes(y = log2(methb_max_dayid), x = day, colour = cyp_cat1),
                dodge.width = 0.65, groupOnX = F, alpha = 0.3, size = 1.25) +
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

## 5,6-orthoquinone
prima0to2 |> 
  filter(timepoint == 0) |> 
  mutate(log_orthoq = log(orthoq)) |> 
  ggplot() +
  geom_point(aes(x = log_orthoq, y = log2(methb_max_dayid), colour = day),
             size = 1.5, alpha = 0.35) +
  geom_smooth(aes(x = log_orthoq, y = log2(methb_max_dayid), colour = day),
              method = "lm", se = F, linewidth = 0.8) +
  scale_x_continuous(limits = c(4.5, 8),
                     breaks = seq(0, 20, by = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = log2(c(0.6, 14)),
                     breaks = log2(c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)), 
                     labels = c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nNatural logarithm of 5,6-orthoquinone',
       y = 'Maximum methaemoglobin (%)\n',
       colour = '') +
  scale_colour_manual(
    values = c('Day 0' = '#536E85',
               'Day 1' = '#E3932B',
               'Day 2' = '#C82F46')
  ) +
  theme(legend.key.width = unit(0.5, 'cm'))

## G6PD
prima0to2 |> 
  filter(timepoint == 0) |> 
  ggplot() +
  geom_point(aes(x = g6pd1, y = log2(methb_max_dayid), colour = day),
             size = 1.5, alpha = 0.35) +
  geom_smooth(aes(x = g6pd1, y = log2(methb_max_dayid), colour = day),
              method = "lm", se = F, linewidth = 0.8) +
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
  scale_colour_manual(
    values = c('Day 0' = '#536E85',
               'Day 1' = '#E3932B',
               'Day 2' = '#C82F46')
  ) +
  theme(legend.key.width = unit(0.5, 'cm'))

# Maximum methaemoglobin on day 2, by some factor
## Age


## Sex


## CYP2D6


## G6PD


## Primaquine daily dose


# Models ------------------------------------------------------------------
























