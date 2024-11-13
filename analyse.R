
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
prima <- read_rds(file = here('data', 'prima.rds'))

# Graphical, non-summary --------------------------------------------------
# Methaemoglobin over timepoints, by day
prima %>%
  drop_na(methb) %>% # Do better here!
  ggplot() +
  geom_line(aes(x = timepoint,
                y = methb,
                group = patid),
            alpha = 0.3, linewidth = 0.8) +
  facet_grid(~day) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 1),
                     expand = expansion(mult = c(0, 0)),
                     # sec.axis = sec_axis(~., name = "Days since enrolment\n",
                     #                     labels = NULL)
                     ) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoint 0')

# Methaemoglobin over timepoints, connected
prima %>%
  drop_na(methb) %>% # Do better here!
  ggplot() +
  annotate("rect", xmin = c(8, 17), xmax = c(9, 18), ymin = -Inf, ymax = Inf,
           fill = '#918580', alpha = 0.3) +
  geom_vline(xintercept = c(0, 9, 18),
             linetype = 'dotted') +
  geom_line(aes(x = timepoint_cont,
                y = methb,
                group = patid),
            alpha = 0.3, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 26),
                     breaks = seq(0, 26, by = 3),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 20, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = '\nTimepoint',
       y = 'Methaemoglobin (%)\n',
       caption = 'Daily pre-dose measurements at timepoints 0, 9, and 18')

# Last timepoint but base timepoint?
# Important for overall time-series

# Graphical, summary ------------------------------------------------------


# Models ------------------------------------------------------------------


