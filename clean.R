
# Code author : Ihsan Fadilah
# Project     : Levels of 5,6-orthoquinone derivative of primaquine in urine
#               correlate with CYP2D6 genotype-predicted metabolic activity
#               score
# About       : Tidy data for analysis of methaemoglobin data in the first 3
#               days of enrolment

# Setup -------------------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(here)
library(janitor)
library(lubridate)

raw_dat <- read_excel(path = here('data', 'PRIMA_MetHb Mario.xlsx'),
                      sheet = 'MetHb')

# Demographics ------------------------------------------------------------
raw_demog <- raw_dat[, 1:9]
raw_demog <- raw_demog[3:nrow(raw_demog), ]
raw_demog <- raw_demog |> row_to_names(1)

clean_demog <- raw_demog |> 
  clean_names() |> 
  select(subject_number, group_is_participant_randomised, visit_date,
         day, month, year, what_is_the_sex_of_the_subject) |> 
  rename(patid = subject_number,
         group = group_is_participant_randomised,
         sex = what_is_the_sex_of_the_subject,
         birth_day = day,
         birth_month = month,
         birth_year = year) |> 
  mutate(clinid = NA,
         age = NA,
         weight = NA) |> 
  select(patid, clinid, age, sex, weight, group, birth_day, birth_month,
         birth_year, visit_date) |> 
  mutate(birth_date = mdy(paste(birth_month, birth_day, birth_year)),
         visit_date = ymd(visit_date),
         age = as.integer(interval(birth_date, visit_date) / years(1))) |> 
  select(-c(birth_day, birth_month, birth_year,
            visit_date, birth_date)) |> 
  mutate(sex = as.factor(sex),
         group = as.factor(group)) |> 
  arrange(patid)

# Methaemoglobin ----------------------------------------------------------
# Visit 1 (day 0)
raw_methb_day0 <- raw_dat[, c(1, 10:36)]
raw_methb_day0 <- raw_methb_day0[3:nrow(raw_methb_day0), ]
raw_methb_day0_wide <- raw_methb_day0 |> row_to_names(1) |> clean_names()

raw_methb_day0_long_timepoint <- raw_methb_day0_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('test_date_test_time'),
    names_to = 'timepoint',
    values_to = 'datetime'
  ) |>
  select(-c(starts_with('met_hb_in_percent'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 0) |> 
  select(subject_number, day, datetime, timepoint)

raw_methb_day0_long_methb <- raw_methb_day0_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('met_hb_in_percent'),
    names_to = 'timepoint',
    values_to = 'methb'
  ) |>
  select(-c(starts_with('test_date_test_time'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 0) |> 
  select(subject_number, day, timepoint, methb)

raw_methb_day0_long <- full_join(raw_methb_day0_long_timepoint,
                                 raw_methb_day0_long_methb,
                                 by = c('subject_number',
                                        'day',
                                        'timepoint'))

# Visit 2 (day 1)
raw_methb_day1 <- raw_dat[, c(1, 37:63)]
raw_methb_day1 <- raw_methb_day1[3:nrow(raw_methb_day1), ]
raw_methb_day1_wide <- raw_methb_day1 |> row_to_names(1) |> clean_names()

raw_methb_day1_long_timepoint <- raw_methb_day1_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('test_date_test_time'),
    names_to = 'timepoint',
    values_to = 'datetime'
  ) |>
  select(-c(starts_with('met_hb_in_percent'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 1) |> 
  select(subject_number, day, datetime, timepoint)

raw_methb_day1_long_methb <- raw_methb_day1_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('met_hb_in_percent'),
    names_to = 'timepoint',
    values_to = 'methb'
  ) |>
  select(-c(starts_with('test_date_test_time'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 1) |> 
  select(subject_number, day, timepoint, methb)

raw_methb_day1_long <- full_join(raw_methb_day1_long_timepoint,
                                 raw_methb_day1_long_methb,
                                 by = c('subject_number',
                                        'day',
                                        'timepoint'))

# Visit 3 (day 2)
raw_methb_day2 <- raw_dat[, c(1, 64:90)]
raw_methb_day2 <- raw_methb_day2[3:nrow(raw_methb_day2), ]
raw_methb_day2_wide <- raw_methb_day2 |> row_to_names(1) |> clean_names()

raw_methb_day2_long_timepoint <- raw_methb_day2_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('test_date_test_time'),
    names_to = 'timepoint',
    values_to = 'datetime'
  ) |>
  select(-c(starts_with('met_hb_in_percent'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 2) |> 
  select(subject_number, day, datetime, timepoint)

raw_methb_day2_long_methb <- raw_methb_day2_wide |> 
  rename(time_point_on_day_1 = time_point_on_day,
         test_date_test_time_1 = test_date_test_time,
         met_hb_in_percent_1 = met_hb_in_percent) |> 
  pivot_longer(
    cols = starts_with('met_hb_in_percent'),
    names_to = 'timepoint',
    values_to = 'methb'
  ) |>
  select(-c(starts_with('test_date_test_time'),
            starts_with('time_point_on_day'))) |> 
  mutate(time_extract = substr(timepoint, nchar(timepoint), nchar(timepoint))) |> 
  select(-timepoint) |> 
  rename(timepoint = time_extract) |> 
  mutate(day = 2) |> 
  select(subject_number, day, timepoint, methb)

raw_methb_day2_long <- full_join(raw_methb_day2_long_timepoint,
                                 raw_methb_day2_long_methb,
                                 by = c('subject_number',
                                        'day',
                                        'timepoint'))

# Merged: days 0-2
clean_methb <- bind_rows(raw_methb_day0_long,
                         raw_methb_day1_long,
                         raw_methb_day2_long) |> 
  rename(patid = subject_number) |> 
  mutate(timepoint = as.integer(timepoint)) |> 
  arrange(patid, day, timepoint) |> 
  select(patid, datetime, day, timepoint, methb)


# Primaquine treatment ----------------------------------------------------
# Visit 1 (day 0)
raw_pq_day0 <- raw_dat[, c(1, 91:106)]
raw_pq_day0 <- raw_pq_day0[3:nrow(raw_pq_day0), ]
raw_pq_day0_wide <- raw_pq_day0 |> row_to_names(1) |> clean_names()

# Visit 2 (day 1)
raw_pq_day1 <- raw_dat[, c(1, 107:122)]
raw_pq_day1 <- raw_pq_day1[3:nrow(raw_pq_day1), ]
raw_pq_day1_wide <- raw_pq_day1 |> row_to_names(1) |> clean_names()

# Visit 3 (day 2)
raw_pq_day2 <- raw_dat[, c(1, 123:160)]
raw_pq_day2 <- raw_pq_day2[3:nrow(raw_pq_day2), ]
raw_pq_day2_wide <- raw_pq_day2 |> row_to_names(1) |> clean_names()

# ...

# Other lab. measures (CYP2D6, 5,6-orthoquinone, G6PD, ...) ---------------
  

# Overall dataset ---------------------------------------------------------
prima <- full_join(clean_demog,
                   clean_methb,
                   by = 'patid') |> 
  mutate(datetime = ymd_hm(datetime),
         methb = as.numeric(methb),
         day = as.factor(day),
         timepoint = timepoint - 1,
         day = case_when(day == '0' ~ 'Day 0',
                         day == '1' ~ 'Day 1',
                         day == '2' ~ 'Day 2'),
         timepoint_cont = case_when(
           day == 'Day 0' ~ timepoint,
           day == 'Day 1' ~ timepoint + 9,
           day == 'Day 2' ~ timepoint + 18
         ))
  
write_rds(prima,
          file = here('data', 'prima.rds'))
  
  
  
  
  
  
  
  
  
  
  
  