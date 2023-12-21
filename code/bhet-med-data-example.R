#  program:  bhet-med-analysis-ex.R
#  task:     mediation models for respiratory outcomes
#  input:    bhet-master; indoor temp files
#  output:   
#  project:  BHET
#  author:   sam harper \ 2023-12-21


##  0 Load needed packages ----
library(here)
library(tidyverse)
library(tidybayes)
library(haven)
library(readxl)
library(kableExtra)
library(modeldb)
library(brms)
library(cmdstanr)
library(modelr)
library(osfr)
library(modelsummary)
library(bayesplot)
library(patchwork)
library(marginaleffects)
library(labelled)

# Use the cmdstanr backend for Stan
# You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
options(mc.cores = 4,
        brms.backend = "cmdstanr")

## download data files from OSF (data-clean component)
# dir.create("data-clean")
# bhet_project <- osf_retrieve_node("b4wze")
# bhet_project %>%
#   osf_ls_files("Master Dataset (Season 4)",
#                pattern = "dta") %>%
#   osf_download(path = here("data-clean"),
#     conflicts = "overwrite")


## 1 Read in dataset, limit to resp vars ----
d <- read_dta(here("data-clean", 
                   "BHET_master_data_15Dec2023.dta"), 
  col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, 
                ban_status_2019, ban_status_2020, 
                ban_status_2021, ban_status_no, 
                ban_status_composite, smoking,
                freq_cough, freq_phlegm,
                freq_wheezing, freq_breath,
                freq_no_chest, health_selfreport,
                freq_exercising, hours_sleep)) 

# define main outcome (as binary)
d1 <- d %>% drop_na() %>%
  mutate(resp = if_else(
    freq_cough < 3 |
    freq_phlegm < 3 |
    freq_wheezing < 3 |
    freq_breath < 3 |
    freq_no_chest < 3, 1, 0),
    fphealth = if_else(health_selfreport > 2, 1, 0, 
                       missing = NA),
    srh = health_selfreport,
    fex = freq_exercising,
    hsleep = hours_sleep,
    year = if_else(wave==1, 2018, 
      if_else(wave==2, 2019,
        if_else(wave==4, 2021, 0))),
    cohort_year = if_else(
      ban_status_composite==1, 2019, 
      if_else(ban_status_composite==2, 2020, 
              if_else(ban_status_composite==3, 2021, 2022))),
    treat = ifelse(year >= cohort_year, 1, 0),
    cohort_year = ifelse(cohort_year == 2022,-Inf, 
                         cohort_year)) %>%
  # relabel last cohort year 
  # treatment cohort dummies
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2019,2020,2021), 
    remove_original = F) %>%
  # wave dummies
  add_dummy_variables(year, 
    values=c(2018,2019,2021), remove_original = F)

dt <- readxl::read_excel(here("data-clean", 
  "HH_indoor_temperature_HEATING_SEASON.xlsx")) %>%
  select(ptc_id, wave, min_h, med_h, max_h,
         sd_h, iqr_h)

# merge indoor temp data with resp data
d2 <- d1 %>%
  left_join(dt, by = join_by(ptc_id,wave))

# limit sample

d3 <- d2 %>% 
  select(starts_with(c("year","cohort")), 
  "ID_VILLAGE","resp","treat", "med_h",
  "min_h", "max_h", "sd_h", "iqr_h", 
  "smoking", "srh", "fphealth", "fex",
  "hsleep") %>%
  
  # limit to complete cases
  drop_na() %>% 
  
  # create unique continuous village Id
  group_by(ID_VILLAGE) %>%
  mutate(v_id = cur_group_id()) %>%
  ungroup()

# remove labels
val_labels(d3) <- NULL

# write to data-clean
saveRDS(d3, file = here("data-clean", "d3"))

