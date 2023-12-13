library(brms)
library(cmdstanr)
options(mc.cores = 4,
  brms.backend = "cmdstanr")
library(mediator)
library(tidyverse)
library(tidybayes)

# mediator model
m_model <- bf(m_01 ~ 1 + x + c)

# outcome model
y_model <- bf(y ~ 1 + x + m_01 + c + x*m_01)

# brms model
bmodel <- brm(
  data = mediation_example, 
  family = bernoulli(),
  mvbf(m_model, y_model, rescor = FALSE),
  cores = 4)

# set up data frame for predictions
ndm <- tibble(x = 0:1, m_01 = 0.504, c = 1.015)

draws <- as_draws_df(bmodel)

# add predicted estimates
pd <- ndm %>% 
  add_epred_draws(bmodel)

# summarize
pd |> 
  filter(.category == "y") |> 
  pivot_wider(id_cols = .draw, 
    names_from = x, values_from = .epred) |> 
  mutate(diff = `1` - `0`) |> 
  median_qi(diff)

# Could also try by hand
effects <- draws %>% 
  mutate(
    cde = `b_y_x` + (0.504 * `b_y_x:m_01`),
    nie = (`b_y_m_01` * `b_m01_x`) + 
      (`b_y_x:m_01` * `b_m01_x`) ) %>%
  pivot_longer(cde:nie) %>% 
  group_by(name) %>% 
  median_qi(value)
effects

%>%
  pivot_longer(cde:nde)