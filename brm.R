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
# mean value of mediator and covariate
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

# Could also try by hand on the log odds scale
effects <- draws %>% 
  mutate(
    cde = `b_y_x` + (0.504 * `b_y_x:m_01`),
    nde = `b_y_x` + 
      log(1 + exp(`b_m01_Intercept` + `b_m01_c` + 
                     `b_y_m_01` + `b_y_x:m_01`)) -
      log(1 + exp(`b_m01_Intercept` + `b_m01_c` +
                    `b_y_m_01`)),
    nie = (`b_y_m_01` * `b_m01_x`) + 
      (`b_y_x:m_01` * `b_m01_x`) ) %>%
  pivot_longer(cde:nie) %>% 
  group_by(name) %>% 
  median_qi(value)
effects

+
  log(1 + exp(beta0 + beta1*0 + beta2*(1) + beta3*1*0 + theta2 + theta3*1 + theta6)) -
  log(1 + exp(beta0 + beta1*0 + beta2*(1) + beta3*1*0 + theta2 + theta3*0 + theta6))

%>%
  pivot_longer(cde:nde)


mreg <- glm(m_01 ~ x + c, family = "binomial",
  data = mediation_example)

yreg <- glm(y ~ x + m_01 + c + x*m_01, 
  family = "binomial", data = mediation_example)

regmedint_obj <- regmedint(data = mediation_example,
                           yvar = "y", avar = "x", mvar = "m_01", cvar = c("c"),
                           ## Values at which effects are evaluated
                           a0 = 0, a1 = 1, m_cde = 0.504, c_cond = 1,
                           ## Model types
                           mreg = "logistic", yreg = "logistic",
                           ## Additional specification
                           interaction = TRUE, casecontrol = FALSE)

beta0 = mreg$coefficients[['(Intercept)']]
beta1 = mreg$coefficients[['x']]
beta2 = mreg$coefficients[['c']]
beta3 = 0
theta0 = yreg$coefficients[['(Intercept)']]
theta1 = yreg$coefficients[['x']]
theta2 = yreg$coefficients[['m_01']]
theta3 = yreg$coefficients[['x:m_01']]
theta4 = yreg$coefficients[['c']]
theta5 = 0
theta6 = 0

pnde <- ((theta1) * (1)) +
  log(1 + exp(beta0 + beta1*0 + beta2*(1) + beta3*1*0 + theta2 + theta3*1 + theta6)) -
  log(1 + exp(beta0 + beta1*0 + beta2*(1) + beta3*1*0 + theta2 + theta3*0 + theta6))
