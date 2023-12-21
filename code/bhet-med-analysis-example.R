#  program:  bhet-med-analysis-ex.R
#  task:     mediation models for respiratory outcomes
#  input:    bhet-master; indoor temp files
#  output:   
#  project:  BHET
#  author:   sam harper \ 2023-12-16


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

# read in data
d3 <- readRDS(here("data-clean", "d3"))

te_med <-
  brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + 
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
       prior(normal(0, 1), class = b),
       prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 3857,
      file = "code/fits/bhet-med-te_med")

# load if model already exists
te_med <- readRDS("code/fits/bhet-med-te_med.rds")

## check the chains
mcmc_trace(te_med, pars=c("b_Intercept", "b_year_2021")) +
  theme_classic() + ylab("")

## extract prior draws
pd <- prior_draws(te_med) %>%
  # rescale to absolute probabilities
  mutate(pd0 = inv_logit_scaled(Intercept),
         pd1 = inv_logit_scaled(Intercept + b)) %>%
  # treatment effect
  mutate(diff = pd1 - pd0)

## plot for intercept prior
plot_pd1 <- pd %>%
  ggplot(aes(x = pd1)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Probability of poor lung health", 
    y = NULL, subtitle = "Prior for Intercept") +
  theme_classic()

## plot for treatment effect prior
plot_pdd <- pd %>%
  ggplot(aes(x = diff)) + 
  stat_halfeye(color= "black", fill =  '#d95f02') +
  labs(x = "Difference in treated vs. control", 
    y = NULL, subtitle = "Prior for treatment effect") +
  theme_classic()

# Combined plot for priors
priors_lung <- (plot_pd1 | plot_pdd) +
  plot_annotation(
    title = "Priors for poor lung symptoms",
                  theme = theme_classic())

# Marginal predictions for Bayesian ETWFE (simple)

## load brms model
# te_med <- readRDS(here("code/fits", 
 # "bhet-med-te_med.rds"))

bme_pred <- predictions(
  te_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# plot of predicted probabilities by treatment
bme_pred_p <- bme_pred |>
  posterior_draws() |>
  ggplot(aes(x = draw, fill=factor(treat))) +
    stat_halfeye(slab_alpha = .5) + 
    annotate("text", x = 0.61, y = 0.7, 
           label="Control", color='#1b9e77') +
    annotate("text", x = 0.49, y = 0.95, 
           label="Treated", color='#d95f02') +
  scale_x_continuous(
    "Probability of poor respiratory symptoms", 
    limits=c(0.4,0.8)) +
  scale_y_continuous("Posterior Density") +
  scale_fill_manual(values = c('#1b9e77','#d95f02')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=12))

bme_avg <- slopes(
  te_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  ) 

# plot of treatment effect
bme_avg_p <- bme_avg |>
  posterior_draws() |>
  ggplot(aes(x = draw)) +
    stat_halfeye(slab_alpha = .5, fill = "#7570b3") +
    annotate("text", x = -0.106, y = 0.95, 
           label="Difference", color = '#7570b3') +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60") +
    scale_x_continuous("Marginal Effect", limits=c(-0.3,0.1)) +
    scale_y_continuous("") +
    theme_classic() + 
    theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=12))


f2 <- bme_pred_p + bme_avg_p + 
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Posterior distributions of marginal predictions: poor respiratory symptoms")
f2


# wrangle aggregate ATTs for model summary table
bmp <- data.frame(
  term = bme_pred$treat,
  estimate = bme_pred$estimate,
  conf.low = bme_pred$conf.low,
  conf.high = bme_pred$conf.high,
  std.error = abs(bme_pred$conf.high - 
                    bme_pred$conf.low) / (2 * 1.96)
)

bti <- data.frame(
  term = 2,
  estimate = bme_avg$estimate,
  conf.low = bme_avg$conf.low,
  conf.high = bme_avg$conf.high,
  std.error = abs(bme_avg$conf.high - 
                    bme_avg$conf.low) / (2 * 1.96)
)

bta <- bind_rows(bmp,bti) %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

gl <- data.frame()

betwfe_me_avg <- list(tidy = bta, glance = gl)
class(betwfe_me_avg) <- "modelsummary_list"

modelsummary(list("Total Effect" = betwfe_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


# by cohort
bme_c <- slopes(
  te_med, 
  newdata   = subset(d3, treat & cohort_year),
  variables = "treat", 
  by        = "cohort_year"
  )

bti_c <- data.frame(
  term = paste("ATT(g", bme_c$cohort_year, ")", sep=""), 
  estimate = bme_c$estimate,
  conf.low = bme_c$conf.low,
  conf.high = bme_c$conf.high,
  std.error = abs(bme_c$conf.high - 
                    bme_c$conf.low) / (2 * 1.96)
)

betwfe_me_c <- list(tidy = bti_c, glance = gl)
class(betwfe_me_c) <- "modelsummary_list"

modelsummary(list("Bayesian Cohort Average" = betwfe_me_c),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


# mediation model with interaction
cde_med <-
  brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + hsleep +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        treat:cohort_year_2019:year_2019:hsleep + 
        treat:cohort_year_2019:year_2021:hsleep +
        treat:cohort_year_2020:year_2021:hsleep +
        treat:cohort_year_2021:year_2021:hsleep +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 2326,
      file = "code/fits/bhet-med-cde_med")

# marginal effects
bcde_pred <- predictions(
  cde_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  )

bcde_avg <- slopes(
  cde_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  ) 

# wrangle aggregate ATTs for model summary table
bmpc <- data.frame(
  term = bcde_pred$treat,
  estimate = bcde_pred$estimate,
  conf.low = bcde_pred$conf.low,
  conf.high = bcde_pred$conf.high,
  std.error = abs(bcde_pred$conf.high - 
                    bcde_pred$conf.low) / (2 * 1.96)
)

btic <- data.frame(
  term = 2,
  estimate = bcde_avg$estimate,
  conf.low = bcde_avg$conf.low,
  conf.high = bcde_avg$conf.high,
  std.error = abs(bcde_avg$conf.high - 
                    bcde_avg$conf.low) / (2 * 1.96)
)

btca <- bind_rows(bmpc,btic) %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

gl <- data.frame()

bcde_me_avg <- list(tidy = btca, glance = gl)
class(bcde_me_avg) <- "modelsummary_list"

modelsummary(list("CDE" = bcde_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')

# both estimates
modelsummary(list("Total Effect" = betwfe_me_avg,
  "CDE" = bcde_me_avg),
  shape = term ~ model + statistic, 
  statistic = "conf.int",
  gof_omit ='._*')

# mediation model with interaction
cde_med_temp <-
  brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + min_h +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        treat:cohort_year_2019:year_2019:min_h + 
        treat:cohort_year_2019:year_2021:min_h +
        treat:cohort_year_2020:year_2021:min_h +
        treat:cohort_year_2021:year_2021:min_h +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 2326,
      file = "code/fits/bhet-med-cde_med_temp")



# does the policy affect self-reported health?
te_hsleep <-
  brm(data = d3, 
      family = gaussian(),
      hsleep ~ 1 + (1 | v_id) +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 265,
      file = "code/fits/te_hsleep")

# load if already run
te_temp_min <- readRDS("code/fits/te_temp_min.rds")

ndc <- subset(d3, treat==1) %>%
  mutate(mean_h=20)


# mediator model
m_model <- bf(hsleep ~ 1 + (1 | v_id) +
                treat:cohort_year_2019:year_2019 + 
                treat:cohort_year_2019:year_2021 +
                treat:cohort_year_2020:year_2021 +
                treat:cohort_year_2021:year_2021 +
                cohort_year_2019 + cohort_year_2020 +
                cohort_year_2021 + year_2019 + year_2021) +
                gaussian()

# outcome model
y_model <- bf(resp ~ 1 + (1 | v_id) + min_h + 
                treat:cohort_year_2019:year_2019 + 
                treat:cohort_year_2019:year_2021 +
                treat:cohort_year_2020:year_2021 +
                treat:cohort_year_2021:year_2021 +
                cohort_year_2019 + cohort_year_2020 +
                cohort_year_2021 + year_2019 + year_2021) +
                bernoulli()

# priors
priormed <- c(
  prior(normal(10, 5), class = Intercept, resp = minh),
  prior(normal(0, 5), class = b, resp = minh),
  prior(exponential(1), class = sd, resp = minh),
  prior(normal(0, 1), class = Intercept, resp = resp),
  prior(normal(0, 1), class = b, resp = resp),
  prior(exponential(1), class = sd, resp = resp))

medfit <- brm(
  m_model + y_model + set_rescor(FALSE),
  data = d3, iter = 2000, warmup = 1000,
  prior = priormed,
  chains = 4, cores = 4,
  sample_prior = "yes")


brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + min_h +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 265,
      file = "code/fits/bhet-resp-cde_med_min_ni")


ndm <- tibble(treat = 0:1, hsleep = 7.74)

draws <- as_draws_df(cde_med)

# add predicted estimates
pd <- ndm %>% 
  add_epred_draws(cde_med, allow_new_levels=TRUE)

# summarize
pd |> 
  filter(.category == "y") |> 
  pivot_wider(id_cols = .draw, 
    names_from = x, values_from = .epred) |> 
  mutate(diff = `1` - `0`) |> 
  median_qi(diff)