# Bayesian model for joystick reseponses
# Makes use of a Hamiltonian Monte Carlo sampler algorithm (MCMC) to approximate
# the posterior (distribution)
# https://www.rensvandeschoot.com/tutorials/brms-started/

library(here)
library(rstan)
library(brms)
library(ggdist)
library(modelr)
library(tidybayes)
library(ggdist)
library(tidyverse)

#for less aliased plots on windows
trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

#read data from online file
df <- read_csv('https://www.dropbox.com/s/klfyv2hk78ko8a4/joystick.csv?dl=1',
               col_types = cols("f", "f", "i", "f", "f", "f", "f", "f", "i", "i", "f", "i", "d", "d", "d"))
#df <- read_csv('D:\\Dropbox\\MyScience\\MyPostdoc\\exp\\global_local\\2-data\\joystick\\3-preprocessed\\joystick.csv',
#                col_types = cols("f", "f", "i", "f", "f", "f", "f", "f", "i", "i", "f", "i", "d", "d", "d"))
               

df_mdl <- df %>%
    filter(lateness > 0) #%>%
    #mutate(au = round(au, 1)) # round to nearest decimal point

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

zoib_model <- bf(
  au ~ lateness + (1+lateness|id),  # Beta dist mean
  phi ~ lateness + (1+lateness|id),  # Beta dist precision (1/variance : higher values = less variation)
  zoi ~ lateness + (1+lateness|id),  # Zero-one inflation (alpha): probability of binary rating
  coi ~ lateness + (1+lateness|id),  # Conditional one inflation (gamma): probability of 1 if binary rating
  family = zero_one_inflated_beta()
)

fit <- brm(zoib_model, 
           data   = df_mdl, 
           #warmup = 1000, 
           #iter   = 3000, 
           #chains = 2, 
           #inits  = "random",
           cores  = 4,  #the cores function tells STAN to make use more CPU cores simultaneously instead of just 1.
           file = here::here("zoib-brms")
           )

summary(fit)

#try monotonic
zoib_model_mono <- bf(
  au ~ mo(lateness) + (1+mo(lateness)|id),  # Beta dist mean
  phi ~ mo(lateness) + (1+mo(lateness)|id),  # Beta dist precision (1/variance : higher values = less variation)
  zoi ~ mo(lateness) + (1+mo(lateness)|id),  # Zero-one inflation (alpha): probability of binary rating
  coi ~ mo(lateness) + (1+mo(lateness)|id),  # Conditional one inflation (gamma): probability of 1 if binary rating
  family = zero_one_inflated_beta()
)

fit2 <- brm(zoib_model_mono, 
           data   = df_mdl, 
           #warmup = 1000, 
           #iter   = 3000, 
           #chains = 2, 
           #inits  = "random",
           cores  = 4,  #the cores function tells STAN to make use more CPU cores simultaneously instead of just 1.
           file = here::here("zoib-brms-mono")
)

get_variables(fit)

# plot lateness
fit %>% 
  gather_draws(b_lateness) %>% 
  ggplot(aes(y=.variable, x = .value)) + 
  stat_histinterval()

#plot more from Vuorre's blog
# Transform-then-summarise
posterior_samples(fit, pars = "b_")[,1:4] %>%
  mutate_at(c("b_phi_Intercept"), exp) %>% #precision was modeled on the log scale
  mutate_at(vars(-"b_phi_Intercept"), plogis) %>% #others modeled on logit scale, plogit is the inverse
  posterior_summary() %>%
  as.data.frame() %>%
  rownames_to_column("Parameter")

posterior_samples(fit, pars = "b_")[,5:8] %>%
  mutate_at(c("b_phi_lateness"), exp) %>% #precision was modeled on the log scale
  mutate_at(vars(-"b_phi_lateness"), plogis) %>% #others modeled on logit scale, plogit is the inverse
  posterior_summary() %>%
  as.data.frame() %>%
  rownames_to_column("Parameter")

#plot
plot(
  conditional_effects(fit, dpar = 'mu'),
  points = TRUE,
  points_args = list(width=.05, shape=5)
)

#b_Intercept is the global mean, and the r_id variables are offsets from the mean for each id
fit %>%
  spread_draws(r_id[id], term) %>%  # or gather_draws() for a long-list
  median_qi()


####
grid = df_mdl %>%
  modelr::data_grid(lateness, id)

fits = grid %>%
  add_fitted_draws(fit)

preds = grid %>%
  add_predicted_draws(fit)

df_mdl %>%
  ggplot(aes(y = au, x = lateness)) +
  stat_slab() + #real data or preds
  stat_interval(aes(y = .prediction), data = preds) +
  stat_pointinterval(aes(y = .value), data = fits, .width = c(.66, .95), position = position_nudge(x = 0.1)) +
  #geom_point() +
  facet_wrap(id~.) +
  scale_color_brewer()


###2

fits2 = grid %>%
  add_fitted_draws(fit2)

preds2 = grid %>%
  add_predicted_draws(fit2)

df_mdl %>%
  ggplot(aes(y = au, x = lateness)) +
  stat_interval(data = preds2, aes(y = .prediction), .width = c(.99, .95, .8, .5)) +
  #geom_jitter(size=.1, width=.2, height=.0, alpha=.25) +
  #stat_halfeye(data = fits2,aes(y = .value, fill=as.factor(lateness)), size = .1, scale = .8, position = position_nudge(x = 0.175)) +
  stat_pointinterval(data = fits2, aes(y = .value), .width = c(.66, .95), position = position_nudge(x = .375), size=.1) +
  facet_wrap(id~.) +
  scale_color_brewer() +
  theme_linedraw()+
  NULL

# check that line_ribbon is based on correct model fit
df_mdl %>%
  #mutate(id = reorder(id, au, mean)) %>%
  ggplot(aes(y = au, x = lateness)) +
    #stat_lineribbon(data = fits, aes(y = .value), .width = c(.66, .95), alpha=.5) +
    stat_lineribbon(data = preds, aes(y = .prediction), .width = c(.99, .95, .8, .5), size = .75, alpha=.375) +
    geom_jitter(size=.2, width=.2, height=.0, alpha=.5) +
    scale_fill_brewer() +
    facet_wrap(id~.) +
    theme_linedraw() +
    NULL
  

# compare fits for best model

# 1. estimate slope effect per participant
#4>3>2>1 or just 4>1

# 1.1 Use the transform-then-summarise method from Vuorre

# TODO: check mutate transformation for eaach type of parameter (maybe r_ do not
# need to be transformed?)
posterior_samples(fit2, pars = c("b_", "r_")) %>%
  mutate_at(c("b_phi_Intercept"), exp) %>% #precision was modeled on the log scale
  mutate_at(vars(-"b_phi_Intercept"), plogis) %>% #others modeled on logit scale, plogit is the inverse
  posterior_summary(probs = c(0.025, 0.975)) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter")

#todo: how to access monotonic effects in model summary?
#seems to be "simo_" variables, but not per id?

# 1.1b This gives more info for other parameters (phi, zoi, coi...)
#posterior_summary(fit2, pars = c("^b_", "^r_"), probs = c(0.025, 0.975) )

# 1.2 

# Get individual slope estimates
#https://vuorre.netlify.app/post/2017/03/21/bayes-factors-with-brms/

h <- hypothesis(fit, "lateness > 0", scope = "coef", group = "id")
table <- as.data.frame(h$hypothesis)

# Create a data frame with subject IDs and coefficients
coefs <- table %>%
  transmute(id = Group, 
            slope = Estimate,
            star = Star)
# Join to main data frame by Subject ID
df_mdl <- left_join(df_mdl, coefs, by="id")

# add star to preds
preds <- left_join(preds, coefs, by='id')

df_mdl %>%
  mutate(id = reorder(id, slope)) %>%
  ggplot(aes(y = au, x = lateness)) +
  #stat_lineribbon(data = fits, aes(y = .value), .width = c(.66, .95), alpha=.5) +
  stat_lineribbon(data = preds, aes(y = .prediction, linetype=star, color=star), .width = c(.99, .95, .8, .5), size = .75, alpha=.375) +
  geom_jitter(size=.2, width=.2, height=.0, alpha=.25) +
  scale_fill_brewer() +
  scale_linetype_manual(values=c("blank", "solid")) +
  scale_color_manual(values=c('black', 'black')) +
  facet_wrap(id~.) +
  theme_linedraw() +
  NULL


# Results of hypothesis() in a data.frame
#h$hypothesis %>% 
  # Obtain group indicators from original data
  #left_join(distinct(dat, Group = id, treatment)) %>% 
  ## Rename Group to id and reverse order for figure
  #mutate(id = factor(Group, levels = rev(1:50))) %>%
  ## Draw a forest plot with ggplot2
  #ggplot(aes(Estimate, id, col = treatment)) +
  #geom_errorbarh(aes(xmin = CI.Lower, xmax = CI.Upper)) +
  #geom_point()

#also lineribbon linetype by sig or not

summary(fit2)