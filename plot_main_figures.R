#==============================================================================#
#                                                                              #
#                   Plot main figures for the phenology paper                  #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(ggdist)
library(patchwork)
Sys.setenv(LANG = "en")


# Functions used later in the code
get_median_and_CI <- function(x, lower, upper) {
  output <- data.frame(y = median(x), 
                       ymin = quantile(x, probs = lower),
                       ymax = quantile(x, probs = upper))
  return(output)
}

get_probability_direction <- function(posterior) {
  pd <- ifelse(median(posterior) > 0, sum(posterior > 0) / length(posterior),
               sum(posterior < 0) / length(posterior))
  return(pd)
}

# Fig. 3: Sea-ice decline -------------------------------------------------

# ~ a. Run models ------------------------------------------------------------

sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F) %>%
  filter(year >= 1988)

day_retreat_s <- as.vector(scale(sea_ice_data$day_retreat))
ice_free_days_s <- as.vector(scale(sea_ice_data$ice_free_days))
time_s <- as.vector(scale(sea_ice_data$year))



# Fit the linear model of sea-ice availability against time (< 2min)
library(nimble)

model_sea_ice <- nimbleCode ({

  for (t in 1:T) {
    # sea-ice: day retreat
    mu_retreat[t] <- alpha[1] +
      alpha[2] * time_s[t]

    day_retreat_s[t] ~ dnorm(mu_retreat[t], sd = sigma_day_retreat)

    mu_ice_free_fays[t] <- beta[1] +
      beta[2] * time_s[t]

    ice_free_days_s[t] ~ dnorm(mu_ice_free_fays[t], sd = sigma_ice_free_days)
  }

  for (k in 1:2) {
    alpha[k] ~ dnorm(0, sd = 3)
    beta[k] ~ dnorm(0, sd = 3)
  }
  sigma_ice_free_days ~ dunif(0, 5)
  sigma_day_retreat ~ dunif(0, 5)

})

dat <- list(day_retreat_s = day_retreat_s,
            ice_free_days_s = ice_free_days_s)

my.constants <- list(T = length(time_s),
                     time_s = time_s)

inits <- function() list(alpha = rnorm(n = 2, mean = 0, sd = 0.5),
                         sigma_day_retreat = runif(1, min = 0, max = 1),
                         beta = rnorm(n = 2, mean = 0, sd = 0.5),
                         sigma_ice_free_days = runif(1, min = 0, max = 1))

inits_values <- list(inits(), inits())

# Parameters monitored
params <- c("alpha", "sigma_day_retreat", "beta", "sigma_ice_free_days")


start <- Sys.time() ; start
fit_model_sea_ice <- nimbleMCMC(code = model_sea_ice,
                                constants = my.constants,
                                data = dat,
                                inits = inits,
                                monitors = params,
                                niter = 20000,
                                nburnin = 5000,
                                thin = 5,
                                nchains = 2)
save(fit_model_sea_ice,
     file = "02_outpus/fit_sea_ice.RData")
end <- Sys.time() ; end - start


# ~ b. Plot ------------------------------------------------------------------

# load("02_outpus/fit_sea_ice.RData")
# res <- rbind(fit_model_sea_ice$chain1,
#              fit_model_sea_ice$chain2)

# Day retreat
day_retreat <- matrix(data = NA, nrow = dim(res)[1], ncol = length(time_s),
                 dimnames = list(1:dim(res)[1], sea_ice_df$year))
years <- sea_ice_df$year
for (t in 1:length(years)) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    day_retreat[k, t] <- res[k, "alpha[1]"] +
      res[k, "alpha[2]"] * time_s[t]
  }                               
}

df_plot_day_retreat_1 <- as.data.frame.table(day_retreat) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  mutate(value = value * sd(sea_ice_df$day_retreat) + mean(sea_ice_df$day_retreat)) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))

n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

day_retreat <- matrix(data = NA, nrow = n_iterations, ncol = length(time_s),
                      dimnames = list(1:n_iterations, sea_ice_df$year))
years <- sea_ice_df$year
for (t in 1:length(years)) {                 # year
  for (k in 1:n_iterations) {         # MCMC iteration
    day_retreat[k, t] <- res[iterations[k], "alpha[1]"] +
      res[iterations[k], "alpha[2]"] * time_s[t]
  }                               
}

df_plot_day_retreat_2 <- as.data.frame.table(day_retreat) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value = value * sd(sea_ice_df$day_retreat) + mean(sea_ice_df$day_retreat),
         var = as.numeric(as.character(var)))


(main_plot_day_retreat <- ggplot() +
    geom_line(data = sea_ice_df,
              aes(x = year, y = day_retreat)) +
    geom_point(data = sea_ice_df,
              aes(x = year, y = day_retreat), size = 1) +
    geom_line(data = df_plot_day_retreat_2, aes(x = var, y = value, group = iteration), color = "grey35", linewidth = 0.3, alpha = 0.15) +
    geom_line(data = df_plot_day_retreat_1, aes(x = var, y = median), linetype = "solid") +
    scale_y_continuous(breaks = yday(c("1991-03-01", "1991-04-01", "1991-05-01", "1991-06-01", "1991-07-01", "1991-08-01")),
                       labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug")) +
    theme_bw() +
    labs(x = "year", y = "date of sea-ice retreat")
)

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(slope = alpha_2*sd(sea_ice_df$day_retreat)/sd(sea_ice_df$year)) %>%
  dplyr::select(slope) %>%
  mutate(parameter = "slope")


caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = slope , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = slope, y = parameter, fill = after_stat(x < 0)), 
                   color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 
  # scale_y_continuous(breaks = c(-0.15, -0.1, -0.05, 0)) +

inset <- inset_element(
  caterpillar_plot,
  left = -0.07,
  bottom = -0.09,
  right = 0.48,
  top = 0.44
)

plot_day_retreat <- main_plot_day_retreat + inset


# Ice-free days
ice_free_days <- matrix(data = NA, nrow = dim(res)[1], ncol = length(time_s),
                      dimnames = list(1:dim(res)[1], sea_ice_df$year))
years <- sea_ice_df$year
for (t in 1:length(years)) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    ice_free_days[k, t] <- res[k, "beta[1]"] +
      res[k, "beta[2]"] * time_s[t]
  }                               
}

df_plot_ice_free_days_1 <- as.data.frame.table(ice_free_days) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  mutate(value = value * sd(sea_ice_df$ice_free_days) + mean(sea_ice_df$ice_free_days)) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))

n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

ice_free_days <- matrix(data = NA, nrow = n_iterations, ncol = length(time_s),
                      dimnames = list(1:n_iterations, sea_ice_df$year))
years <- sea_ice_df$year
for (t in 1:length(years)) {                 # year
  for (k in 1:n_iterations) {         # MCMC iteration
    ice_free_days[k, t] <- res[iterations[k], "beta[1]"] +
      res[iterations[k], "beta[2]"] * time_s[t]
  }                               
}

df_plot_ice_free_days_2 <- as.data.frame.table(ice_free_days) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value = value * sd(sea_ice_df$ice_free_days) + mean(sea_ice_df$ice_free_days),
         var = as.numeric(as.character(var)))

(main_plot_ice_free_days <- ggplot() +
    geom_line(data = sea_ice_df,
              aes(x = year, y = ice_free_days)) +
    geom_point(data = sea_ice_df,
               aes(x = year, y = ice_free_days), size = 1) +
    geom_line(data = df_plot_ice_free_days_2, aes(x = var, y = value, group = iteration), 
              color = "grey35", linewidth = 0.3, alpha = 0.15) +
    geom_line(data = df_plot_ice_free_days_1, aes(x = var, y = median), linetype = "solid") +
    theme_bw() +
    labs(x = "year", y = "number of ice-free days")
)

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(slope = beta_2*sd(sea_ice_df$ice_free_days)/sd(sea_ice_df$year)) %>%
  dplyr::select(slope) %>%
  mutate(parameter = "slope")

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = slope , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = slope, y = parameter, fill = after_stat(x > 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 


inset <- inset_element(
  caterpillar_plot,
  left = -0.07,
  bottom = 0.55,
  right = 0.48,
  top = 1.08
)

plot_ice_free_days <- main_plot_ice_free_days + inset

plot_sea_ice <- plot_day_retreat + plot_ice_free_days +
  plot_annotation(tag_levels = list(c("A", "", "B", ""))) 

ggsave(plot_sea_ice, filename = "02_outputs/Figure 3.png",
       width = 18, height = 7.5, units = "cm", dpi = 600)
  


# Fig. 4: departure VS sea-ice --------------------------------------------

# Sea ice data
sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F) 
sea_ice <- NULL
years <- 1988:2024
for (t in 1:length(years)) {
  sea_ice[t] <- sea_ice_data$ice_free_days[which(sea_ice_data$year == years[t])] 
}
sea_ice_s <- as.vector(scale(sea_ice))

# Denning phenology data
denning_dates <- read_csv("01_inputs/denning_dates.csv", show_col_types = F) %>%
  left_join(x = .,
            y = data.frame(year = 1989:2025,  # shifted by one year so that for den departure in year t, sea ice in t-1 is associated
                           ice_free_days = sea_ice),
            by = "year")


mean_date_out <- mean(denning_dates$doy_out, na.rm = T) 
sd_date_out <- sd(denning_dates$doy_out, na.rm = T)


load("02_outputs/fit_PA_CR_phenology_1.RData")
load("02_outputs/fit_PA_CR_phenology_2.RData")
res <- rbind(fit_PA_CR_phenology_1$samples, 
             fit_PA_CR_phenology_2$samples)

lengthgrid <- 100
grid <- seq(min(sea_ice), max(sea_ice), length = lengthgrid) 
grid_scaled <- seq(min(sea_ice_s), max(sea_ice_s), length = lengthgrid) 
den_in <- den_out <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                            dimnames = list(1:dim(res)[1], grid))
for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    den_in[k, t] <- res[k, "tau[1]"] +
      res[k, "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[k, "kappa[1]"] +
      res[k, "kappa[2]"] * den_in[k, t] +
      res[k, "kappa[4]"] * grid_scaled[t]
  }                               
}

df_plot_1 <- as.data.frame.table(den_out) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value_bt = value * sd_date_out + mean_date_out) %>%
  group_by(var) %>%
  summarize(median = median(value_bt)) %>%
  mutate(var = as.numeric(as.character(var)))

n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

den_in <- den_out <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                            dimnames = list(1:n_iterations, grid))
for (t in 1:lengthgrid) {                 # year
  for (k in 1:n_iterations) {         # MCMC iteration
    den_in[k, t] <- res[iterations[k], "tau[1]"] +
      res[iterations[k], "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[iterations[k], "kappa[1]"] +
      res[iterations[k], "kappa[2]"] * den_in[k, t] +
      res[iterations[k], "kappa[4]"] * grid_scaled[t]
  }                               
}

df_plot_2 <- as.data.frame.table(den_out) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value_bt = value * sd_date_out + mean_date_out) %>%
  mutate(var = as.numeric(as.character(var)))


denning_dates_plot <- denning_dates %>%
  mutate(offset = rnorm(n = nrow(denning_dates), mean = 0, sd = 1),
         ice_free_days = ice_free_days + offset) %>%
  pivot_longer(cols = c("doy_out_lower", "doy_out_upper"), 
               names_to = "boundary", values_to = "doy_out_interval") 

main_plot <- ggplot() +
  geom_point(data = denning_dates_plot,
             aes(x = ice_free_days, y = doy_out_interval), shape = 4, size = 1) +
  geom_line(data = denning_dates_plot,
            aes(x = ice_free_days, y = doy_out_interval, group = ID_year), linewidth = 0.25) +
  geom_jitter(data = denning_dates,
              aes(x = ice_free_days, y = doy_out), width = 1.5, size = 1) +
  geom_line(data = df_plot_2, aes(x = var, y = value_bt, group = iteration), color = "grey35", linewidth = 0.25, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  theme_bw() +
  scale_y_continuous(limits = c(355, NA),
                     breaks = 364 + yday(c("1991-01-01", "1991-02-01", "1991-03-01", "1991-04-01", "1991-05-01", "1991-06-01")),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun")) +
  theme_bw() +
  labs(x = "number of ice-free days (previous year)", y = "date of den departure")


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(slope = (tau_3 * kappa_2 + kappa_4)*sd_date_out/sd(sea_ice)) %>%
  dplyr::select(slope) %>%
  mutate(parameter = "slope")

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = slope , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = slope, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, linewidth = 0.4,
               position = position_dodge(0.5),
               orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.4, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  coord_cartesian(xlim = c(NA, 0)) +
  scale_x_continuous(breaks = c(-0.15, -0.1, -0.05, 0)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.03,
  bottom = -0.08,
  right = 0.35,
  top = 0.27
)

main_plot + inset

ggsave("02_outputs/Figure 4.png",
       width = 15, height = 9, units = "cm", dpi = 600)




# Fig. 5: early repr success VS departure --------------------------------------

denning_dates <- read_csv("01_inputs/denning_dates.csv", show_col_types = F)

mean_date_out <- mean(denning_dates$doy_out, na.rm = T) 
sd_date_out <- sd(denning_dates$doy_out, na.rm = T)

load("02_outputs/fit_PA_CR_phenology_1.RData")
load("02_outputs/fit_PA_CR_phenology_2.RData")
res <- rbind(fit_PA_CR_phenology_1$samples, 
             fit_PA_CR_phenology_2$samples)

lengthgrid <- 100
grid <- seq(min(denning_dates$doy_out, na.rm = T),
            max(denning_dates$doy_out, na.rm = T),
            length = lengthgrid)
grid_scaled <- (grid - mean_date_out)/sd_date_out


# Early litter survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

probability <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                      dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:dim(res)[1]) {                # MCMC iteration
    probability[k, t] <- plogis(res[k, "beta_beta[2]"] +
                                  res[k, "beta_beta[5]"] * grid_scaled[t])
  }                               
}

df_plot_1 <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))

n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

probability <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                      dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    probability[k, t] <- plogis(res[iterations[k], "beta_beta[2]"] +
                                  res[iterations[k], "beta_beta[5]"] * grid_scaled[t])
  }                               
}

df_plot_2 <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


set.seed(1)
offset <- rnorm(n = nrow(denning_dates), mean = 0, sd = 0.75)

plot_early_litter_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = denning_dates, aes(x = doy_out + offset), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  theme_bw() +
  scale_x_continuous(breaks = 365 + c(0, 32, 60, 91, 121),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May")) +
  labs(x = "date of departure", y = "early litter survival")


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  dplyr::select(slope = beta_beta_5) %>%
  mutate(parameter = "beta")

label_pd <- paste0("p[d] == ", round(get_probability_direction(res[, "beta_beta[5]"]), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = slope , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = slope, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5),
               orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(res[, "beta_beta[5]"]), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1, 5)) +
  # scale_x_continuous(breaks = c(-0.15, -0.1, -0.05, 0)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 


inset <- inset_element(
  caterpillar_plot,
  left = 0.43,
  bottom = -0.03,
  right = 0.98,
  top = 0.50
)

plot_early_litter_survival <- plot_early_litter_survival + inset


# Twinning probability +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

probability <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                      dimnames = list(1:dim(res)[1], grid))
for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:dim(res)[1]) {                # MCMC iteration
    probability[k, t] <- plogis(res[k, "beta_gamma[2]"] +
                                  res[k, "beta_gamma[6]"] * grid_scaled[t])
  }                               
}

df_plot_1 <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))

n_iterations <- 100
set.seed(2)
iterations <- sample(1:dim(res)[1], size = n_iterations)

probability <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                      dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    probability[k, t] <- plogis(res[iterations[k], "beta_gamma[2]"] +
                                  res[iterations[k], "beta_gamma[6]"] * grid_scaled[t])
  }                               
}

df_plot_2 <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))

set.seed(1)
offset <- rnorm(n = nrow(denning_dates), mean = 0, sd = 0.75)

plot_twinning <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = denning_dates, aes(x = doy_out + offset), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  theme_bw() +
  scale_x_continuous(breaks = 365 + c(0, 32, 60, 91, 121),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May")) +
  labs(x = "date of departure", y = "twinning probability")


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  dplyr::select(beta_gamma_6) %>%
  mutate(parameter = "beta")

label_pd <- paste0("p[d] == ", round(get_probability_direction(res[, "beta_gamma[6]"]), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = beta_gamma_6, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = beta_gamma_6, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(res[, "beta_gamma[6]"]), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1, 5)) +
  # scale_x_continuous(breaks = c(-0.15, -0.1, -0.05, 0)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 


inset <- inset_element(
  caterpillar_plot,
  left = 0.43,
  bottom = -0.03,
  right = 0.98,
  top = 0.50
)

plot_twinning <- plot_twinning + inset


(plot_early_litter_survival) + (plot_twinning) +
  plot_annotation(tag_levels = list(c("A", "", "B", ""))) 

ggsave("02_outputs/Figure 5.png",
       width = 18, height = 8, units = "cm", dpi = 600)



# Fig. 6: early repr success VS sea-ice -----------------------------------

denning_dates <- read_csv("01_inputs/denning_dates.csv", show_col_types = F)

mean_date_out <- mean(denning_dates$doy_out, na.rm = T) 
sd_date_out <- sd(denning_dates$doy_out, na.rm = T)

doy_out_s <- (denning_dates$doy_out - mean_date_out)/sd_date_out
mean_doy_out_s <- mean(doy_out_s, na.rm = T)

sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F) 
sea_ice <- NULL
years <- 1988:2024
for (t in 1:length(years)) {
  sea_ice[t] <- sea_ice_data$ice_free_days[which(sea_ice_data$year == years[t])] 
}
sea_ice_s <- as.vector(scale(sea_ice))


load("02_outputs/fit_PA_CR_phenology_1.RData")
load("02_outputs/fit_PA_CR_phenology_2.RData")
res <- rbind(fit_PA_CR_phenology_1$samples, 
             fit_PA_CR_phenology_2$samples)

lengthgrid <- 100
grid <- seq(min(sea_ice), max(sea_ice), length = lengthgrid) 
grid_scaled <- seq(min(sea_ice_s), max(sea_ice_s), length = lengthgrid) 



# ~~~ a. Indirect effect -------------------------------------------------------


# Early litter survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
         dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    den_in[k, t] <- res[k, "tau[1]"] +
      res[k, "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[k, "kappa[1]"] +
      res[k, "kappa[2]"] * den_in[k, t] +
      res[k, "kappa[4]"] * grid_scaled[t]
    
    early_litter_survival[k, t] <- plogis(res[k, "beta_beta[2]"] +
                                            res[k, "beta_beta[5]"] * den_out[k, t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
         dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    den_in[k, t] <- res[iterations[k], "tau[1]"] +
      res[iterations[k], "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[iterations[k], "kappa[1]"] +
      res[iterations[k], "kappa[2]"] * den_in[k, t] +
      res[iterations[k], "kappa[4]"] * grid_scaled[t]
    
    early_litter_survival[k, t] <- plogis(res[iterations[k], "beta_beta[2]"] +
                                            res[iterations[k], "beta_beta[5]"] * den_out[k, t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


set.seed(1)
ice_free_days <- sea_ice_data %>% filter(year > 1986)
offset <- rnorm(n = nrow(ice_free_days), mean = 0, sd = 0.75)

plot_early_litter_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(x = "ice-free days (previous year)", y = "early litter survival")


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_beta_5 * (kappa_2 * tau_3 + kappa_4)) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1.45, 0.45)) +
  scale_x_continuous(breaks = c(-1, -0.5, 0)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.05,
  bottom = -0.04,
  right = 0.45,
  top = 0.47
)

plot_early_litter_survival_indirect <- plot_early_litter_survival + inset


# Twinning prob ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
         dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    den_in[k, t] <- res[k, "tau[1]"] +
      res[k, "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[k, "kappa[1]"] +
      res[k, "kappa[2]"] * den_in[k, t] +
      res[k, "kappa[4]"] * grid_scaled[t]
    
    twinning[k, t] <- plogis(res[k, "beta_gamma[2]"] +
                               res[k, "beta_gamma[6]"] * den_out[k, t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(twinning) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
         dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    den_in[k, t] <- res[iterations[k], "tau[1]"] +
      res[iterations[k], "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[iterations[k], "kappa[1]"] +
      res[iterations[k], "kappa[2]"] * den_in[k, t] +
      res[iterations[k], "kappa[4]"] * grid_scaled[t]
    
    twinning[k, t] <- plogis(res[iterations[k], "beta_gamma[2]"] +
                               res[iterations[k], "beta_gamma[6]"] * den_out[k, t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(twinning) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))

plot_twinning <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  scale_size(range = c(0.2, 3)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days (previous year)", y = "twinning probability")

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_gamma_6 * (kappa_2 * tau_3 + kappa_4)) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice , y = parameter)) +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1.45, 0.45)) +
  scale_x_continuous(breaks = c(-1, -0.5, 0)) +
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")

inset <- inset_element(
  caterpillar_plot,
  left = -0.05,
  bottom = -0.04,
  right = 0.45,
  top = 0.47
)

plot_twinning_indirect <- plot_twinning + inset


# ~~~ b. All effects -------------------------------------------------------

# Early litter survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
         dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    den_in[k, t] <- res[k, "tau[1]"] +
      res[k, "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[k, "kappa[1]"] +
      res[k, "kappa[2]"] * den_in[k, t] +
      res[k, "kappa[4]"] * grid_scaled[t]
    
    early_litter_survival[k, t] <- plogis(res[k, "beta_beta[2]"] +
                                            res[k, "beta_beta[5]"] * den_out[k, t] +
                                            res[k, "beta_beta[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
         dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    den_in[k, t] <- res[iterations[k], "tau[1]"] +
      res[iterations[k], "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[iterations[k], "kappa[1]"] +
      res[iterations[k], "kappa[2]"] * den_in[k, t] +
      res[iterations[k], "kappa[4]"] * grid_scaled[t]
    
    early_litter_survival[k, t] <- plogis(res[iterations[k], "beta_beta[2]"] +
                                            res[iterations[k], "beta_beta[5]"] * den_out[k, t] +
                                            res[iterations[k], "beta_beta[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


set.seed(1)
ice_free_days <- sea_ice_data %>% filter(year > 1986)
offset <- rnorm(n = nrow(ice_free_days), mean = 0, sd = 0.75)

plot_early_litter_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(x = "ice-free days (previous year)", y = "early litter survival")


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_beta_6 + beta_beta_5 * (kappa_2 * tau_3 + kappa_4)) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice , y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = stat(x > 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1.45, 0.45)) +
  scale_x_continuous(breaks = c(-1, -0.5, 0)) +
  scale_fill_manual(values = c("grey75", "grey90")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")
#   coord_flip(ylim = c(-1.35, 0.35), clip = 'off')


inset <- inset_element(
  caterpillar_plot,
  left = -0.05,
  bottom = -0.04,
  right = 0.45,
  top = 0.47
)

plot_early_litter_survival_all <- plot_early_litter_survival + inset



# Twinning prob ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
         dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    den_in[k, t] <- res[k, "tau[1]"] +
      res[k, "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[k, "kappa[1]"] +
      res[k, "kappa[2]"] * den_in[k, t] +
      res[k, "kappa[4]"] * grid_scaled[t]
    
    twinning[k, t] <- plogis(res[k, "beta_gamma[2]"] +
                               res[k, "beta_gamma[6]"] * den_out[k, t] +
                               res[k, "beta_gamma[7]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(twinning) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

den_in <- den_out <- early_litter_survival <- twinning <- 
  matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
         dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    den_in[k, t] <- res[iterations[k], "tau[1]"] +
      res[iterations[k], "tau[3]"] * grid_scaled[t]
    
    den_out[k, t] <- res[iterations[k], "kappa[1]"] +
      res[iterations[k], "kappa[2]"] * den_in[k, t] +
      res[iterations[k], "kappa[4]"] * grid_scaled[t]
    
    twinning[k, t] <- plogis(res[iterations[k], "beta_gamma[2]"] +
                               res[iterations[k], "beta_gamma[6]"] * den_out[k, t] +
                               res[iterations[k], "beta_gamma[7]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(twinning) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))

plot_twinning <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days), alpha = 0.75, linewidth = 0.3) +
  ylim(c(0, 1)) +
  scale_size(range = c(0.2, 3)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days (previous year)", y = "twinning probability")

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_gamma_7 + beta_gamma_6 * (kappa_2 * tau_3 + kappa_4)) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice , y = parameter)) +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = stat(x < 0)), 
               color = NA, alpha = 1) + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.3, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5), orientation = "y") +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, label = label_pd, 
           parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(-1.45, 0.45)) +
  scale_x_continuous(breaks = c(-1, -0.5, 0)) +
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")

inset <- inset_element(
  caterpillar_plot,
  left = -0.05,
  bottom = -0.04,
  right = 0.45,
  top = 0.47
)

plot_twinning_all <- plot_twinning + inset

(plot_early_litter_survival_indirect + plot_twinning_indirect) /
  (plot_early_litter_survival_all + plot_twinning_all) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", "D", ""))) 

ggsave("02_outputs/Figure 6.png",
       width = 18, height = 15, units = "cm", dpi = 600)



# Fig. 7. late summer body condition  VS Sea-ice -------------------------------

load("02_outputs/fit_PA_girth.RData")
res <- rbind(fit_PA_girth$chain1,
             fit_PA_girth$chain2)

girth_data <- read_csv("01_inputs/girth_late_summer.csv", show_col_types = F) %>%
  left_join(x = ., 
            y = read_csv("01_inputs/sea_ice_metrics.csv", 
                         show_col_types = F),
            by = "year") %>%
  drop_na(day_retreat)
  
  
length_s <- as.vector(scale(girth_data$s_length))
doy_s <- as.vector(scale(girth_data$doy))

lengthgrid <- 100
grid <- seq(min(girth_data$day_retreat),
            max(girth_data$day_retreat),
            length = lengthgrid)
grid_scaled <- (grid - mean(girth_data$day_retreat))/sd(girth_data$day_retreat)

girth_s <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                  dimnames = list(1:dim(res)[1], grid))
for (l in 1:lengthgrid) {         # iteration
  for (k in 1:dim(res)[1]) {          # MCMC iteration
    girth_s[k, l] <- res[k, "gamma[1]"] +
      res[k, "gamma[4]"] * grid_scaled[l] 
  }                               
}

df_plot_1 <- as.data.frame.table(girth_s) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value = value * sd(girth_data$girth) + mean(girth_data$girth)) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

girth_s <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                  dimnames = list(1:n_iterations, grid))
for (l in 1:lengthgrid) {              # date of departure
  for (k in 1:n_iterations) {                # MCMC iteration
    girth_s[k, l] <- res[iterations[k], "gamma[1]"] +
      res[iterations[k], "gamma[4]"] * grid_scaled[l] 
  }                               
}

df_plot_2 <- as.data.frame.table(girth_s) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(value = value * sd(girth_data$girth) + mean(girth_data$girth),
         var = as.numeric(as.character(var)))


# Length- and date-corrected mass
girth_s <- as.vector(scale(girth_data$girth))
residuals_girth  <- matrix(data = NA, nrow = dim(res)[1], ncol = nrow(girth_data), 
                           dimnames = list(1:dim(res)[1], 1:nrow(girth_data))) 
for (k in 1:nrow(girth_data)) { 
  for (l in 1:dim(res)[1]) { 
    residuals_girth[l, k] <- girth_s[k] - 
      
      res[l, "gamma[2]"] * res[l, paste0("length_s[", k,"]")]  -
      res[l, "gamma[3]"] * doy_s[k]
  } 
} 

residuals_girth_df <- as.data.frame.table(residuals_girth) %>% 
  rename(iteration = Var1, individual = Var2, value = Freq) %>% 
  mutate(value = value * sd(girth_data$girth) + mean(girth_data$girth)) %>%
  group_by(individual) %>% 
  summarize(median = median(value), 
            ci_2.5 = quantile(value, probs = 0.025), 
            ci_97.5 = quantile(value, probs = 0.975)) %>% 
  mutate(sex = girth_data$sex,
         day_retreat = girth_data$day_retreat) 



plot_body_condition <- ggplot() +
  geom_jitter(data = residuals_girth_df,
              aes(x = day_retreat, y = median)) +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), color = "grey35", linewidth = 0.3, alpha = 0.15) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  scale_x_continuous(breaks = yday(ymd(c("2021-05-01", "2021-06_01", "2021-07-01",
                                         "2021-08-01"))),
                     labels = c("May", "Jun", "July", "Aug")) +
  theme_bw() +
  theme( panel.background = element_rect(fill = "transparent", colour = NA)) +
  ylim(c(NA, 170)) +
  labs(x = "date of sea-ice retreat", y = "girth in late-summer")

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(slope = gamma_4*sd(girth_data$girth)/sd(girth_data$day_retreat)) %>%
  dplyr::select(slope) %>%
  mutate(parameter = "slope")
quantile(caterpillar$slope, prob = c(0.025, 0.5, 0.975))

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$slope), 2))


caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = slope, y = parameter)) +
  stat_halfeye(aes(x = slope, y = parameter, fill = stat(x > 0)), 
               color = NA, alpha = 1) + 
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, linewidth = 0.4,
               position = position_dodge(0.5), orientation = "y") +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.4, linewidth = 1, 
               position = position_dodge(0.5),  orientation = "y") +
  annotate("text", y = 1.3, x = median(caterpillar$slope), 
           label = label_pd, parse = TRUE, size = 3.25) +
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  coord_cartesian(xlim = c(-0.075, 0.275)) +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.03,
  bottom = 0.55,
  right = 0.40,
  top = 1.01
)

plot_body_condition + inset

ggsave("02_outputs/Figure 7.png",
       width = 11, height = 8, units = "cm", dpi = 600)


#___________________________________________________________________________----










