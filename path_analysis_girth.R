library(tidyverse)
library(nimble)


# ~ 1. Build and run model ------------------------------------------------------

model_PA_girth <- nimbleCode ({
  
  for (t in 1:T) {
    # sea-ice: day retreat
    mu_retreat[t] <- beta[1] + 
      beta[2] * time_s[t]
    
    sea_ice_s[t] ~ dnorm(mu_retreat[t], sd = sigma_retreat)
  }
  
  for (k in 1:N) {
    
    # doy
    mu_doy[k] <- alpha[1] + 
      alpha[2] * time_s[index_year[k]]
    
    doy_s[k] ~ dnorm(mu_doy[k], sd = sigma_doy)
    
    
    
    # Length
    mu_length[k] <- delta[1] +
      delta[2] * age_s[k] 
    
    length_s[k] ~ dnorm(mu_length[k], sd = sigma_length)
    
    # girth
    mu_girth[k] <- gamma[1] + 
      gamma[2] * length_s[k] +
      gamma[3] * doy_s[k] +
      gamma[4] * sea_ice_s[index_year[k]] 
    
    girth_s[k] ~ dnorm(mu_girth[k], sd = sigma_girth)
  }
  
  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++
  
  for (k in 1:2) {
    alpha[k] ~ dnorm(0, sd = 3)  # coefs for doy
    beta[k] ~ dnorm(0, sd = 3)  # coefs for sea-ice
    delta[k] ~ dnorm(0, sd = 3)  # coefs for length
  }
  sigma_doy ~ dunif(0, 5)
  sigma_retreat ~ dunif(0, 5)
  sigma_length ~ dunif(0, 5)
  
  for (k in 1:4) {
    gamma[k] ~ dnorm(0, sd = 3)  # coefs for girth
  }
  sigma_girth ~ dunif(0, 5)
})




girth_data <- read_csv("01_inputs/girth_late_summer.csv", show_col_types = F) %>%
  left_join(x = ., 
            y = read_csv("01_inputs/sea_ice_metrics.csv", 
                         show_col_types = F),
            by = "year") %>%
  drop_na(day_retreat)

env_data <- girth_data %>%
  distinct(year, day_retreat)

sea_ice_s <- as.vector(scale(env_data$day_retreat))
time_s <- as.vector(scale(env_data$year))


index_year <- NULL
for (k in 1:nrow(girth_data)) {
  index_year[k] <- which(env_data$year == girth_data$year[k])
}

age_s <- as.vector(scale(girth_data$age_for_analyses))
length_s <- as.vector(scale(girth_data$s_length))
doy_s <- as.vector(scale(girth_data$doy))
girth_s <- as.vector(scale(girth_data$girth))

table(girth_data$year)


dat <- list(doy_s = doy_s,
            sea_ice_s = sea_ice_s,
            length_s = length_s,
            girth_s = girth_s)

my.constants <- list(T = length(time_s),
                     N = length(age_s),
                     time_s = time_s,
                     age_s = age_s,
                     index_year = index_year)

inits <- function() list(alpha = rnorm(n = 2, mean = 0, sd = 0.5),
                         sigma_doy = runif(1, min = 0, max = 1),
                         beta = rnorm(n = 2, mean = 0, sd = 0.5),
                         sigma_retreat = runif(1, min = 0, max = 1),
                         delta = rnorm(n = 2, mean = 0, sd = 0.5),
                         sigma_length = runif(1, min = 0, max = 1),
                         gamma = rnorm(n = 4, mean = 0, sd = 0.5),
                         sigma_girth = runif(1, min = 0, max = 1))


# Parameters monitored 
params <- c("alpha", "sigma_doy", "beta", "sigma_retreat",
            "delta", "sigma_length", "gamma", "sigma_girth",
            "length_s")





fit_PA_girth <- nimbleMCMC(code = model_PA_girth,
                           constants = my.constants,
                           data = dat,
                           inits = inits,
                           monitors = params,
                           niter = 40000,
                           nburnin = 10000,
                           thin = 10,
                           nchains = 2)
save(fit_PA_girth,
     file = "02_outputs/fit_PA_girth.RData")


# ~ 2. Check results ------------------------------------------------------------

check_convergence <- function(params.plot, nimble_output) {
  # Process Nimble output into dataframe
  n_chains <- length(nimble_output)
  chains <- data.frame(do.call(rbind, nimble_output)) %>%
    dplyr::select(any_of(params.plot)) %>%
    mutate(chain = rep(1:n_chains, each = nrow(nimble_output$chain1)),
           iteration = rep(1:nrow(nimble_output$chain1), times = n_chains))
  chains_l <- pivot_longer(chains, cols = any_of(params.plot), names_to = "parameter") %>%
    mutate(parameter = factor(parameter, levels = params.plot),
           chain = factor(chain, levels = c("1", "2", "3"))) 
  
  param.mean <- chains_l %>%
    group_by(parameter, chain) %>%
    summarise(m = mean(value))
  
  param.running.mean <- chains_l %>%
    arrange(parameter, iteration) %>%
    group_by(parameter, chain) %>%
    mutate(rm = cumsum(value)/iteration)
  
  trace.plots <- ggplot(data = chains_l, 
                        aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.5) +
    labs(y = "trace") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free",
                ncol = 1)
  
  density.plots <- ggplot(data = chains_l, 
                          aes(x = value, color = chain, fill = chain)) +
    geom_vline(xintercept = 0) +
    geom_density(alpha = 0.25) +
    
    labs(x = "density") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap( ~ parameter,
                scales = "free_y",
                ncol = 1)
  
  running.mean.plot <- ggplot(param.running.mean, 
                              aes(x = iteration, y = rm, color = chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept = m), param.mean,
               colour = "black", alpha = 0.5) + 
    theme_bw() +
    theme(legend.position = "none") +
    ylab("running mean") +
    facet_grid(parameter ~ ., scales = "free")
  
  # Plot all the plots together
  diagnostic_plot <- trace.plots + density.plots + running.mean.plot
  
  return(diagnostic_plot)
}


load("02_outputs/fit_PA_girth.RData")
params.plot <- c()
for (k in 1:2) {
  params.plot <- c(params.plot, paste0("alpha.", k, "."))
}
params.plot <- c(params.plot, "sigma_doy")
for (k in 1:2) {
  params.plot <- c(params.plot, paste0("beta.", k, "."))
}
params.plot <- c(params.plot, "sigma_retreat")
for (k in 1:2) {
  params.plot <- c(params.plot, paste0("delta.", k, "."))
}
params.plot <- c(params.plot, "sigma_length")
for (k in 1:4) {
  params.plot <- c(params.plot, paste0("gamma.", k, "."))
}
params.plot <- c(params.plot, "sigma_girth")


diagnostic_plot <-  check_convergence(nimble_output = fit_PA_girth,
                                      params.plot = params.plot)
diagnostic_plot

