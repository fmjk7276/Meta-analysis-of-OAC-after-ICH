library(rstan)
library(ggplot2)
library(posterior)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -------------------------
# data
# -------------------------
data <- read.csv("data.csv")

trial <- data$trial

treatment_events <- data$treatment_events
treatment_total  <- data$treatment_total
control_events   <- data$control_events
control_total    <- data$control_total

# -------------------------
# model
# -------------------------
model <- stan_model("model.stan")

# -------------------------
# priors
# -------------------------
tau_priors <- list(
  Low = 0.2,
  High = 0.5
)

logRR_priors <- list(
  Neutral = c(0, 0.5),
  WeakFav = c(-0.3285, 0.5),
  StrongFav = c(-0.3285, 0.2),
  WeakUnfav = c(0.3285, 0.5),
  StrongUnfav = c(0.3285, 0.2)
)

if (!dir.exists("results")) dir.create("results")

fits <- list()

# -------------------------
# MCMC runs
# -------------------------
for (t in names(tau_priors)) {
  for (p in names(logRR_priors)) {
    
    mu <- logRR_priors[[p]][1]
    sd <- logRR_priors[[p]][2]
    
    dat <- list(
      N = length(trial),
      treatment_events = treatment_events,
      treatment_total = treatment_total,
      control_events = control_events,
      control_total = control_total,
      logRR_prior_mu = mu,
      logRR_prior_sigma = sd,
      tau_prior_sigma = tau_priors[[t]]
    )
    
    name <- paste(t, p, sep = "_")
    
    fits[[name]] <- sampling(
      model,
      data = dat,
      iter = 4000,
      chains = 4,
      warmup = 2000,
      seed = 2026,
      control = list(adapt_delta = 0.95)
    )
  }
}

# -------------------------
# diagnostics
# -------------------------
diag <- data.frame()

for (n in names(fits)) {
  
  s <- summary(fits[[n]])$summary
  
  diag <- rbind(diag, data.frame(
    model = n,
    max_rhat = max(s[, "Rhat"], na.rm = TRUE),
    min_ess = min(s[, "n_eff"], na.rm = TRUE)
  ))
}

write.csv(diag, "results/diagnostics_table.csv", row.names = FALSE)

print(diag)

# -------------------------
# posterior summary (theta)
# -------------------------
out <- data.frame()

for (n in names(fits)) {
  
  theta <- extract(fits[[n]])$theta
  
  out <- rbind(out, data.frame(
    model = n,
    mean = mean(theta),
    low = quantile(theta, 0.025),
    high = quantile(theta, 0.975)
  ))
}

write.csv(out, "results/posterior_results_table.csv", row.names = FALSE)

# -------------------------
# simple plot
# -------------------------
df <- data.frame(theta = extract(fits[[1]])$theta)

p <- ggplot(df, aes(theta)) +
  geom_density(fill = "red", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal()

ggsave("results/posterior.png", p)