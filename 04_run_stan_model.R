library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(bayesplot)


###############################################################################
compound_levels <- c("ent-Paroxol", "3-Fluoro-5-hydroxybenzonitrile")

pal <- c("ent-Paroxol"                    = "#4A98E0",
         "3-Fluoro-5-hydroxybenzonitrile" = "#F47C7C")


raw_data <- read_csv("data/fig4f.csv")

summary_data <- raw_data |>
  filter(conc_um > 0) |>  
  group_by(compound, conc_um) |>
  summarise(
    mean = mean(viability), 
    sd   = sd(viability), 
    .groups = "drop"
  ) |>
  mutate(ymin = mean - sd, ymax = mean + sd)
###############################################################################



###############################################################################
# ---- 1. Prepare data ----
d <- raw_data |>
  mutate(
    compound    = factor(compound, levels = c("ent-Paroxol", "3-Fluoro-5-hydroxybenzonitrile")),
    compound_id = as.integer(compound)
  ) |>
  arrange(compound_id, conc_um, replicate) |>
  mutate(obs_id = row_number())



grid <- crossing(
  compound_id = seq_along(compound_levels),
  conc_um = 10^seq(log10(5), log10(2000), length.out = 600)
) |>
  mutate(
    compound = factor(compound_levels[compound_id], levels = compound_levels),
    grid_id  = row_number()
  )

# 进入模型的是原始数据，包括conc_um =0，
# 而图中的点是不是进入模型的数据，只是汇总后的数据，它删除了0
stan_data <- list(
  N                  = nrow(d),
  J                  = length(compound_levels),
  compound           = d$compound_id,
  concentration      = d$conc_um,
  y                  = d$viability,
  M                  = nrow(grid),
  compound_grid      = grid$compound_id,
  concentration_grid = grid$conc_um
)
###############################################################################



###############################################################################
# ---- 2. Compile and sample ----
model <- cmdstan_model("Stan/dose_response1.stan", pedantic = TRUE)
# model <- cmdstan_model("Stan/dose_response2.stan", pedantic = TRUE)


# 每条链用不同的 log_hill 初值以帮助探索.
init <- map(seq_len(4), function(k) {
  list(
    bottom   = rep(-0.2, 2),
    log_gap  = rep(log(0.9), 2),      # top = bottom + exp(log_gap) ≈ 0.92
    log_ec50 = log(c(400, 600)),
    log_hill = log(rep(c(1, 1.2, 2, 4)[k], 2)),
    sigma    = rep(0.2, 2)
  )
})

fit <- model$sample(
  data            = stan_data,
  seed            = 1024,
  chains          = 4,
  parallel_chains = 4,
  init            = init,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  adapt_delta     = 0.99,
  max_treedepth   = 12
)


fit$cmdstan_diagnose()

fit$draws(
  variables = c("ec50", "hill"),
  format = "draws_array"
) |>
  mcmc_trace()
###############################################################################




###############################################################################
# ---- 3. Parameter summary ----
parameter_summary <- fit$draws(
  variables = c("ec50", "hill", "top", "bottom", "sigma")
) |>
  summarise_draws(median, ~quantile(.x, c(0.025, 0.975)), rhat, ess_bulk)

parameter_summary
###############################################################################




###############################################################################
# ---- 4. Posterior curves with 90% credible band ----
# 用 grid_id 做 join 对齐, 不依赖行顺序

pred_df <- fit |>
  gather_draws(mu_grid[grid_id]) |>
  group_by(grid_id) |>
  summarise(
    mean = mean(.value),
    median = median(.value),
    q025 = quantile(.value, 0.025),
    q975 = quantile(.value, 0.975),
    .groups = "drop"
  ) |>
  left_join(grid, by = join_by(grid_id))


###############################################################################



###############################################################################
#  5. Reproduce Fig. 4f with posterior uncertainty
pred_df |> 
  ggplot(aes(x = conc_um, y = median, colour = compound)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = compound), color = NA, alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  geom_errorbar(
    data = summary_data,
    aes(y = mean, ymin = mean - sd, ymax = mean + sd),
    width = 0.055,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = summary_data,
    aes(y = mean),
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_log10(
    breaks = c(10, 100, 1000),
    labels = scales::label_log()
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.2, 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  coord_cartesian(
    xlim = c(3, 3000),
    ylim = c(0, 1.2)
  ) +
  labs(x = "Concentration/\u00b5M", y = "Cell Viability",
       colour = NULL, fill = NULL) +
  theme_test(base_size = 15) +
  theme(axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        axis.ticks.length = unit(4, "pt"),
        legend.position = c(0.02, 0.14),
        legend.justification = c(0, 0.5),
        legend.text = element_text(size = 11))
###############################################################################

