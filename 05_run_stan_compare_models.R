# 04_compare_models.R ---------------------------------------------------------
# 比较四个版本的 4PL 拟合曲线: nls, brms (wide), stan_tight, stan_wide
#   - 实线: 点估计 (nls 最小二乘 / 贝叶斯 MAP)
#   - 虚线: 贝叶斯逐点后验中位数曲线
# 第三个面板为两化合物曲线之差 (3-F − ent-P), 直接比较"间距".
# ----------------------------------------------------------------------------

library(tidyverse)
library(brms)
library(cmdstanr)
library(posterior)
library(tidybayes)


compound_levels <- c("ent-Paroxol", "3-Fluoro-5-hydroxybenzonitrile")

raw_data <- read_csv("data/fig4f.csv", show_col_types = FALSE) |>
  mutate(compound    = factor(compound, levels = compound_levels),
         compound_id = as.integer(compound))

summary_data <- raw_data |>
  filter(conc_um > 0) |>
  group_by(compound, conc_um) |>
  summarise(mean = mean(viability), sd = sd(viability), .groups = "drop")

grid <- crossing(
  compound_id = seq_along(compound_levels),
  conc_um     = 10^seq(log10(5), log10(2000), length.out = 300)
) |>
  mutate(compound = factor(compound_levels[compound_id], levels = compound_levels),
         grid_id  = row_number())

f4 <- function(x, bottom, top, ec50, hill) bottom + (top - bottom) / (1 + (x / ec50)^hill)

# 由参数表 (每化合物一行: bottom, top, ec50, hill) 生成曲线
curve_from_pars <- function(pars) {
  grid |>
    left_join(pars, by = join_by(compound)) |>
    mutate(value = f4(conc_um, bottom, top, ec50, hill)) |>
    select(compound, conc_um, value)
}


# 1. nls (最小二乘) ------------------------------------------------
# 01 脚本用 minpack::nlsLM; 这里用 base nls 的 port 算法, 估计值相同.
fit_nls <- function(data) {
  nls(viability ~ bottom + (top - bottom) / (1 + (conc_um / ec50)^hill),
      data      = data,
      start     = list(bottom = -0.2, top = 1, ec50 = 500, hill = 1),
      lower     = c(bottom = -5, top = 0, ec50 = 1, hill = 0.05),
      upper     = c(bottom = 1, top = 3, ec50 = 1e5, hill = 20),
      algorithm = "port",
      control   = list(maxiter = 1000))
}

nls_pars <- raw_data |>
  group_nest(compound) |>
  mutate(coef = map(data, ~ as_tibble_row(coef(fit_nls(.x))))) |>
  select(compound, coef) |>
  unnest(coef)

nls_curve <- curve_from_pars(nls_pars) |>
  mutate(version = "nls", summary = "point (nls / MAP)")


# 2. Stan: tight vs wide 先验 --------------------------------------------------
stan_model_4pl <- cmdstan_model("Stan/dose_response_compare.stan")

stan_base_data <- list(
  N                  = nrow(raw_data),
  J                  = length(compound_levels),
  compound           = raw_data$compound_id,
  concentration      = raw_data$conc_um,
  y                  = raw_data$viability,
  M                  = nrow(grid),
  compound_grid      = grid$compound_id,
  concentration_grid = grid$conc_um
)

stan_priors <- list(
  stan_tight = list(prior_bottom_sd = 0.2, 
                    prior_gap_sd = 0.25,
                    prior_hill_mu = log(1.5), prior_hill_sd = 0.8),
  stan_wide  = list(prior_bottom_sd = 1,   
                    prior_gap_sd = 1,
                    prior_hill_mu = 0,        prior_hill_sd = 1.5)
)

stan_init <- list(bottom   = rep(-0.1, 2),
                  log_gap  = rep(0, 2),
                  log_ec50 = log(c(500, 500)),
                  log_hill = log(c(1.5, 1.5)),
                  sigma    = rep(0.1, 2))

run_stan <- function(priors, version) {
  data <- c(stan_base_data, priors)

  fit <- stan_model_4pl$sample(
    data = data, 
    seed = 1024, 
    chains = 4,
    parallel_chains = 4,
    init = rep(list(stan_init), 4),
    iter_warmup = 2000, 
    iter_sampling = 2000,
    adapt_delta = 0.99, 
    max_treedepth = 12, 
    refresh = 0
  )

  fit_map <- stan_model_4pl$optimize(
    data = data, seed = 1024, init = list(stan_init),
    jacobian = FALSE, algorithm = "lbfgs", refresh = 0
  )

  map_pars <- fit_map$summary(c("bottom", "top", "ec50", "hill")) |>
    separate_wider_regex(variable, c(par = "\\w+", "\\[", compound_id = "\\d+", "\\]")) |>
    mutate(compound = factor(compound_levels[as.integer(compound_id)], levels = compound_levels)) |>
    select(-compound_id) |>
    pivot_wider(names_from = par, values_from = estimate)

  median_curve <- fit |>
    spread_draws(mu_grid[grid_id]) |>
    group_by(grid_id) |>
    summarise(value = median(mu_grid), .groups = "drop") |>
    left_join(grid, by = join_by(grid_id)) |>
    select(compound, conc_um, value)

  list(
    pars  = map_pars,
    curve = bind_rows(
      curve_from_pars(map_pars) |> mutate(summary = "point (nls / MAP)"),
      median_curve              |> mutate(summary = "posterior median")
    ) |>
      mutate(version = version)
  )
}

stan_res <- imap(stan_priors, run_stan) # map2(x, names(x), ...) 


# 3. brms (wide 先验, 与 02 脚本一致; 有缓存则直接读取) -------------------------
bform <- bf(
  viability ~ bottom + (top - bottom) / (1 + (conc_um / exp(logic50))^hill),
  bottom + top + logic50 + hill ~ 0 + compound,
  sigma ~ 0 + compound,
  nl = TRUE
)

bprior <- c(
  prior(normal(0, 1),    nlpar = "bottom"),
  prior(normal(1, 0.5),  nlpar = "top"),
  prior(normal(6.21, 1), nlpar = "logic50"),
  prior(normal(1, 5),    nlpar = "hill", lb = 0)
)

dir.create("fits", showWarnings = FALSE)

fit_brms <- brm(bform,
                data       = raw_data,
                prior      = bprior,
                backend    = "cmdstanr",
                chains     = 4,
                cores      = 4,
                iter       = 2000,
                seed       = 1024,
                control    = list(adapt_delta = 0.98),
                file       = "fits/brms_fit",
                file_refit = "on_change")

brms_map <- cmdstan_model(write_stan_file(stancode(fit_brms)))$optimize(
  data      = standata(fit_brms),
  seed      = 1024,
  jacobian  = FALSE,
  algorithm = "lbfgs",
  refresh   = 0,
  init      = list(list(b_bottom  = c(-0.1, -0.1),
                        b_top     = c(0.95, 0.95),
                        b_logic50 = log(c(500, 500)),
                        b_hill    = c(1.5, 1.5),
                        b_sigma   = log(c(0.1, 0.1))))
)

# b_<nlpar>[k] 的 k 与 compound 因子水平顺序一致
brms_pars <- brms_map$summary(c("b_bottom", "b_top", "b_logic50", "b_hill")) |>
  separate_wider_regex(variable, c("b_", par = "\\w+", "\\[", k = "\\d+", "\\]")) |>
  mutate(compound = factor(compound_levels[as.integer(k)], levels = compound_levels)) |>
  select(-k) |>
  pivot_wider(names_from = par, values_from = estimate) |>
  mutate(ec50 = exp(logic50), .keep = "unused")

brms_median <- grid |>
  select(compound, conc_um) |>
  add_epred_draws(fit_brms) |>
  summarise(value = median(.epred), .groups = "drop") |>
  select(compound, conc_um, value)

brms_curve <- bind_rows(
  curve_from_pars(brms_pars) |> mutate(summary = "point (nls / MAP)"),
  brms_median                |> mutate(summary = "posterior median")
) |>
  mutate(version = "brms")


# 4. 汇总: 点估计参数表 --------------------------------------------------------
version_levels <- c("nls", "brms", "stan_wide", "stan_tight")

param_table <- bind_rows(
  nls        = nls_pars,
  brms       = brms_pars,
  stan_wide  = stan_res$stan_wide$pars,
  stan_tight = stan_res$stan_tight$pars,
  .id = "version"
) |>
  mutate(version = factor(version, levels = version_levels)) |>
  select(compound, version, bottom, top, ec50, hill) |>
  arrange(compound, version)

print(param_table, n = Inf)


# 5. 曲线与间距 ----------------------------------------------------------------
curve_all <- bind_rows(nls_curve, brms_curve,
                       map(stan_res, "curve") |> list_rbind()) |>
  mutate(version = factor(version, levels = version_levels))

gap_curve <- curve_all |>
  pivot_wider(names_from = compound, values_from = value) |>
  mutate(value = `3-Fluoro-5-hydroxybenzonitrile` - `ent-Paroxol`,
         panel = "Gap: 3-F − ent-P") |>
  select(-all_of(compound_levels))

panel_levels <- c(compound_levels, "Gap: 3-F − ent-P")

plot_curves <- bind_rows(
  curve_all |> mutate(panel = as.character(compound)) |> select(-compound),
  gap_curve
) |>
  mutate(panel = factor(panel, levels = panel_levels))

plot_points <- summary_data |>
  mutate(panel = factor(as.character(compound), levels = panel_levels))

# 选定浓度下的间距数值表 (在 log 浓度上线性插值)
gap_at <- c(50, 100, 200, 300, 500)

gap_table <- gap_curve |>
  group_by(version, summary) |>
  reframe(value   = approx(log(conc_um), value, xout = log(gap_at))$y,
          conc_um = gap_at) |>
  pivot_wider(names_from = conc_um, names_prefix = "gap_", values_from = value) |>
  arrange(summary, version)

print(gap_table, width = Inf)


# 6. 出图 ----------------------------------------------------------------------
version_pal <- c(nls        = "grey70",
                 brms       = "#E69F00",
                 stan_wide  = "#0072B2",
                 stan_tight = "#CC3311")

ggplot(plot_curves, aes(conc_um, value, colour = version, linetype = summary)) +
  geom_hline(data = tibble(panel = factor(panel_levels[3], levels = panel_levels)),
             aes(yintercept = 0), colour = "grey70", inherit.aes = FALSE) +
  geom_errorbar(data = plot_points,
                aes(x = conc_um, ymin = mean - sd, ymax = mean + sd),
                inherit.aes = FALSE, width = 0.05, colour = "grey55") +
  geom_point(data = plot_points, aes(conc_um, mean),
             inherit.aes = FALSE, colour = "grey35", size = 2) +
  # nls 画成粗灰底线: 与 brms / stan_wide 的 MAP 重合时仍可见
  geom_line(aes(linewidth = version)) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_x_log10(breaks = c(10, 100, 1000), labels = scales::label_log()) +
  scale_colour_manual(values = version_pal) +
  scale_linewidth_manual(values = c(nls = 2.4, brms = 0.8, stan_wide = 0.8, stan_tight = 0.8),
                         guide = "none") +
  scale_linetype_manual(values = c("point (nls / MAP)" = "solid",
                                   "posterior median"  = "22")) +
  labs(x = "Concentration/µM", y = "Cell Viability  (gap panel: difference)",
       colour = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
