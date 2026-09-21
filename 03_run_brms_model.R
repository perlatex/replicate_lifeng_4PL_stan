library(tidyverse)
library(brms)
library(posterior)
library(tidybayes)


raw_data <- read_csv("data/fig4f.csv") |>
  mutate(compound = fct_relevel(compound, "ent-Paroxol", "3-Fluoro-5-hydroxybenzonitrile"))



summary_data <- raw_data |>
  filter(conc_um > 0) |>
  group_by(compound, conc_um) |>
  summarise(
    mean = mean(viability), 
    sd   = sd(viability), 
    .groups = "drop"
  ) |>
  mutate(ymin = mean - sd, ymax = mean + sd)




# 非线性公式: IC50 放在对数尺度 (exp(logic50)) 保证为正;
# 每个参数 ~ 0 + compound 让两化合物各自有一套系数. --------------------------
# sigma ~ 0 + compound: 每个化合物各自的残差 sd (log 链接), 与 Stan 版一致.

bform <- bf(
  viability ~ bottom + (top - bottom) / (1 + (conc_um / exp(logic50))^hill),
  bottom + top + logic50 + hill ~ 0 + compound,
  sigma ~ 0 + compound,
  nl = TRUE
)

bprior <- c(
  prior(normal(0, 0.1),   nlpar = "bottom"),
  prior(normal(1, 0.15),  nlpar = "top"),
  prior(normal(6.21, 1),  nlpar = "logic50"),   # log(500) ≈ 6.21
  prior(normal(1.5, 0.8), nlpar = "hill", lb = 0)
)



fit <- brm(bform,
           data       = raw_data,
           prior      = bprior,
           backend    = "cmdstanr",
           chains     = 4,
           cores      = 4,
           iter       = 2000,
           seed       = 1024,
           control    = list(adapt_delta = 0.98),
           file       = "fits/brms_fit",
           file_refit = "on_change")          # 公式或先验改动后自动重拟合


fit

fit |> 
  as_draws_df(variable = "^b_(bottom|top|logic50|hill)_", regex = TRUE) |>
  as_tibble() |>
  pivot_longer(starts_with("b_"), names_to = "variable") |>
  mutate(par = str_match(variable, "^b_([a-z0-9]+)_")[, 2]) |>
  group_by(par) |>
  mutate(compound = levels(raw_data$compound)[match(variable, unique(variable))]) |>
  ungroup() |>
  mutate(value = if_else(par == "logic50", exp(value), value),
         par   = if_else(par == "logic50", "ic50", par)) |>
  group_by(compound, par) |>
  median_qi(value) |>
  select(compound, par, median = value, .lower, .upper)



pal <- c("ent-Paroxol"                    = "#4A98E0",
         "3-Fluoro-5-hydroxybenzonitrile" = "#F47C7C")

# 后验拟合曲线: 每个化合物在浓度网格上取 add_epred_draws 的中位数 + 95% 区间.
# newdata 必须含模型两个预测变量 (conc_um, compound). ------------------------
curve_data <- expand_grid(
  conc_um  = 10^seq(log10(5), log10(2000), length.out = 300),
  compound = factor(levels(raw_data$compound), levels = levels(raw_data$compound))
) |>
  add_epred_draws(fit) |>
  median_qi(.epred, .width = 0.95) |>
  ungroup()


curve_data |> 
  ggplot(aes(colour = compound)) +
  geom_ribbon(
    aes(x = conc_um, ymin = .lower, ymax = .upper, fill = compound),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    aes(x = conc_um, y = .epred),
    linewidth = 0.8
  ) +
  geom_errorbar(
    data = summary_data,
    aes(x = conc_um, ymin = mean - sd, ymax = mean + sd),
    width = 0.06, linewidth = 0.6
  ) +
  geom_point(
    data = summary_data,
    aes(x = conc_um, y = mean),
    size = 2.6
  ) +
  scale_x_log10(breaks = c(10, 100, 1000),
                labels = scales::label_math(10^.x)(c(1, 2, 3)),
                limits = c(4, 2200)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  coord_cartesian(xlim = c(4, 2200), ylim = c(0, 1.2)) + # 只裁剪视图
  labs(x = "Concentration/\u00b5M", y = "Cell Viability",
       colour = NULL, fill = NULL) +
  theme_classic(base_size = 15) +
  theme(axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        axis.ticks.length = unit(4, "pt"),
        legend.position = c(0.02, 0.14),
        legend.justification = c(0, 0.5),
        legend.text = element_text(size = 11))
