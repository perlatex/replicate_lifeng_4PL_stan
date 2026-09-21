# run_ggplot2.R --------------------------------------------------------------
# 复现 Nature Communications (2024) 15:8396 图 4f
# Hep G2 细胞在两种化合物下的剂量-反应曲线 (CCK-8, mean ± sd, n = 6)
# 用 nls 拟合 4PL 曲线, 用 ggplot2 出图
# ----------------------------------------------------------------------------

library(tidyverse)
library(minpack.lm)


# 1. 读入逐孔数据, 化合物设为有序因子 -----------------------------------------
raw_data <- read_csv("data/fig4f.csv") |>
  mutate(compound = fct_relevel(compound, "ent-Paroxol", "3-Fluoro-5-hydroxybenzonitrile"))



# 2. 汇总 mean ± sd; 0 uM 对照无法上对数轴, 绘图时滤除 -------------------------
summary_data <- raw_data |>
  filter(conc_um > 0) |>
  group_by(compound, conc_um) |>
  summarise(
    mean = mean(viability), 
    sd   = sd(viability), 
    .groups = "drop"
  ) |>
  mutate(ymin = mean - sd, ymax = mean + sd)


# 3. 逐化合物拟合 4PL, 生成平滑曲线 -------------------------------------------
fit_4pl <- function(data) {
  nlsLM(
    viability ~ bottom + (top - bottom) / (1 + (conc_um / ec50)^hill),
    data = data, 
    start = list(bottom = -0.2, top = 1, ec50 = 500, hill = 1),
    control = nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10)
  )
}


grid <- tibble(conc_um = 10^seq(log10(5), log10(2000), length.out = 600))

pred_data <- raw_data |>
  group_nest(compound) |>
  mutate(
    mod = map(data, fit_4pl)
  ) |> 
  mutate(pred = map(mod, ~ mutate(grid, mean = predict(.x, newdata = grid)))) |>
  select(compound, pred) |>
  unnest(pred)



# 4. 出图 --------------------------------------------------------------------
summary_data |> 
  ggplot(aes(x = conc_um, y = mean, colour = compound)) +
  geom_line(data = pred_data, linewidth = 0.8) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.06, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_x_log10(
    breaks = c(10, 100, 1000),
    labels = scales::label_log()
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.2, 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_colour_manual(
    values = c("ent-Paroxol"                    = "#4A98E0",
               "3-Fluoro-5-hydroxybenzonitrile" = "#F47C7C")
  ) +
  labs(x = "Concentration/\u00b5M", y = "Cell Viability", colour = NULL) +
  theme_test(base_size = 15) +
  theme(axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        axis.ticks.length = unit(4, "pt"),
        legend.position = c(0.02, 0.14),
        legend.justification = c(0, 0.5),
        legend.text = element_text(size = 11))

