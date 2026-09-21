library(tidyverse)

rawdata <- read_csv("./data/fig4f.csv")
rawdata

# 按照（化合物 × 浓度）分组
d <- rawdata |>
  filter(conc_um > 0) |>
  summarise(
    mean = mean(viability),
    sd   = sd(viability),
    .by  = c(compound, conc_um)
  ) 



# 建模的是 rawdata, 画图的是 d
rawdata |> 
  ggplot(aes(x = conc_um, y = viability, color = compound)) +
  geom_point(data = d, aes(x = conc_um, y = mean)) +
  geom_smooth(
    method = "nls",
    formula = y ~ bottom + (top - bottom) / (1 + (x / ec50)^hill),
    method.args = list(
      start     = c(bottom = -0.2, top = 1, ec50 = 500, hill = 1),
      algorithm = "port",
      control   = list(maxiter = 1000)
    ),
    se = FALSE,
    linewidth = 0.9
  ) 
  

# 下一步：x 轴放在 log10 尺度
# 因为 (x / ec50)^hill 里，负数做非整数次幂会得到 NaN。
# 所以调整为(10^x / ec50)^hill，这只是一个变通的方法

rawdata |> 
  mutate(conc_um = conc_um + 0.01) |> 
  ggplot(aes(x = conc_um, y = viability, color = compound)) +
  geom_point(data = d, aes(x = conc_um + 0.01, y = mean)) +
  geom_smooth(
    method = "nls",
    formula = y ~ bottom + (top - bottom) / (1 + (10^x / ec50)^hill),
    method.args = list(
      start     = c(bottom = -0.2, top = 1, ec50 = 500, hill = 1),
      algorithm = "port",
      control   = list(maxiter = 1000)
    ),
    se = FALSE,
    linewidth = 0.9
  ) +
  scale_x_log10()

