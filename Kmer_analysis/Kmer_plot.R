library(ggplot2)
library(tidyr)
library(ggsci)

kmer <- read.csv("kmer_plot.txt",header=T,sep="\t")

kmer_long <- kmer %>%
  pivot_longer(cols = -x, names_to = "对象", values_to = "y")

ggplot(kmer_long, aes(x = x, y = y, color = 对象)) +
  geom_line(size = 1) +
#  geom_smooth(method = "loess", se = FALSE, span = 0.5) +
  scale_color_npg() +   # 核心：应用Nature (NPG) 配色方案[reference:2]
  labs(title = "Nature风格平滑折线图", x = "X 坐标", y = "Y 坐标") +
  theme_minimal()
