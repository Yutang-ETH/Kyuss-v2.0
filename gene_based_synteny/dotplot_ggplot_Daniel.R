
rm(list = ls())

setwd("P:/Yutangchen/Kyuss_data/kyuss_nextpolish/liftofftools/liftofftools_output")

library(ggplot2)
library(dplyr)

df <- read.delim("gene_order", skip = 2, header = F)

mutate(df, V4 = factor(V4, levels = (paste0('chr', 1:7))),
       V5 = factor(V5, levels = rev(paste0('chr', 1:7)))) %>%
  ggplot(aes(x = V2, y = V3, color = V6)) +
  geom_point(size = .25) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  scale_color_viridis_c(direction = -1) +
  labs(x = "Kyuss v1.0", y = "Kyuss v2.0", color = "Sequence\nidentity") +
  facet_grid(V5 ~ V4, scales = 'free')

ggsave(filename = 'dotplot.svg', width = 170, height = 130, units = 'mm', dpi = 300, bg = 'transparent')

