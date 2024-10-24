#####figure1
rm(list=ls())

library(data.table)
library(readxl)
library(tidyverse)
library(scales)
library(patchwork)
library(colorspace)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)

mse_real <- read.csv("result/time.csv")
mse_real1 <- mse_real[,-1]

datmse1 <- data.frame(N = c(1:6),
                      DMAP_SA = mse_real1[, 1],
                      DMAP_SL = mse_real1[, 4],
                      gMAP= mse_real1[, 7])

datP_long <- reshape2::melt(datmse1, id.vars = "N", variable.name = "Method", value.name = "time")
datP_long$Method <- factor(datP_long$Method, 
                           levels = c("DMAP_SA", "DMAP_SL", "gMAP"), 
                           labels = c("DMAP-SA", "DMAP-SL", "gMAP"))

p1 <- ggplot(data = datP_long) +
  geom_point(aes(x = N, y = time, color = Method, shape = Method), size = 2) +
  geom_line(aes(x = N, y = time, color = Method, linetype = Method), linewidth = 1) +
  labs(title = "(CM1)", x = "Number of machines", y = "Time") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("DMAP-SA" = "#F0E442", "DMAP-SL" = "#C699E7", "gMAP" = "#28BBD7")) +
  scale_shape_manual(values = c("DMAP-SA" = 16, "DMAP-SL" = 17, "gMAP" = 18)) +
  scale_linetype_manual(values = c("DMAP-SA" = 1, "DMAP-SL" = 2, "gMAP" = 3)) +
  scale_x_continuous(breaks = 1:6, labels = c(2, 4, 8, 16,32,64)) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank(),  # 移除网格线
        panel.border = element_rect(colour = "black", size = 1.2))


datmse2 <- data.frame(N = c(1:6),
                      DMAP_SA = mse_real1[, 2],
                      DMAP_SL = mse_real1[, 5],
                      gMAP = mse_real1[, 8])

datP_long <- reshape2::melt(datmse2, id.vars = "N", variable.name = "Method", value.name = "time")
datP_long$Method <- factor(datP_long$Method, 
                           levels = c("DMAP_SA", "DMAP_SL", "gMAP"), 
                           labels = c("DMAP-SA", "DMAP-SL", "gMAP"))

p2 <- ggplot(data = datP_long) +
  geom_point(aes(x = N, y = time, color = Method, shape = Method), size = 2) +
  geom_line(aes(x = N, y = time, color = Method, linetype = Method), linewidth = 1) +
  labs(title = "(CM2)", x = "Number of machines", y = "Time") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("DMAP-SA" = "#F0E442", "DMAP-SL" = "#C699E7", "gMAP" = "#28BBD7")) +
  scale_shape_manual(values = c("DMAP-SA" = 16, "DMAP-SL" = 17, "gMAP" = 18)) +
  scale_linetype_manual(values = c("DMAP-SA" = 1, "DMAP-SL" = 2, "gMAP" = 3)) +
  scale_x_continuous(breaks = 1:6, labels = c(2, 4, 8, 16,32,64)) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank(),  # 移除网格线
        panel.border = element_rect(colour = "black", size = 1.2))

datmse3 <- data.frame(N = c(1:6),
                      DMAP_SA =mse_real1[, 3],
                      DMAP_SL = mse_real1[, 6],
                      gMAP = mse_real1[, 9])

datP_long <- reshape2::melt(datmse3, id.vars = "N", variable.name = "Method", value.name = "time")
datP_long$Method <- factor(datP_long$Method, 
                           levels = c("DMAP_SA", "DMAP_SL", "gMAP"), 
                           labels = c("DMAP-SA", "DMAP-SL", "gMAP"))

p3 <- ggplot(data = datP_long) +
  geom_point(aes(x = N, y = time, color = Method, shape = Method), size = 2) +
  geom_line(aes(x = N, y = time, color = Method,linetype = Method), linewidth = 1) +
  labs(title = "(CM3)", x = "Number of machines", y = "Time") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("DMAP-SA" = "#F0E442", "DMAP-SL" = "#C699E7", "gMAP" = "#28BBD7")) +
  scale_shape_manual(values = c("DMAP-SA" = 16, "DMAP-SL" = 17, "gMAP" = 18)) +
  scale_linetype_manual(values = c("DMAP-SA" = 1, "DMAP-SL" = 2, "gMAP" = 3)) +
  scale_x_continuous(breaks = 1:6, labels = c(2,4, 8, 16,32,64)) +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank(),  # 移除网格线
        panel.border = element_rect(colour = "black", size = 1.2))+  # 加粗外框
  guides(color = guide_legend(title = NULL),   # 移除 color 图例标题
         shape = guide_legend(title = NULL),   # 移除 shape 图例标题
         linetype = guide_legend(title = NULL))  # 移除 linetype 图例标题

legend <- get_legend(p3 + theme(legend.position = "bottom"))

gt <- gtable(widths = unit(c(4, 4, 4), "null"), heights = unit(c(4, 1), "null"))

gt <- gtable_add_grob(gt, grobs = list(ggplotGrob(p1)), t = 1, l = 1)
gt <- gtable_add_grob(gt, grobs = list(ggplotGrob(p2)), t = 1, l = 2)
gt <- gtable_add_grob(gt, grobs = list(ggplotGrob(p3)), t = 1, l = 3)

gt <- gtable_add_grob(gt, grobs = list(legend), t = 2, l = 2)

grid.draw(gt)
