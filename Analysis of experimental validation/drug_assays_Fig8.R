---
title: "Analysis of experimental validation, Figure 8"
author: "Arnau Montagud, arnau.montagud@gmail.com"
date: "February 2022"
---
  
  # Note: All the scripts and data can be downloaded from https://github.com/ArnauMontagud/PROFILE_v2
  
  rm(list=ls(all=TRUE))

# load libraries ----
if (!require("pacman"))
  install.packages("pacman")
list.of.packages <-
  c(
    "tidyverse",
    "magrittr",
    "patchwork",
    "wesanderson",
    "RColorBrewer",
    "devtools",
    "scales"
  )
pacman::p_load(list.of.packages, character.only = TRUE)
devtools::install_github('erocoar/gghalves')
devtools::install_github('road2stat/ggsci')
library(gghalves)
library(ggsci)


# Real-time cell electronic sensing (RT-CES) cytotoxicity assay ----
b1pre <-
  read.table(
    "./LNCAP_RT.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
b1_ids <-
  b1pre %>% .[1:2, ] %>% .[-c(1:2)] %>% t() %>% as.data.frame() %>% rename(., drug =
                                                                             1, uM = 2) %>% mutate(cell = row.names(.))
b1pre$min <- NA
b1pre$min[-(1:2)] <-
  c(as.matrix(read.table(
    text = b1pre$Time.Interval[-(1:2)], sep = ":"
  )) %*% c(1, 1 / 60, 1 / 3660))
b1 <-
  b1pre %>% select(., min, everything(), -contains("Time")) %>% t() %>% as.data.frame()
colnames(b1) <- b1[1, ]
b1 %<>% .[-1, ]
colnames(b1)[1] <- "drug"
colnames(b1)[2] <- "uM"
b1 %<>% mutate(cell = row.names(.)) %>% select(cell, everything())
b1long <-
  b1 %>% pivot_longer(., !(c(drug, uM, cell))) %>% rename(., time = name, CI =
                                                            value) %>% as.data.frame()
b1long[which(b1long$drug == "Control"), ]$uM <- 0
b1long$uM <- factor(b1long$uM, levels = c("0", "3.3", "10", "30"))
b1long$time <- as.numeric(b1long$time)
b1long$CI <- as.numeric(b1long$CI)

times <- c("25:41:08", "49:18:09", "74:16:37", "97:17:03")
times_min <-
  c(as.matrix(read.table(text = times, sep = ":")) %*% c(1, 1 / 60, 1 / 3660))
times_hour <- c(times_min[1], 48.05246, 72.02869, 96.03388)
times_hour <-
  c(times_min[1], times_min[1] + 24, times_min[1] + 48, times_min[1] + 72)
times_hour2 <- times_hour - times_min[1]

b1long %<>% mutate(time2 = time - times_min[1])
b1long_mean <- b1long %>% filter(grepl("Y.", cell))
b1long_sd <- b1long %>% filter(grepl("SD.", cell))

b2 <- left_join(
  b1long %>% filter(grepl("Y.", cell)) %>% rename(mean = CI) %>% .[, -1],
  b1long %>% filter(grepl("SD.", cell)) %>% rename(sd = CI) %>% .[, -1]
)

summary_24 <- left_join(
  b1long %>% filter(time > 49 &
                      time < 49.1) %>% filter(grepl("Y.", cell)) %>% rename(mean = CI) %>% .[, -1],
  b1long %>% filter(time > 49 &
                      time < 49.1) %>% filter(grepl("SD.", cell)) %>% rename(sd = CI) %>% .[, -1]
)


summary_48 <- left_join(
  b1long %>% filter(time > 73 &
                      time < 73.1) %>% filter(grepl("Y.", cell)) %>% rename(mean = CI) %>% .[, -1],
  b1long %>% filter(time > 73 &
                      time < 73.1) %>% filter(grepl("SD.", cell)) %>% rename(sd = CI) %>% .[, -1]
)

summary_72 <- left_join(
  b1long %>% filter(time > 97 &
                      time < 97.1) %>% filter(grepl("Y.", cell)) %>% rename(mean = CI) %>% .[, -1],
  b1long %>% filter(time > 97 &
                      time < 97.1) %>% filter(grepl("SD.", cell)) %>% rename(sd = CI) %>% .[, -1]
)

pal <- pal_uchicago("default")(9)

# 17-DMAG:
cols <- c(
  "0" = pal[2],
  "3.3" = "#fa0000",
  "10" = "#8f0000",
  "30" = "#350000"
)

subdata1 <- b2 %>% filter(drug == "17DMAG" | drug == "Control")
`17-DMAGt_ribbon` <- ggplot() +
  theme_bw() +
  geom_ribbon(
    data = subdata1,
    aes(
      x = time2,
      y = subdata1$mean,
      ymin = subdata1$mean - subdata1$sd,
      ymax = subdata1$mean + subdata1$sd,
      fill = uM
    ),
    alpha = 0.2
  ) +
  geom_point(data = subdata1, aes(time2, mean, colour = uM)) +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  scale_fill_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`17-DMAGt` <-
  ggplot(
    data = b1long %>% filter(!grepl("SD.", cell)) %>% filter(drug == "17DMAG" |
                                                               drug == "Control"),
    aes(time2, CI, colour = uM)
  ) +
  theme_bw() +
  geom_point() +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "17-DMAG", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`17-DMAG24` <-
  ggplot(summary_24 %>% filter(drug == "17DMAG" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "24 hours", x = "Drug dose (uM)", y = "Cell Index (a.u.)") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`17-DMAG48` <-
  ggplot(summary_48 %>% filter(drug == "17DMAG" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "48 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`17-DMAG72` <-
  ggplot(summary_72 %>% filter(drug == "17DMAG" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "72 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# NMS-E973:
cols <- c(
  "0" = pal[2],
  "3.3" = "#c1c960",
  "10" = "#6e7337",
  "30" = "#292b14"
)

subdata2 <- b2 %>% filter(drug == "NMS-E973" | drug == "Control")

`NMS-E973t_ribbon` <- ggplot() +
  theme_bw() +
  geom_ribbon(
    data = subdata2,
    aes(
      x = time2,
      y = subdata2$mean,
      ymin = subdata2$mean - subdata2$sd,
      ymax = subdata2$mean + subdata2$sd,
      fill = uM
    ),
    alpha = 0.2
  ) +
  geom_point(data = subdata2, aes(time2, mean, colour = uM)) +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  scale_fill_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`NMS-E973t` <-
  ggplot(
    data = b1long %>% filter(!grepl("SD.", cell)) %>% filter(drug == "NMS-E973" |
                                                               drug == "Control"),
    aes(time2, CI, colour = uM)
  ) +
  theme_bw() +
  geom_point() +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "NMS-E973", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`NMS-E97324` <-
  ggplot(summary_24 %>% filter(drug == "NMS-E973" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "24 hours", x = "Drug dose (uM)", y = "Cell Index (a.u.)") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`NMS-E97348` <-
  ggplot(summary_48 %>% filter(drug == "NMS-E973" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "48 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`NMS-E97372` <-
  ggplot(summary_72 %>% filter(drug == "NMS-E973" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "72 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

patchwork_HSP2 = `17-DMAGt_ribbon` / (`17-DMAG24` + `17-DMAG48` + `17-DMAG72`) / `NMS-E973t_ribbon` / (`NMS-E97324` + `NMS-E97348` + `NMS-E97372`)
patchwork_HSP2[[1]] = patchwork_HSP2[[1]] + theme(plot.title = element_blank())
patchwork_HSP2[[2]][[1]] = patchwork_HSP2[[2]][[1]] + theme(plot.title = element_blank())
patchwork_HSP2[[2]][[2]] = patchwork_HSP2[[2]][[2]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_HSP2[[2]][[3]] = patchwork_HSP2[[2]][[3]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_HSP2[[3]] = patchwork_HSP2[[3]] + theme(plot.title = element_blank())
patchwork_HSP2[[4]][[1]] = patchwork_HSP2[[4]][[1]] + theme(plot.title = element_blank())
patchwork_HSP2[[4]][[2]] = patchwork_HSP2[[4]][[2]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_HSP2[[4]][[3]] = patchwork_HSP2[[4]][[3]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
# Save Figure S40:
png(
  filename = "FigS40_HSP_RT.png",
  units = "in",
  width = 15,
  height = 20,
  res = 300
)
patchwork_HSP2 + plot_annotation(tag_levels = "A") + plot_layout(heights = c(.7, .3, .7, .3))
dev.off()

patchwork_HSP3 = `17-DMAGt_ribbon` / `NMS-E973t_ribbon`
patchwork_HSP3[[1]] = patchwork_HSP3[[1]] + theme(plot.title = element_blank())
patchwork_HSP3[[2]] = patchwork_HSP3[[2]] + theme(plot.title = element_blank())

# Save Figure 8:
png(
  filename = "Fig8_HSP_RT.png",
  units = "in",
  width = 15,
  height = 15,
  res = 300
)
patchwork_HSP3 + plot_annotation(tag_levels = "A")
dev.off()

# PI-103:
cols <- c(
  "0" = pal[2],
  "3.3" = "#28baff",
  "10" = "#176a92",
  "30" = "#082736"
)

subdata3 <- b2 %>% filter(drug == "PI-103" | drug == "Control")

`PI-103t_ribbon` <- ggplot() +
  theme_bw() +
  geom_ribbon(
    data = subdata3,
    aes(
      x = time2,
      y = subdata3$mean,
      ymin = subdata3$mean - subdata3$sd,
      ymax = subdata3$mean + subdata3$sd,
      fill = uM
    ),
    alpha = 0.2
  ) +
  geom_point(data = subdata3, aes(time2, mean, colour = uM)) +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  scale_fill_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`PI-103t` <-
  ggplot(
    data = b1long %>% filter(!grepl("SD.", cell)) %>% filter(drug == "PI-103" |
                                                               drug == "Control"),
    aes(time2, CI, colour = uM)
  ) +
  theme_bw() +
  geom_point() +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "PI-103", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`PI-10324` <-
  ggplot(summary_24 %>% filter(drug == "PI-103" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "24 hours", x = "Drug dose (uM)", y = "Cell Index (a.u.)") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`PI-10348` <-
  ggplot(summary_48 %>% filter(drug == "PI-103" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "48 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`PI-10372` <-
  ggplot(summary_72 %>% filter(drug == "PI-103" | drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "72 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

# Pictilisib:

cols <- c(
  "0" = pal[2],
  "3.3" = "#ffc641",
  "10" = "#cc7125",
  "30" = "#4c2a0e"
)

subdata4 <- b2 %>% filter(drug == "Pictilisib" | drug == "Control")
`Pictilisibt_ribbon` <- ggplot() +
  theme_bw() +
  geom_ribbon(
    data = subdata4,
    aes(
      x = time2,
      y = subdata4$mean,
      ymin = subdata4$mean - subdata4$sd,
      ymax = subdata4$mean + subdata4$sd,
      fill = uM
    ),
    alpha = 0.2
  ) +
  geom_point(data = subdata4, aes(time2, mean, colour = uM)) +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  scale_fill_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`Pictilisibt` <-
  ggplot(
    data = b1long %>% filter(!grepl("SD.", cell)) %>% filter(
      drug == "Pictilisib" |
        drug == "Control"),
    aes(time2, CI, colour = uM)
  ) +
  theme_bw() +
  geom_point() +
  scale_x_continuous(breaks = seq(-24, 108, 12)) +
  scale_colour_manual(name = "Drug dose (uM)", values = cols) +
  geom_vline(
    xintercept = times_hour2,
    colour = c(pal[3], pal[7], pal[7], pal[7]),
    size = 1,
    linetype = "dashed"
  ) +
  ylim(NA, 4) +
  labs(title = "Pictilisib", x = "Time (hours)", y = "Cell Index (a.u.)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.15, 0.8),
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`Pictilisib24` <-
  ggplot(summary_24 %>% filter(drug == "Pictilisib" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "24 hours", x = "Drug dose (uM)", y = "Cell Index (a.u.)") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`Pictilisib48` <-
  ggplot(summary_48 %>% filter(drug == "Pictilisib" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "48 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

`Pictilisib72` <-
  ggplot(summary_72 %>% filter(drug == "Pictilisib" |
                                 drug == "Control")) +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_point(aes(uM, mean, colour = uM), size = 8) +
  geom_errorbar(
    aes(
      x = uM,
      y = mean,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = uM
    ),
    width = .4,
    size = 2
  ) +
  scale_y_continuous(limits = c(-1, 4)) +
  labs(title = "72 hours", x = "Drug dose (uM)", y = "Cell Index") + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 28),
    axis.text = element_text(colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank()
  )

patchwork_PI3K = `PI-103t_ribbon` / (`PI-10324` + `PI-10348` + `PI-10372`) / `Pictilisibt_ribbon` / (`Pictilisib24` + `Pictilisib48` + `Pictilisib72`)
patchwork_PI3K[[1]] = patchwork_PI3K[[1]] + theme(plot.title = element_blank())
patchwork_PI3K[[2]][[1]] = patchwork_PI3K[[2]][[1]] + theme(plot.title = element_blank())
patchwork_PI3K[[2]][[2]] = patchwork_PI3K[[2]][[2]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_PI3K[[2]][[3]] = patchwork_PI3K[[2]][[3]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_PI3K[[3]] = patchwork_PI3K[[3]] + theme(plot.title = element_blank())
patchwork_PI3K[[4]][[1]] = patchwork_PI3K[[4]][[1]] + theme(plot.title = element_blank())
patchwork_PI3K[[4]][[2]] = patchwork_PI3K[[4]][[2]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork_PI3K[[4]][[3]] = patchwork_PI3K[[4]][[3]] + theme(
  plot.title = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank()
)

# Save Figure S41:
png(
  filename = "FigS41_PI3K_RT.png",
  units = "in",
  width = 15,
  height = 20,
  res = 300
)
patchwork_PI3K + plot_annotation(tag_levels = "A") + plot_layout(heights = c(.7, .3, .7, .3))
dev.off()

patchwork_PI3K2 = `PI-103t_ribbon` / `Pictilisibt_ribbon`
patchwork_PI3K2[[1]] = patchwork_PI3K2[[1]] + theme(plot.title = element_blank())
patchwork_PI3K2[[2]] = patchwork_PI3K2[[2]] + theme(plot.title = element_blank())


# Save Figure 9:
png(
  filename = "Fig9_PI3K_RT.png",
  units = "in",
  width = 15,
  height = 15,
  res = 300
)
patchwork_PI3K2 + plot_annotation(tag_levels = "A")
dev.off()
