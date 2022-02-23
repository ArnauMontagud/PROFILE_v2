---
title: "Analysis of experimental validation, Figure 7"
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

# Endpoint cell viability measurements ----
end1pre <-
  read.table(
    "./LNCAPendpoint.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
end_ids <-
  read.table(
    "./LNCAPendpoint_id.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  ) %>% pivot_longer(-X) %>% mutate(., id = as.integer(gsub("X", "", .$name))) %>% rename(., drug =
                                                                                            X, nM = value) %>% select(-name)

end1 <- full_join(end1pre, end_ids)
end1$nM <- as.integer(end1$nM)
end1$count_noise <-
  end1$count - (end1 %>% filter(drug == "NoCell") %>% .$count %>% mean())
end1$count_norm <-
  end1$count_noise / (end1 %>% filter(drug == "kontroll") %>% .$count_noise %>% mean())
end2 <- end1 %>% filter(!(drug=="NoCell")) %>% mutate(id=factor(id, levels=c(0,5,4,3,2,1))) %>% filter(!(is.na(count)))
end2[is.na(end2$nM),]$id <- 0
end2[is.na(end2$nM),]$nM <- 0.0

cols <- c(
  "17-DMAG" = pal_uchicago("default")(6)[1],
  "NMS-E973" = pal_uchicago("default")(6)[4],
  "PI-103" = pal_uchicago("default")(6)[5],
  "Pictilisib" = pal_uchicago("default")(6)[6]
)

`17-DMAGa` <-
  ggplot() +
  theme_bw() +
  scale_y_continuous(limits = c(0.25, 1.26)) +
  scale_x_discrete(labels=c("0" = "0", "5" = "1", "4" = "5", "3" = "25", "2" = "125", "1" = "625")) +
  geom_half_point(data=end2 %>% filter(drug=="17-DMAG" | drug=="kontroll"), aes(x = id, y = count_norm), colour = cols[1], transformation = position_jitter(width = 0.1,height = 0,seed = 4), size = 3) +
  geom_half_boxplot(data=end2 %>% filter(drug=="17-DMAG" | drug=="kontroll"), aes(x = id, y = count_norm), fill = cols[1], colour = cols[1], alpha = .3) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(title = "17-DMAG", x = "Dose (nM)", y = "Normalised cell viability") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",text=element_text(size=16),axis.text =element_text(colour = "black"),panel.background = element_blank(),panel.grid.minor = element_blank())

`NMS-E973a` <-
  ggplot() +
  theme_bw() +
  scale_y_continuous(limits = c(0.25, 1.26)) +
  scale_x_discrete(labels=c("0" = "0", "5" = "2", "4" = "8", "3" = "32", "2" = "128", "1" = "512")) +
  geom_half_point(data=end2 %>% filter(drug=="NMS-E973" | drug=="kontroll"), aes(x = id, y = count_norm), colour = cols[2], transformation = position_jitter(width = 0.1,height = 0,seed = 4), size = 3) +
  geom_half_boxplot(data=end2 %>% filter(drug=="NMS-E973" | drug=="kontroll"), aes(x = id, y = count_norm), fill = cols[2], colour = cols[2], alpha = .3) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(title = "NMS-E973", x = "Dose (nM)", y = "Normalised cell viability") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none",text=element_text(size=16),axis.text =element_text(colour = "black"),panel.background = element_blank(),panel.grid.minor = element_blank())

`PI-103a` <-
  ggplot() +
  theme_bw() +
  scale_y_continuous(limits = c(0.25, 1.25)) +
  scale_x_discrete(labels=c("0" = "0", "5" = "1", "4" = "5", "3" = "25", "2" = "125", "1" = "625")) +
  geom_half_point(data=end2 %>% filter(drug=="PI-103" | drug=="kontroll"), aes(x = id, y = count_norm), colour = cols[3], transformation = position_jitter(width = 0.12,height = 0,seed = 4), size = 3) +
  geom_half_boxplot(data=end2 %>% filter(drug=="PI-103" | drug=="kontroll"), aes(x = id, y = count_norm), fill = cols[3], colour = cols[3], alpha = .3) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(title = "PI-103", x = "Dose (nM)", y = "Normalised cell viability") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none",text=element_text(size=16),axis.text =element_text(colour = "black"),panel.background = element_blank(),panel.grid.minor = element_blank())

Pictilisiba <-
  ggplot() +
  theme_bw() +
  scale_y_continuous(limits = c(0.25, 1.25)) +
  scale_x_discrete(labels=c("0" = "0", "5" = "1", "4" = "5", "3" = "25", "2" = "125", "1" = "625")) +
  geom_half_point(data=end2 %>% filter(drug=="Pictilisib" | drug=="kontroll"), aes(x = id, y = count_norm), colour = cols[4], transformation = position_jitter(width = 0.12,height = 0,seed = 4), size = 3) +
  geom_half_boxplot(data=end2 %>% filter(drug=="Pictilisib" | drug=="kontroll"), aes(x = id, y = count_norm), fill = cols[4], colour = cols[4], alpha = .3) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(title = "Pictilisib", x = "Dose (nM)", y = "Normalised cell viability") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none",text=element_text(size=16),axis.text =element_text(colour = "black"),panel.background = element_blank(),panel.grid.minor = element_blank())

patchwork = (`17-DMAGa` + `PI-103a`) / (`NMS-E973a` + Pictilisiba)

# Remove titles from subplots
patchwork[[1]][[1]] = patchwork[[1]][[1]] + theme(axis.title.x = element_blank())
patchwork[[1]][[2]] = patchwork[[1]][[2]] + theme(
  axis.title.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)
patchwork[[2]][[2]] = patchwork[[2]][[2]] + theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

# Save Figure 7:
png(
  filename = "Fig7_endpoints.png",
  units = "in",
  width = 10,
  height = 10,
  res = 300
)
patchwork + plot_annotation(tag_levels = "A")
dev.off()
