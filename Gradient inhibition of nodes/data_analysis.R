---
title: "Analysis of gradient inhibition simulation"
author: "Arnau Montagud, arnau.montagud@gmail.com"
date: "June 2021"
---

rm(list = ls(all = TRUE))

# import and load packages ----
if (!require("pacman"))
  install.packages("pacman")
list.of.packages <-
  c("tidyverse",
    "magrittr",
    "reshape2",
    "dplyr",
    "RColorBrewer",
    "patchwork",
    "ggpubr")
pacman::p_load(list.of.packages, character.only = TRUE)
options(scipen = 999)

# load data ----
load("drugs_single.RData")
load("drugs_double.RData")

# Plotting WT ----
WT_pheno %<>% mutate(activ = replace(activ, is.na(activ) == T, "rand")) %<>%
  mutate(RNA = replace(RNA, is.na(RNA) == T, "PC20191203")) %<>%
  mutate(RNA = replace(RNA, RNA == "mutRNA", "LNCAP")) %>%
  rename(WT = RNA)

dropcol <- c("Time", "drug", "dose", "RNA", "TH")
WT_pheno_large <-
  WT_pheno  %>% select(-one_of(dropcol)) %>% gather(., key = Phenotype, value, -c(WT, activ))
WT_pheno_large %<>% mutate(activ = replace(activ, activ == "AR", "Androgen")) %>% mutate(activ =
                                                                                           replace(activ, activ == "AR_EGF", "Androgen_EGF"))

WT_pheno_large$activ <-
  factor(WT_pheno_large$activ,
         levels = c("rand", "00", "EGF", "Androgen_EGF", "Androgen"))

phenos <-
  c("Proliferation",
    "Invasion",
    "Metastasis",
    "Migration",
    "DNA_Repair",
    "Apoptosis")
WT_pheno_large %<>% filter(Phenotype %in% phenos)

#Save Figure S31
png(
  filename = "FigS31_WT_pheno.png",
  units = "px",
  width = 1000,
  height = 900,
  res = 100
)
ggplot(WT_pheno_large, aes(factor(0), activ)) +
  geom_tile(aes(fill = value), colour = "black") +
  facet_grid(WT ~ Phenotype) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white",
                      high = "firebrick",
                      limits = c(0, 1)) +
  labs(title = "WT phenotypes' probabilities",
       x = "Phenotype",
       y = "Growth additive",
       fill = "Cell line\nscore") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.text.x = element_blank()
  )
dev.off()

# Plotting single drug ----

phenos <-
  c("Proliferation",
    "Invasion",
    "Metastasis",
    "Migration",
    "DNA_Repair",
    "Apoptosis")
single_large %<>% filter(Phenotype %in% phenos)
single_large %>% .$Phenotype %>% unique()

single_large1 <-
  single_large %>% mutate (activ = fct_recode(activ, Androgen = "AR", Androgen_EGF = "AR_EGF"))

#Save Figure S30
png(
  filename = "FigS30_Single_all.png",
  units = "px",
  width = 900,
  height = 900,
  res = 100
)
ggplot(single_large1, aes(dose, activ)) +
  geom_tile(aes(fill = value), colour = "black") +
  facet_grid(drug ~ Phenotype) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     labels = c(0, 20, 40, 60, 80, 100)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "blue",
    high = "firebrick",
    mid = "white",
    limits = c(-1, 1)
  ) +
  labs(title = "Phenotype variations upon node inhibition",
       x = "Node inhibition (%)",
       y = "Growth additive",
       fill = "Treated -\nuntreated\ncell line\nscore") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.text.x = element_text(angle = 90),
    strip.text.y.right = element_text(angle = 0)
  )
dev.off()

# Single only EGF
single_EGF <- single_large %>% filter(activ == "EGF")
single_EGF$drug_ord <-
  with(single_EGF, factor(drug, levels = rev(sort(unique(
    drug
  )))))

#Save Figure S29
png(
  filename = "FigS29_Single_EGF.png",
  units = "px",
  width = 900,
  height = 700,
  res = 100
)
ggplot(single_EGF, aes(dose, drug_ord)) +
  geom_tile(aes(fill = value), colour = "black") +
  facet_grid(~ Phenotype) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     labels = c(0, 20, 40, 60, 80, 100)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "blue",
    high = "firebrick",
    mid = "white",
    limits = c(-1, 1)
  ) +
  labs(title = "Phenotype variations upon node inhibition",
       x = "Node inhibition (%)",
       y = "Growth additive",
       fill = "Treated -\nuntreated\ncell line\nscore") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.text.x = element_text(angle = 90)
  )
dev.off()

drugnames <- single_large %>% .$drug %>% unique()

pdf("Single drugs.pdf", onefile = T)
for (j in 1:length(drugnames)) {
  print(paste0(drugnames[j]))
  print(
    ggplot(single_large %>% filter(drug == drugnames[j]), aes(dose, activ)) +
      geom_tile(aes(fill = value), colour = "black") +
      facet_grid(~ Phenotype) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_gradient2(low = "blue", high = "firebrick", mid = "white") +
      labs(
        title = paste0("Phenotype variations upon ", drugnames[j], " inhibition"),
        x = paste0(drugnames[j], " inhibition"),
        y = "AR or EGF presence",
        fill = "mutant -\ncell line"
      ) +
      theme(
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5)
      )
  )
}
dev.off()


# Plotting double drugs ----

tot<-single_large %>% filter(activ =="EGF") %>% mutate(drug1=drug) %>% 
  rename(drug2=drug) %>% mutate(dose1=dose/10) %>% rename(dose2=dose) %>% 
  mutate(dose2=dose2/10) %>% 
  select(drug1,drug2,dose1,dose2,RNA,activ,Phenotype,value,value2)
tot <- rbind(tot, double_large %>% filter(activ =="EGF"))

double_drugs <- double_large
drugnames <-
  c(double_drugs %>% .$drug1, double_drugs %>% .$drug2) %>% unique()
pdf("Double drugs.pdf", onefile = T)
for (j in 1:length(drugnames)) {
  for (k in 1:length(drugnames)) {
    if (k > j) {
      print(paste0(drugnames[j], " ", drugnames[k]))
      print(
        ggplot(
          double_drugs %>% filter(drug1 == drugnames[j], drug2 == drugnames[k]),
          aes(dose1, dose2)
        ) +
          geom_tile(aes(fill = value), colour = "black") +
          facet_grid(activ ~ Phenotype) +
          scale_fill_gradient2(
            low = "blue",
            high = "firebrick",
            mid = "white",
            limits = c(-1, 1)
          ) +
          labs(
            title = paste0(
              "Phenotype variations upon ",
              drugnames[j],
              " and ",
              drugnames[k],
              " inhibition"
            ),
            x = paste0(drugnames[j], " inhibition"),
            y = paste0(drugnames[k], " inhibition"),
            fill = "mutant -\ncell line\nscore"
          ) +
          theme(
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(fill = NA),
            plot.title = element_text(hjust = 0.5)
          )
      )
    }
  }
}
dev.off()

# Synergies: Highest Single Agent i Bliss independence ----

drugnames <-
  c(double_large %>% .$drug1, double_large %>% .$drug2) %>% unique()
phenos <-
  c("Proliferation",
    "Invasion",
    "Metastasis",
    "Migration",
    "DNA_Repair",
    "Apoptosis")
double_Bliss <- double_large %>% filter(activ == "EGF")
single_Bliss <- single_large %>% filter(activ == "EGF")

desnorm <-
  WT_pheno %>% filter(WT == "LNCAP" &
                        activ == "EGF") %>% .$Proliferation

tot_Bliss <- data.frame()
for (j in 1:length(drugnames)) {
  for (k in 1:length(drugnames)) {
    if (k > j) {
      print(paste0(drugnames[j], " ", drugnames[k]))
      for (l in 1:length(phenos)) {
        z1 <- double_Bliss %>% filter(Phenotype == phenos[l]) %>%
          filter(drug1 == drugnames[j] & drug2 == drugnames[k]) %>%
          select(-c(RNA, activ, value2)) %>% rename(., both = value)
        z2 <- single_Bliss %>% filter(Phenotype == phenos[l]) %>%
          filter(drug == drugnames[j] | drug == drugnames[k]) %>%
          select(-c(RNA, activ, value2)) %>% spread(., key = drug, value) %>%
          rename(dose1 = dose,
                 drug1v = drugnames[j],
                 drug2v = drugnames[k]) %>%
          mutate(dose1 = dose1 / 10, dose2 = dose1)
        z3 <-
          merge(z1, z2 %>% select(.,-c(dose2, drug2v)) %>% .[complete.cases(.), ], by = c("dose1"))
        z3 <-
          merge(z3, z2 %>% select(.,-c(dose1, drug1v)) %>% .[complete.cases(.), ], by = c("dose2"))
        z3 %<>% select(-c(Phenotype, Phenotype.y)) %>% rename(Phenotype =
                                                                Phenotype.x)
        z3$HSA <-
          (apply(abs(z3 %>% select(
            c(drug1v, drug2v)
          )), 1, max)) / abs(z3$both)
        z3$drug1v2 <- z3$drug1v + desnorm
        z3$drug2v2 <- z3$drug2v + desnorm
        z3$both2 <- z3$both + desnorm
        z3$Bliss <-
          (z3$drug1v2 + z3$drug2v2 - (z3$drug1v2 * z3$drug2v2)) / z3$both2
        tot_Bliss <- rbind(tot_Bliss, z3)
        rm(z1, z2, z3)
      }
    }
  }
}

tot_Bliss %<>% mutate(Bliss1a=replace(Bliss, Bliss>2, 2)) %>% mutate(Bliss1a=replace(Bliss1a, dose1==0.0, NA)) %>% 
                            mutate(Bliss1a=replace(Bliss1a, dose2==0.0, NA)) %>% 
                            mutate(Bliss = Bliss1a, Bliss1a = NULL)
tot_Bliss_prolif <- tot_Bliss %>% filter(Phenotype=="Proliferation")
tot_Bliss_apop <- tot_Bliss %>% filter(Phenotype=="Apoptosis")

# Plotting figures ----

load("drugs_figures_datasets.RData")

#Save Figure S31
png(
  filename = "Fig_S31_Double_EGF_Proliferation.png",
  units = "in",
  width = 11,
  height = 10,
  res = 100
)
ggplot(tot %>% filter(Phenotype == "Proliferation"), aes(dose1, dose2)) +
  facet_grid(drug2 ~ drug1, switch = "both") +
  geom_tile(aes(fill = value), colour = "black") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradient2(
    low = "blue",
    high = "firebrick",
    mid = "white",
    limits = c(-1, 1)
  ) +
  labs(
    title = paste0("Proliferation phenotype variations upon drugs inhibition"),
    x = "Node inhibition (%)",
    y = "Node inhibition (%)",
    fill = "Treated -\nuntreated\ncell line\nscore"
  ) +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    axis.text.x = element_text(angle = 90)
  )
dev.off()

#Save Figure S32
png(
  filename = "FigS32_Double_EGF_Apoptosis.png",
  units = "in",
  width = 11,
  height = 10,
  res = 100
)
ggplot(tot %>% filter(Phenotype == "Apoptosis"), aes(dose1, dose2)) +
  facet_grid(drug2 ~ drug1, switch = "both") +
  geom_tile(aes(fill = value), colour = "black") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradient2(low = "blue", high = "firebrick", mid = "white") +
  labs(
    title = paste0("Apoptosis phenotype variations upon drugs inhibition"),
    x = "Node inhibition (%)",
    y = "Node inhibition (%)",
    fill = "Treated -\nuntreated\ncell line\nscore"
  ) +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    axis.text.x = element_text(angle = 90)
  )
dev.off()

bliss_pro <- ggplot(tot_Bliss_prolif, aes(dose1, dose2)) +
  geom_tile(aes(fill = Bliss1), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradientn(
    colours = c("#8c5ba9", "white", "#7fbf7b"),
    breaks = (c(0.0, 0.5, 1, 1.5, 2)),
    labels = (c(0.0, 0.5, 1, 1.5, ">2"))
  ) +
  labs(x = "Node inhibition (%)",
       y = "Node inhibition (%)",
       fill = "Bliss score") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    axis.text.x = element_text(angle = 90)
  )

# Save Figure S33
png(
  filename = "FigS33_Bliss_proliferation.png",
  units = "px",
  width = 1150,
  height = 1000,
  res = 100
)
bliss_pro
dev.off()

bliss_apop <- ggplot(tot_Bliss_apop, aes(dose1, dose2)) +
  geom_tile(aes(fill = Bliss1), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradientn(
    colours = c("#8c5ba9", "white", "#7fbf7b"),
    breaks = (c(0, 0.5, 1, 1.5, 2)),
    labels = (c(0, 0.5, 1, 1.5, ">2"))
  ) +
  labs(x = "Node inhibition (%)",
       y = "Node inhibition (%)",
       fill = "Bliss score") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    axis.text.x = element_text(angle = 90)
  )

# Save Figure S34
png(
  filename = "FigS34_Bliss_apoptosis.png",
  units = "px",
  width = 1150,
  height = 1000,
  res = 100
)
bliss_apop
dev.off()

interesants_prolif <-
  tot %>% filter(drug1 == "ERK" & drug2 == "MYC_MAX")
interesants_bliss <-
  tot_Bliss_prolif %>% filter(drug1 == "ERK" & drug2 == "MYC_MAX")

prolif <-
  ggplot(interesants_prolif %>% filter(Phenotype == "Proliferation"),
         aes(dose1, dose2)) +
  geom_tile(aes(fill = value), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradient2(low = "blue",
                       high = "white",
                       limits = c(-1, 0.1)) +
  labs(x = "ERK node inhibition (%)",
       y = "MYC_MAX node\ninhibition (%)",
       fill = "Treated -\nuntreated\ncell line\nscore") +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5)
  )
prolif

prolif_bliss <- ggplot(interesants_bliss, aes(dose1, dose2)) +
  geom_tile(aes(fill = Bliss1), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradientn(
    colours = c("#8c5ba9", "white", "#7fbf7b"),
    breaks = (c(0, 0.5, 1, 1.5, 2)),
    labels = (c(0, 0.5, 1, 1.5, ">2")),
    limits = c(0, 2)
  ) +
  labs(x = "ERK node inhibition (%)",
       y = "MYC_MAX node\ninhibition (%)",
       fill = "Bliss\nscore") +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5)
  )
prolif_bliss

interesants_prolif <-
  tot %>% filter(drug1 == "HSPs" & drug2 == "PI3K")
interesants_bliss <-
  tot_Bliss_prolif %>% filter(drug1 == "HSPs" & drug2 == "PI3K")

prolif2 <-
  ggplot(interesants_prolif %>% filter(Phenotype == "Proliferation"),
         aes(dose1, dose2)) +
  geom_tile(aes(fill = value), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradient2(low = "blue",
                       high = "white",
                       limits = c(-1, 0.1)) +
  labs(x = "HSPs node inhibition (%)",
       y = "PI3K node inhibition (%)",
       fill = "Treated -\nuntreated\ncell line\nscore") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5)
  )
prolif2

prolif_Bliss1 <- ggplot(interesants_bliss, aes(dose1, dose2)) +
  geom_tile(aes(fill = Bliss1), colour = "black") +
  facet_grid(drug2 ~ drug1, switch = "both") +
  scale_x_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_y_continuous(
    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c(0, 20, 40, 60, 80, 100)
  ) +
  scale_fill_gradientn(
    colours = c("#8c5ba9", "white", "#7fbf7b"),
    breaks = (c(0, 0.5, 1, 1.5, 2)),
    labels = (c(0, 0.5, 1, 1.5, ">2")),
    limits = c(0, 2)
  ) +
  labs(x = "HSPs node inhibition (%)",
       y = "PI3K node inhibition (%)",
       fill = "Bliss\nscore") +
  theme(
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.text = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(hjust = 0.5)
  )
prolif_Bliss1

# Save Figure 5
png(
  filename = "Fig5_ERK,MYC_MAX,HSPs,PI3K.png",
  units = "in",
  width = 10,
  height = 10,
  res = 100
)
prolif + prolif2 + prolif_bliss  + prolif_Bliss1 + plot_annotation(tag_levels = 'A')
dev.off()
