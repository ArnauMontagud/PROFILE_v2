---
title: "Gathering gradient inhibition simulations for single drugs"
author: "Arnau Montagud, arnau.montagud@gmail.com"
date: "June 2021"
---

rm(list = ls(all = TRUE))

# import and load packages ----
if (!require("pacman"))
  install.packages("pacman")
list.of.packages <-
  c("tidyr",
    "magrittr",
    "reshape2",
    "dplyr",
    "RColorBrewer",
    "stringr")
lapply(list.of.packages, require, character.only = TRUE)
options(scipen = 999)

# clean probtraj ----
cleaning_prob <- function(opt) {
  ## Parse required parameters
  if (!("input" %in% names(opt))) {
    stop(
      "No input file given! \n---> Specify an input file (.probtraj file from MaBoSS simulation) using '-i input.probtraj' or '--input input.protraj'"
    )
  }
  if (!("conditions" %in% names(opt))) {
    stop(
      "No conditions given to analyse! Specify at least one node with '-n \"conditions1,conditions,...\"' or '--conditions \"conditions1,conditions,...\"'"
    )
  }
  ## format parameters
  conditions <- unlist(strsplit(opt$conditions, split = ","))
  if (!("states" %in% names(opt))) {
    states <-
      as.character(rep(1, length(unlist(
        strsplit(opt$conditions, ",")
      ))))
  } else {
    states <- unlist(strsplit(opt$states, ","))
  }
  # data <- read.csv(opt$input,sep="\t", header=T)
  data <- opt$input
  ## Function to grep values for a specific node
  ##############################################################################
  averageProb <- function(data, conditions, states) {
    #Verify input data
    if (length(conditions) != length(states)) {
      stop(
        "Length of states should be the same as the total number of conditions in conditions,
           considering combined conditions separately."
      )
    }
    ## get columns that define states probabilities
    data_probCol <-
      !(grepl('Err', colnames(data))) & grepl('Prob', colnames(data))
    data_probCol_colnames <- colnames(data)[data_probCol]
    ## get columns of states fullfiling the conditions
    conditions_proba_list <-
      lapply(
        1:length(conditions),
        FUN = function(conditionN) {
          node <- unlist(strsplit(conditions[conditionN], '&'))
          value <- as.numeric(unlist(strsplit(states[conditionN], '&')))
          ## for each node
          data_probCol_colnamesFullfillCondition_list <-
            lapply(
              1:length(node),
              FUN = function(nodeN) {
                ## for each colname
                sapply(
                  data_probCol_colnames,
                  FUN = function(state) {
                    if (value[nodeN] == 1) {
                      return(node[nodeN] %in% unlist(strsplit(
                        gsub("Prob.", "", gsub(".$", "", state)), split = "\\.\\."
                      )))
                    }
                    if (value[nodeN] == 0) {
                      return(!(node[nodeN] %in% unlist(
                        strsplit(gsub(
                          "Prob.", "", gsub(".$", "", state)
                        ), split = "\\.\\.")
                      )))
                    }
                  }
                )
              }
            )
          data_probCol_colnamesFullfillCondition <-
            data_probCol_colnames[Reduce('*', data_probCol_colnamesFullfillCondition_list) >
                                    0]
          if (length(data_probCol_colnamesFullfillCondition) == 0) {
            condition_proba <- rep(0, length(data[, "Time"]))
          }
          if (length(data_probCol_colnamesFullfillCondition) == 1) {
            condition_proba <- data[, data_probCol_colnamesFullfillCondition]
          }
          if (length(data_probCol_colnamesFullfillCondition) > 1) {
            condition_proba <-
              apply(data[, data_probCol_colnamesFullfillCondition], 1, sum)
          }
          return(condition_proba)
        }
      )
    ## change names
    names(conditions_proba_list) <-
      unlist(lapply(
        1:length(conditions),
        FUN = function(conditionN) {
          node = unlist(strsplit(conditions[conditionN], '&'))
          value = as.numeric(unlist(strsplit(states[conditionN], '&')))
          return(paste0(node, ".", value, collapse = "_"))
        }
      ))
    #conditions_proba_df <- data.frame(Time=data[,"Time"],do.call(cbind, conditions_proba_list))
    conditions_proba_df <-
      data.frame(Time = data[, "Time"],
                 do.call(cbind, conditions_proba_list),
                 TH = data[, "TH"])
    return(conditions_proba_df)
  }
  ## extract trajectories
  conditions_proba_df <- averageProb(data, conditions, states)
  return(conditions_proba_df)
}

# WT ----
WT <-
  list.dirs(path = "./WT", recursive = F) %>% grep("_00|_AR|_EGF", ., value =
                                                     TRUE) %>% sub("./WT/", "", .)

cond <-
  "Proliferation,Angiogenesis,Invasion,Metastasis,Migration,Glycolysis,Hypermethylation,DNA_Repair,Quiescence,Apoptosis"
WT_pheno <- data.frame()
for (j in 1:length(WT)) {
  folder <- WT[j]
  a2 <- data.frame()
  file <-
    list.files(paste0(folder, "/")) %>% grep("_probtraj_table.csv", ., value =
                                               TRUE)
  folder_file <- paste0(folder, "/", file)
  res <-
    read.table(
      folder_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    ) %>% .[nrow(.), , drop = FALSE]
  optiona <- list(
    input = res,
    outprefix = file %>% sub("_probtraj_table.csv", "", .),
    conditions = cond
  )
  a1 <- cleaning_prob(optiona) %>% .[-ncol(.)]
  rownames(a1)[1] <- file %>% sub("_probtraj_table.csv", "", .)
  a1$drug <- NA
  a1$dose <- NA
  a1$RNA <- ifelse(grepl("mutRNA", file) == T,
                   "mutRNA",
                   ifelse(grepl("RNA", file) == T, "RNA", NA))
  a1$activ <- ifelse(grepl("_AR_EGF", file) == T,
                     "AR_EGF",
                     ifelse(
                       grepl("_AR", file) == T,
                       "AR",
                       ifelse(grepl("_EGF", file) == T, "EGF",
                              ifelse(grepl("_00", file) == T, "00", NA))
                     ))
  WT_pheno <- rbind(WT_pheno, a1)
  rm(a1, res, optiona)
}
colnames(WT_pheno) <- gsub(".1", "", colnames(WT_pheno), fixed = T)

# drug loop ----
# single ----
druglist_pre <-
  c(
    "AKT",
    "AR",
    "AR_ERG",
    "Caspase8",
    "cFLAR",
    "EGFR",
    "ERK",
    "GLUT1",
    "HIF1",
    "HSPs",
    "MEK1_2",
    "MYC_MAX",
    "p14ARF",
    "PI3K",
    "ROS",
    "SPOP",
    "TERT"
  )
druglist <- paste(druglist_pre, collapse = "|")
drugs <-
  list.dirs(path = ".", recursive = F) %>% grep("LNCAP", ., value = TRUE) %>% .[grepl(druglist, .)]
cond <-
  "Proliferation,Angiogenesis,Invasion,Metastasis,Migration,Glycolysis,Hypermethylation,DNA_Repair,Quiescence,Apoptosis"
drug_pheno <- data.frame()
for (j in 1:length(drugs)) {
  print(drugs[j])
  file <-
    list.files(paste0(drugs[j], "/")) %>% grep("_probtraj_table.csv", ., value =
                                                 TRUE)
  a2 <- data.frame()
  for (i in 1:length(file)) {
    # print(file[i])
    folder_file <- paste0(drugs[j], "/", file[i])
    res <-
      read.table(
        folder_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      ) %>%
      .[nrow(.), , drop = FALSE]
    optiona <- list(
      input = res,
      outprefix = file[i] %>% sub("_probtraj_table.csv", "", .),
      conditions = cond
    )
    a1 <- cleaning_prob(optiona) %>% .[-ncol(.)]
    rownames(a1)[1] <- file[i] %>% sub("_probtraj_table.csv", "", .)
    a1$drug <- ifelse(
      grepl("_AR_EGF_", drugs[j]) == T,
      gsub("\\.\\/.*\\_AR\\_EGF\\_", "", drugs[j], perl = T),
      ifelse(
        grepl("_EGF_", drugs[j]) == T,
        gsub("\\.\\/.*\\_EGF\\_", "", drugs[j], perl = T),
        gsub("\\.\\/.*\\_00\\_|\\.\\/.*?\\_AR\\_", "", drugs[j], perl = T)
      )
    )
    a1$dose <- ifelse(grepl("_1_0", file[i]) == T,
                      "10",
                      ifelse(
                        grepl("_0_8", file[i]) == T,
                        "08",
                        ifelse(
                          grepl("_0_6", file[i]) == T,
                          "06",
                          ifelse(
                            grepl("_0_4", file[i]) == T,
                            "04",
                            ifelse(grepl("_0_2", file[i]) ==
                                     T, "02",
                                   ifelse(grepl("_0_0", file[i]) ==
                                            T, "00", NA))
                          )
                        )
                      ))
    a1$RNA <- ifelse(grepl("mutRNA", drugs[j]) == T,
                     "mutRNA",
                     ifelse(grepl("RNA", drugs[j]) == T, "RNA", NA))
    a1$activ <-
      ifelse(grepl(paste0(a1$RNA, "_AR_EGF_"), drugs[j]) == T,
             "AR_EGF",
             ifelse(
               grepl(paste0(a1$RNA, "_AR_"), drugs[j]) == T,
               "AR",
               ifelse(
                 grepl(paste0(a1$RNA, "_EGF_"), drugs[j]) == T,
                 "EGF",
                 ifelse(grepl(paste0(a1$RNA, "_00_"), drugs[j]) ==
                          T, "00", NA)
               )
             ))
    a2 <- rbind(a2, a1)
  }
  drug_pheno <- rbind(drug_pheno, a2)
  rm(a1, a2, res, optiona)
}
colnames(drug_pheno) <- gsub(".1", "", colnames(drug_pheno), fixed = T)
single_raw <- drug_pheno

single_noWT <- single_raw
single_noWT %<>% .[, -c(1)] %>% .[, -c(7)]

drug_pheno_norm2 <- merge(
  rbind(
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "00" &
                             drug_pheno$RNA == "RNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                              "00" & 
                                              WT_pheno$RNA == "RNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "00" &
                             drug_pheno$RNA == "mutRNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                              "00" & 
                                              WT_pheno$RNA == "mutRNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR_EGF" &
                             drug_pheno$RNA == "RNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                            "AR_EGF" & 
                                            WT_pheno$RNA == "RNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR_EGF" &
                             drug_pheno$RNA == "mutRNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                          "AR_EGF" & 
                                          WT_pheno$RNA == "mutRNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR" &
                             drug_pheno$RNA == "RNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                "AR" & 
                                                WT_pheno$RNA == "RNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR" &
                             drug_pheno$RNA == "mutRNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                              "AR" & 
                                              WT_pheno$RNA == "mutRNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "EGF" &
                             drug_pheno$RNA == "RNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                "EGF" & 
                                                WT_pheno$RNA == "RNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "EGF" &
                             drug_pheno$RNA == "mutRNA", ][, -c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                            "EGF" & 
                                            WT_pheno$RNA == "mutRNA", ][1, -c(1, 12:ncol(WT_pheno))]), `-`
    )))
  ),
  drug_pheno[, c(12:ncol(drug_pheno))],
  by = "row.names"
)

rownames(drug_pheno_norm2) <- drug_pheno_norm2$Row.names
drug_pheno_norm2 <-
  drug_pheno_norm2 %>% .[, -c(1)] %>% .[, -c(7)] #.[!colSums(is.na(.[,1:10])) > 0]
single <- drug_pheno_norm2

drug_pheno_large <- drug_pheno_norm2 %>%
  gather(., key = Phenotype, value, -c(drug, dose, RNA, activ)) %>%
  mutate(., dose = as.numeric(dose)) %>% mutate(., value2 = scale(.$value))

drug_pheno_large$activ <-
  factor(drug_pheno_large$activ, levels = c("00", "EGF", "AR", "AR_EGF"))
single_large <- drug_pheno_large

drug_pheno_large <- drug_pheno_norm2 %>%
  gather(., key = Phenotype, value, -c(drug, dose, RNA, activ)) %>%
  mutate(., dose = as.numeric(dose)) %>% mutate(., value2 = scale(.$value))

drug_pheno_large$activ <-
  factor(drug_pheno_large$activ, levels = c("00", "EGF", "AR", "AR_EGF"))
single_large_noWT <- drug_pheno_large

save(WT_pheno,
     single_large_noWT,
     single_noWT,
     single_large,
     single,
     single_raw,
     file = "./drugs_single.RData")
