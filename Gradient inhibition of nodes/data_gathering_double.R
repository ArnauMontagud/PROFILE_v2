---
title: "Gathering gradient inhibition simulations for double drugs"
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
    "RColorBrewer")
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
      !(grepl('Err', colnames(data))) &
      grepl('Prob', colnames(data))
    data_probCol_colnames <- colnames(data)[data_probCol]
    ## get columns of states fullfiling the conditions
    conditions_proba_list <-
      lapply(
        1:length(conditions),
        FUN = function(conditionN) {
          node <- unlist(strsplit(conditions[conditionN], '&'))
          value <-
            as.numeric(unlist(strsplit(states[conditionN], '&')))
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
  list.dirs(path = "./WT_sims", recursive = F) %>% grep("_00|_AR|_EGF", ., value =
                                                          TRUE) %>% sub("./WT/", "", .)

cond <-
  "Proliferation,Angiogenesis,Invasion,Metastasis,Migration,Glycolysis,Hypermethylation,DNA_Repair,Quiescence,Apoptosis"
WT_pheno <- data.frame()
for (j in 1:length(WT)) {
  folder <- WT[j]
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
# double ----

drugs <-
  list.dirs(path = "./doubles", recursive = F) %>% grep("LNCAP", ., value =
                                                          TRUE) %>% sub("./", "", .)
dosage <- c("_0_0", "_0_2", "_0_4", "_0_6", "_0_8", "_1_0")
cond <-
  "Proliferation,Angiogenesis,Invasion,Metastasis,Migration,Glycolysis,Hypermetilation,DNA_Repair,Quiescence,Apoptosis"
drugnames <-
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

drug_pheno <- data.frame()
for (j in 1:length(drugnames)) {
  for (k in 1:length(drugnames)) {
    if (k > j) {
      a3 <- data.frame()
      folder <-
        drugs[grepl(paste0("_", drugnames[j], "_", drugnames[k], "$"), drugs)]
      for (f in folder) {
        a2 <- data.frame()
        for (l in 1:length(dosage)) {
          for (m in 1:length(dosage)) {
            pattern <- paste0("_", drugnames[j], "_", drugnames[k], "$")
            prefix <-
              paste0(
                sub("doubles/drugsim_", "", sub(pattern, "", f, perl = T)),
                "_",
                drugnames[j],
                dosage[l],
                "_",
                drugnames[k],
                dosage[m]
              )
            print(prefix)
            res <-
              read.table(
                paste0("./", f, "/", prefix, "_probtraj_table.csv"),
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE
              ) %>% .[nrow(.), , drop = FALSE]
            optiona <- list(
              input = res,
              outprefix = prefix,
              conditions = cond
            )
            a1 <- cleaning_prob(optiona)
            rownames(a1)[1] <- prefix
            a1$drug1 <- drugnames[j]
            a1$drug2 <- drugnames[k]
            a1$dose1 <-
              as.numeric(sub("_", ".", sub("^_", "", dosage[l], perl = T)))
            a1$dose2 <-
              as.numeric(sub("_", ".", sub("^_", "", dosage[m], perl = T)))
            a1$RNA <- ifelse(grepl("mutRNA", prefix) == T,
                             "mutRNA",
                             ifelse(grepl("RNA", prefix) == T, "RNA", NA))
            a1$activ <-
              ifelse(grepl(paste0(a1$RNA, "_AR_EGF_"), prefix) == T,
                     "AR_EGF",
                     ifelse(
                       grepl(paste0(a1$RNA, "_AR_"), prefix) == T,
                       "AR",
                       ifelse(
                         grepl(paste0(a1$RNA, "_EGF_"), prefix) == T,
                         "EGF",
                         ifelse(grepl(
                           paste0(a1$RNA, "_00_"), prefix
                         ) == T, "00", NA)
                       )
                     ))
            
            a2 <- rbind(a2, a1)
          }
        }
        a3 <- rbind(a3, a2)
      }
      drug_pheno <- rbind(drug_pheno, a3)
      rm(a1, a2, a3, res, optiona)
    } else{
      next
    }
  }
}

colnames(drug_pheno) <-
  gsub(".1", "", colnames(drug_pheno), fixed = T)
double_raw <- drug_pheno

drug_pheno_norm2 <- merge(
  rbind(
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "00" &
                             drug_pheno$RNA == "RNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                          "00" &
                                                          WT_pheno$RNA == "RNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "00" &
                             drug_pheno$RNA == "mutRNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                      "00" &
                                                      WT_pheno$RNA == "mutRNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR_EGF" &
                             drug_pheno$RNA == "RNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                          "AR_EGF" &
                                                          WT_pheno$RNA == "RNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR_EGF" &
                             drug_pheno$RNA == "mutRNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                      "AR_EGF" &
                                                      WT_pheno$RNA == "mutRNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR" &
                             drug_pheno$RNA == "RNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                          "AR" &
                                                          WT_pheno$RNA == "RNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "AR" &
                             drug_pheno$RNA == "mutRNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                      "AR" &
                                                      WT_pheno$RNA == "mutRNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "EGF" &
                             drug_pheno$RNA == "RNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                          "EGF" &
                                                          WT_pheno$RNA == "RNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    ))),
    as.data.frame((sweep(
      as.matrix(drug_pheno[drug_pheno$activ == "EGF" &
                             drug_pheno$RNA == "mutRNA",][,-c(1, 12:ncol(drug_pheno))]), 2, as.matrix(WT_pheno[WT_pheno$activ ==
                                                      "EGF" &
                                                      WT_pheno$RNA == "mutRNA",][1,-c(1, 12:ncol(WT_pheno))]), `-`
    )))
  ),
  drug_pheno[, c(12:ncol(drug_pheno))],
  by = "row.names"
)


rownames(drug_pheno_norm2) <- drug_pheno_norm2$Row.names
drug_pheno_norm2 <- drug_pheno_norm2 %>% .[,-c(1)]
double <- drug_pheno_norm2

drug_pheno_large <- drug_pheno_norm2 %>%
  gather(., key = Phenotype, value,-c(drug1, drug2, dose1, dose2, RNA, activ)) %>%
  mutate(., dose1 = as.numeric(dose1), dose2 = as.numeric(dose2)) %>% mutate(., value2 =
                                                                               scale(.$value))

drug_pheno_large$activ <-
  factor(drug_pheno_large$activ, levels = c("00", "EGF", "AR", "AR_EGF"))
double_large <- drug_pheno_large

save(WT_pheno, 
    double, 
    double_raw, 
    double_large, 
    file = "./drugs_double.RData")
