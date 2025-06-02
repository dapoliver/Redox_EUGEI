setwd("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Redox EU-GEI/")

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(glmnet)
library(caret)
library(pROC)
library(dcurves)
library(probably)
library(ggplot2)
library(ggsignif)
library(ggrain)
library(ggpubr)
library(predRupdate)
library(tibble)
library(Metrics)
library(rms)
library(predtools)
library(Hmisc)
library(readxl)
library(rstatix)
library(gtsummary)

df <- read_csv("data_survival.csv")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

df_chr <- df
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)
df_chr$GAF <- rowMeans(df_chr[, c("gafex01", "gafex02")], na.rm = TRUE)

df_chr <- df_chr %>% dplyr::rename(Gender = Gender.x)
tbl <- tbl_summary(
  include = c(Age, Gender, Ethnicity, CAARMS, GAF), data = df_chr, by = Group,
  type = list(
    Age ~ "continuous2",
    CAARMS ~ "continuous2",
    GAF ~ "continuous2"
  ),
  statistic = list(all_continuous() ~ c("{mean}", "{sd}")),
  digits = list(
    all_continuous() ~ c(1, 1),
    all_categorical() ~ c(0, 1)
  )
)

df_cc <- df_chr %>% subset(select = c(Group, MIR132, MIR34A, MIR9, MIR941, MIR137))
df_cc <- df_cc[complete.cases(df_cc), ]
df_cc <- df_cc %>% filter(MIR137 < 75)

data <- df_cc
data$study <- "EU-GEI"

shapiro.test(df_cc$MIR9)
shapiro.test(df_cc$MIR34A)
shapiro.test(df_cc$MIR132)
shapiro.test(df_cc$MIR137)
shapiro.test(df_cc$MIR941)

univ.summary <- data.frame(
  miRNA = c(
    "miR9", "miR34A", "miR132", "miR137", "miR941",
    "miR9", "miR34A", "miR132", "miR137", "miR941",
    "miR9", "miR34A", "miR132", "miR137", "miR941"
  ),
  contrast = c(
    "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT",
    "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls",
    "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls"
  ),
  p.value = c(
    wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value
  ),
  eff.size = c(
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")")
  ),
  magnitude = c(
    paste(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude)
  )
)

univ.summary$p.value <- p.adjust(univ.summary$p.value, method = "BH")
write_csv(univ.summary, "univariate_summary_EUGEI.csv")

wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])
wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])
wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])
wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])
wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])

wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])
wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])
wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])
wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])
wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ])

wilcox.test(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])
wilcox.test(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])
wilcox.test(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])
wilcox.test(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])
wilcox.test(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ])

wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)
wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)
wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)
wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)
wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)

wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)
wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)
wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)
wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)
wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)

wilcox_effsize(MIR9 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)
wilcox_effsize(MIR34A ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)
wilcox_effsize(MIR132 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)
wilcox_effsize(MIR137 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)
wilcox_effsize(MIR941 ~ Group, data = df_cc[df_cc$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)

# Load NAPLS data
df_NAPLS <- read_csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Redox EU-GEI/NAPLS/NAPLS.csv")

df_NAPLS <- df_NAPLS %>%
  mutate(Group = case_when(
    `GROUP (UC = CTRL group)` == "UC" ~ "Control",
    `GROUP (UC = CTRL group)` == "CHR-C" ~ "AtRisk_Trans",
    `GROUP (UC = CTRL group)` == "CHR-NC" ~ "AtRisk_NoTr"
  )) %>%
  rename(MIR9 = `miR-9`, MIR34A = `miR-34`, MIR132 = `miR-132`, MIR137 = `miR-137`, MIR941 = `miR-941`) %>%
  subset(select = c(
    MIR9, MIR34A, MIR132, MIR137, MIR941, Group
  ))

df_NAPLS <- df_NAPLS %>% filter(MIR941 < 100)

shapiro.test(df_NAPLS$MIR9)
shapiro.test(df_NAPLS$MIR34A)
shapiro.test(df_NAPLS$MIR132)
shapiro.test(df_NAPLS$MIR137)
shapiro.test(df_NAPLS$MIR941)


univ.summary <- data.frame(
  miRNA = c(
    "miR9", "miR34A", "miR132", "miR137", "miR941",
    "miR9", "miR34A", "miR132", "miR137", "miR941",
    "miR9", "miR34A", "miR132", "miR137", "miR941"
  ),
  contrast = c(
    "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT", "CHR-T vs CHR-NT",
    "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls", "CHR-T vs Controls",
    "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls", "CHR-NT vs Controls"
  ),
  p.value = c(
    wilcox.test(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$p.value,
    wilcox.test(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$p.value,
    wilcox.test(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value,
    wilcox.test(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$p.value
  ),
  eff.size = c(
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")"),
    paste0(round(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$effsize, 2), " (", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.low, "-", wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ], ci = TRUE)$conf.high, ")")
  ),
  magnitude = c(
    paste(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "AtRisk_NoTr"), ])$magnitude),
    paste(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_Trans", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR9 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR34A ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR132 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR137 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude),
    paste(wilcox_effsize(MIR941 ~ Group, data = df_NAPLS[df_NAPLS$Group %in% c("AtRisk_NoTr", "Control"), ])$magnitude)
  )
)

univ.summary$p.value <- p.adjust(univ.summary$p.value, method = "BH")
write_csv(univ.summary, "univariate_summary.csv")

df_NAPLS$study <- factor(c("NAPLS-3"), levels = c("EU-GEI", "NAPLS-3"))
data$study <- factor(c("EU-GEI"), levels = c("EU-GEI", "NAPLS-3"))

MIR_9_plot <- ggplot(aes(x = Group, y = MIR9, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("***", "**"),
    tip_length = 0.02
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-9", limits = c(0, max(data$MIR9, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_34A_plot <- ggplot(aes(x = Group, y = MIR34A, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("***", "***"),
    tip_length = 0.02
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-34a", limits = c(0, max(data$MIR34A, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_132_plot <- ggplot(aes(x = Group, y = MIR132, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("**", "***"),
    tip_length = 0.02,
    y_position = c(
      max(data$MIR132, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(data$MIR132, na.rm = TRUE) * 1.15
    ) # Second one higher up
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-132", limits = c(0, max(data$MIR132, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_137_plot <- ggplot(aes(x = Group, y = MIR137, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("***", "***"),
    tip_length = 0.02,
    y_position = c(
      max(data$MIR137, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(data$MIR137, na.rm = TRUE) * 1.15
    )
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-137", limits = c(0, max(data$MIR137, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_941_plot <- ggplot(aes(x = Group, y = MIR941, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("**", "***"),
    tip_length = 0.02,
    y_position = c(
      max(data$MIR941, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(data$MIR941, na.rm = TRUE) * 1.15
    ) # Second one higher up
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-941", limits = c(0, max(data$MIR941, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_9_NAPLS_plot <- ggplot(aes(x = Group, y = MIR9, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("***", "***"),
    tip_length = 0.02
  ) +
  theme_classic() +
  scale_color_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-9", limits = c(0, max(df_NAPLS$MIR9, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_34A_NAPLS_plot <- ggplot(aes(x = Group, y = MIR34A, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("***", "**"),
    tip_length = 0.02
  ) +
  theme_classic() +
  scale_color_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-34a", limits = c(0, max(df_NAPLS$MIR34A, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_132_NAPLS_plot <- ggplot(aes(x = Group, y = MIR132, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("***", "***"),
    tip_length = 0.02,
    y_position = c(
      max(df_NAPLS$MIR132, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(df_NAPLS$MIR132, na.rm = TRUE) * 1.15
    ) # Second one higher up
  ) +
  theme_classic() +
  scale_color_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-132", limits = c(0, max(df_NAPLS$MIR132, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_137_NAPLS_plot <- ggplot(aes(x = Group, y = MIR137, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .2, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "Control"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("**", "**"),
    tip_length = 0.02,
    y_position = c(
      max(df_NAPLS$MIR137, na.rm = TRUE) * 1.15, # First significance bar slightly above max
      max(df_NAPLS$MIR137, na.rm = TRUE) * 1.05
    )
  ) +
  theme_classic() +
  scale_color_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-137", limits = c(0, max(df_NAPLS$MIR137, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_941_NAPLS_plot <- ggplot(aes(x = Group, y = MIR941, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("***", "***"),
    tip_length = 0.02,
    y_position = c(
      max(df_NAPLS$MIR941, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(df_NAPLS$MIR941, na.rm = TRUE) * 1.15
    ) # Second one higher up
  ) +
  theme_classic() +
  scale_color_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#c8526a", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-941", limits = c(0, max(df_NAPLS$MIR941, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

combined_plot <- ggarrange(
  MIR_9_plot, NULL, MIR_34A_plot,
  MIR_9_NAPLS_plot, NULL, MIR_34A_NAPLS_plot,
  MIR_132_plot, MIR_137_plot, MIR_941_plot, # NULL for the empty spot
  MIR_132_NAPLS_plot, MIR_137_NAPLS_plot, MIR_941_NAPLS_plot, # NULL for the empty spot
  ncol = 3, nrow = 4
)
ggsave("Figure 4 wide 220525.png", combined_plot, width = 20, height = 22, scale = 0.5)
