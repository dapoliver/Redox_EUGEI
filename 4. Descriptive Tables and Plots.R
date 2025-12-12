setwd("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/")

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

# df <- read_csv("data_survival.csv")
df <- read_excel("Eugei vf 0810 final BATCH 091025.xlsx")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

##### Table 1 #####

df_chr <- df
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)
df_chr$GAF <- rowMeans(df_chr[, c("gafex01", "gafex02")], na.rm = TRUE)

df_chr <- df_chr %>% dplyr::rename(Gender = Gender.x)
df_chr <- df_chr %>% filter(!is.na(MIR132))

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

##### Results Summary Plot #####
results_plot <- data.frame(
  Type = c(
    "CHR-T v CHR-NT",
    "CHR-T v CHR-NT",
    "CHR-T v Controls",
    "CHR-T v Controls"
  ),
  sample = c(
    "EU-GEI",
    "NAPLS-3",
    "EU-GEI",
    "NAPLS-3"
  ),
  C = c(
    0.933,
    1,
    0.908,
    0.868
  ),
  lCI = c(
    0.919,
    0.983,
    0.876,
    0.850
  ),
  uCI = c(
    0.947,
    1,
    0.940,
    0.887
  )
)
results_plot$Type <- factor(results_plot$Type, levels = c("CHR-T v CHR-NT", "CHR-T v Controls"))
results_plot$sample <- factor(results_plot$sample, levels = c("EU-GEI", "NAPLS-3"))

ggplot(data = results_plot, aes(x = Type, group = sample)) +
  geom_rect(
    data = results_plot, aes(), xmin = 0, xmax = 9.5, ymin = 0.5, ymax = 0.7,
    fill = "gray92"
  ) +
  geom_rect(
    data = results_plot, aes(), xmin = 0, xmax = 9.5, ymin = 0.8, ymax = 0.9,
    fill = "gray92"
  ) +
  geom_text(aes(x = 1.5, y = 0.6, label = "Above chance"), stat = "unique", size = 8, color = "gray80", family = "Roboto Condensed") +
  geom_text(aes(x = 1.5, y = 0.45, label = "Below chance"), stat = "unique", size = 8, color = "gray80", family = "Roboto Condensed") +
  geom_text(aes(x = 1.5, y = 0.75, label = "Acceptable"), stat = "unique", size = 8, color = "gray80", family = "Roboto Condensed") +
  geom_text(aes(x = 1.5, y = 0.85, label = "Excellent"), stat = "unique", size = 8, color = "gray80", family = "Roboto Condensed") +
  geom_text(aes(x = 1.5, y = 0.95, label = "Outstanding"), stat = "unique", size = 8, color = "gray80", family = "Roboto Condensed") +
  geom_text(aes(x = 0.75, y = 1.05, label = "0.93"), stat = "unique", size = 8, color = "#599ec4", family = "Roboto Condensed") +
  geom_text(aes(x = 1.25, y = 1.05, label = "1.00"), stat = "unique", size = 8, color = "#c8526a", family = "Roboto Condensed") +
  geom_text(aes(x = 1.75, y = 1.05, label = "0.91"), stat = "unique", size = 8, color = "#599ec4", family = "Roboto Condensed") +
  geom_text(aes(x = 2.25, y = 1.05, label = "0.87"), stat = "unique", size = 8, color = "#c8526a", family = "Roboto Condensed") +
  geom_pointrange(data = results_plot, mapping = aes(x = Type, y = C, ymin = lCI, ymax = uCI, color = sample), size = 2, fatten = 2, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("#599ec4", "#c8526a")) +
  theme_classic() +
  xlab("Model") +
  ylab("C-index") +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
  guides(color = guide_legend(title = "Sample")) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 21), legend.position = "right", legend.title = element_text(size = 23), legend.text = element_text(size = 23))
ggsave("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/Figure 2 111225.png", width = 42, height = 32, units = "cm", scale = 0.65)

##### Univariate analyses #####

df_cc <- df_chr %>%
  subset(select = c(Group, MIR132, MIR34A, MIR9, MIR941, MIR137))
df_cc <- df_cc[complete.cases(df_cc), ]
# df_cc <- df_cc %>% filter(MIR137 < 75)

data <- df_cc
data$study <- "EU-GEI"

shapiro.test(df_cc$MIR9)
shapiro.test(df_cc$MIR34A)
shapiro.test(df_cc$MIR132)
shapiro.test(df_cc$MIR137)
shapiro.test(df_cc$MIR941)

univ.summary <- data.frame(
  study = "EU-GEI",
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
write_csv(univ.summary, "univariate_summary_EUGEI_111025.csv")

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
df_NAPLS <- read_excel("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/NAPLS/NAPLS BATCH info17102025.xlsx")

df_NAPLS <- df_NAPLS %>%
  mutate(
    Group = case_when(
      `GROUP (UC = CTRL group)` == "UC" ~ "Control",
      `GROUP (UC = CTRL group)` == "CHR-C" ~ "AtRisk_Trans",
      `GROUP (UC = CTRL group)` == "CHR-NC" ~ "AtRisk_NoTr"
    ),
    P2_P3_SOPS = case_when(
      P2_SOPS > P3_SOPS ~ P2_SOPS,
      TRUE ~ P3_SOPS
    ),
    P1_CAARMS = case_when(
      P1_SOPS == 0 ~ 0.011,
      P1_SOPS == 1 ~ 0.916,
      P1_SOPS == 2 ~ 1.816,
      P1_SOPS == 3 ~ 2.759,
      P1_SOPS == 4 ~ 3.799,
      P1_SOPS == 5 ~ 4.976,
      P1_SOPS == 6 ~ 6.033
    ),
    P2_CAARMS = case_when(
      P2_P3_SOPS == 0 ~ 0.007,
      P2_P3_SOPS == 1 ~ 0.919,
      P2_P3_SOPS == 2 ~ 1.778,
      P2_P3_SOPS == 3 ~ 2.735,
      P2_P3_SOPS == 4 ~ 3.806,
      P2_P3_SOPS == 5 ~ 4.991,
      P2_P3_SOPS == 6 ~ 6.025
    ),
    P3_CAARMS = case_when(
      P4_SOPS == 0 ~ 0.013,
      P4_SOPS == 1 ~ 1.112,
      P4_SOPS == 2 ~ 2.106,
      P4_SOPS == 3 ~ 3.045,
      P4_SOPS == 4 ~ 4.059,
      P4_SOPS == 5 ~ 5.153,
      P4_SOPS == 6 ~ 6.099
    ),
    P4_CAARMS = case_when(
      P5_SOPS == 0 ~ 0.079,
      P5_SOPS == 1 ~ 1.126,
      P5_SOPS == 2 ~ 2.017,
      P5_SOPS == 3 ~ 2.968,
      P5_SOPS == 4 ~ 3.936,
      P5_SOPS == 5 ~ 4.844,
      P5_SOPS == 6 ~ 5.889
    ),
    CAARMS = P1_CAARMS + P2_CAARMS + P3_CAARMS + P4_CAARMS,
    Transition = ifelse(`GROUP (UC = CTRL group)` == "CHR-C", 1, 0),
    Ethnicity = case_when(
      demo_racial == "European" ~ "White",
      demo_racial == "African" ~ "Black",
      demo_racial == "East Asian" | demo_racial == "South Asian" ~ "Asian",
      demo_racial == "Interracial" ~ "Mixed",
      TRUE ~ "Other"
    )
  ) %>%
  rename(MIR9 = `miR-9`, MIR34A = `miR-34`, MIR132 = `miR-132`, MIR137 = `miR-137`, MIR941 = `miR-941`, day_exit = fudays) %>%
  subset(select = c(
    MIR9, MIR34A, MIR132, MIR137, MIR941, Group, demo_age_ym, demo_sex, Ethnicity, CAARMS, GlobalAssessmentFunction, day_exit
  ))

tbl_NAPLS <- tbl_summary(
  include = c(demo_age_ym, demo_sex, Ethnicity, CAARMS, GlobalAssessmentFunction), data = df_NAPLS, by = Group,
  type = list(
    demo_age_ym ~ "continuous2",
    CAARMS ~ "continuous2",
    GlobalAssessmentFunction ~ "continuous2"
  ),
  statistic = list(all_continuous() ~ c("{mean}", "{sd}")),
  digits = list(
    all_continuous() ~ c(1, 1),
    all_categorical() ~ c(0, 1)
  )
)

shapiro.test(df_NAPLS$MIR9)
shapiro.test(df_NAPLS$MIR34A)
shapiro.test(df_NAPLS$MIR132)
shapiro.test(df_NAPLS$MIR137)
shapiro.test(df_NAPLS$MIR941)


univ.summary_NAPLS <- data.frame(
  study = "NAPLS-3",
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
write_csv(univ.summary, "univariate_summary_NAPLS.csv")

##### Descriptive Figures #####

df_NAPLS$study <- factor(c("NAPLS-3"), levels = c("EU-GEI", "NAPLS-3"))
data$study <- factor(c("EU-GEI"), levels = c("EU-GEI", "NAPLS-3"))

MIR_9_plot <- ggplot(aes(x = Group, y = MIR9, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
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
  scale_y_continuous(name = "MIR-9", limits = c(0, max(data$MIR9, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_34A_plot <- ggplot(aes(x = Group, y = MIR34A, fill = study, colour = study), data = data) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
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
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("***", "***"),
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
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
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
    annotations = c("***", "***"),
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
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
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
  scale_y_continuous(name = "MIR-9", limits = c(0, max(df_NAPLS$MIR9, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_34A_NAPLS_plot <- ggplot(aes(x = Group, y = MIR34A, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
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
  scale_y_continuous(name = "MIR-34a", limits = c(0, max(df_NAPLS$MIR34A, na.rm = TRUE) * 1.25)) +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "Controls")) +
  guides(fill = "none", color = "none")

MIR_132_NAPLS_plot <- ggplot(aes(x = Group, y = MIR132, fill = study, colour = study), data = df_NAPLS) +
  geom_rain(
    alpha = .5, rain.side = "l",
    boxplot.args = list(color = "black", outlier.shape = NA),
    boxplot.args.pos = list(
      position = ggpp::position_dodgenudge(x = .15, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "AtRisk_Trans"), c("AtRisk_Trans", "Control"), c("AtRisk_NoTr", "Control")),
    colour = "black",
    annotations = c("***", "*", "***"),
    tip_length = 0.02,
    y_position = c(
      max(df_NAPLS$MIR132, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(df_NAPLS$MIR132, na.rm = TRUE) * 1.05,
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
      position = ggpp::position_dodgenudge(x = 0.15, width = 0.2), width = 0.15
    )
  ) +
  geom_signif(
    comparisons = list(c("AtRisk_NoTr", "Control"), c("AtRisk_Trans", "Control")),
    colour = "black",
    annotations = c("**", "***"),
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
ggsave("Figure 4 wide 171025.png", combined_plot, width = 20, height = 22, scale = 0.5)

##### Descriptive PI plots #####
df_NAPLS <- df_NAPLS %>% mutate(PI_CHR = -1.4588309 + (-0.7308739 * MIR9) + (1.02934979 * MIR34A) + (-0.2071588 * MIR132) + (0 * MIR137) + (-1.1861267 * MIR941))
df_NAPLS$risk <- 1 / (1 + exp(-df_NAPLS$PI_CHR))
df_NAPLS_chr <- df_NAPLS %>% filter(`GROUP (UC = CTRL group)` != "UC" & !is.na(risk))
df_NAPLS_chr <- df_NAPLS_chr %>% mutate(
  Transition = case_when(`GROUP (UC = CTRL group)` == "CHR-C" ~ 1, TRUE ~ 0),
  Group = case_when(`GROUP (UC = CTRL group)` == "CHR-C" ~ "CHR-T", TRUE ~ "CHR-NT")
)
recal_model <- glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit"))
df_NAPLS_chr$recalibrated_probs <- predict(recal_model, type = "response")

ggplot(data = df_NAPLS_chr, aes(x = recalibrated_probs * 100, fill = Group)) +
  geom_histogram() +
  scale_fill_manual(
    values = c("#c8526a", "#599ec4"),
    labels = c("CHR-NT", "CHR-T")
  ) +
  labs(x = "Risk (%)", y = "Frequency") +
  theme_classic() +
  theme(legend.position = "top")
ggsave("Risk Histogram 171025.png", width = 20, height = 22, scale = 0.5)

##### KM plot ######
data_all <- rbind(data, df_NAPLS)
data_all <- data_all %>%
  filter(Group != "Control") %>%
  mutate(Transition_status = case_when(Group == "AtRisk_Trans" ~ 1, Group == "AtRisk_NoTr" ~ 0))

pdf("KM.pdf")
ggsurvplot(survfit(Surv(day_exit, Transition_status) ~ study, data = data_all),
  fun = "event",
  xscale = 365.25,
  break.time.by = 182.625,
  xlab = "Time (years)",
  palette = c("#599ec4", "#c8526a"),
  risk.table = FALSE,
  cumevents = TRUE,
  conf.int = FALSE
)
dev.off()
ggsave("KM.png", width = 42, height = 32, units = "cm", scale = 0.65)

##### Correlation matrix ####
library(corrplot)

# Compute correlation matrix
corr_matrix <- cor(data[, 2:6])
testRes <- cor.mtest(data[, 2:6], conf.level = 0.95)

# Draw the heatmap
corrplot(corr_matrix,
  method = "color", type = "upper", addCoef.col = "black",
  tl.col = "black", tl.srt = 45, p.mat = testRes$p, sig.level = 0.10
)
