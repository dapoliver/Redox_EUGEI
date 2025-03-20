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
library(ggrain)
library(predRupdate)
library(tibble)
library(Metrics)
library(rms)
library(predtools)
library(Hmisc)
library(readxl)
library(rstatix)

df <- read_csv("data_survival.csv")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

df_chr <- df
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)

df_chr <- df_chr %>% dplyr::rename(Gender = Gender.x)
df_cc <- df_chr %>% subset(select = c(Group, MIR132, MIR34A, MIR9, MIR941, MIR137, Age, Gender, Ethnicity, CAARMS, gafex01, gafex02, Transition_status, day_exit.x, site))
df_cc <- df_cc[complete.cases(df_cc), ]
df_cc <- df_cc %>% filter(MIR137 < 75)
df_cc <- df_cc %>% mutate(
  Gender = as.factor(Gender),
  Ethnicity = as.factor(Ethnicity)
)

data <- df_cc

MIR_9_plot <- ggplot(aes(x = Group, y = MIR9, fill = Group, colour = Group), data = data) +
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
  scale_y_continuous(name = "MIR-9") +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "HC")) +
  guides(fill = "none", color = "none")

MIR_34A_plot <- ggplot(aes(x = Group, y = MIR34A, fill = Group, colour = Group), data = data) +
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
  scale_y_continuous(name = "MIR-34A") +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "HC")) +
  guides(fill = "none", color = "none")

MIR_132_plot <- ggplot(aes(x = Group, y = MIR132, fill = Group, colour = Group), data = data) +
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
  scale_y_continuous(name = "MIR-132") +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "HC")) +
  guides(fill = "none", color = "none")

MIR_137_plot <- ggplot(aes(x = Group, y = MIR137, fill = Group, colour = Group), data = data) +
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
  scale_y_continuous(name = "MIR-137") +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "HC")) +
  guides(fill = "none", color = "none")

MIR_941_plot <- ggplot(aes(x = Group, y = MIR941, fill = Group, colour = Group), data = data) +
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
      max(data$MIR941, na.rm = TRUE) * 1.05, # First significance bar slightly above max
      max(data$MIR941, na.rm = TRUE) * 1.15
    ) # Second one higher up
  ) +
  theme_classic() +
  scale_color_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_fill_manual(values = c("#599ec4", "#c8526a", "gray80")) +
  scale_y_continuous(name = "MIR-941") +
  scale_x_discrete(labels = c("CHR-NT", "CHR-T", "HC")) +
  guides(fill = "none", color = "none")


combined_plot <- ggarrange(
  MIR_9_plot, NULL, MIR_34A_plot,
  MIR_132_plot, MIR_137_plot, MIR_941_plot, # NULL for the empty spot
  ncol = 3, nrow = 2
)
ggsave("Figure 4 wide.png", combined_plot, width = 20, height = 22, scale = 0.5)

shapiro.test(df_cc$MIR9)
shapiro.test(df_cc$MIR34A)
shapiro.test(df_cc$MIR132)
shapiro.test(df_cc$MIR137)
shapiro.test(df_cc$MIR941)


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
