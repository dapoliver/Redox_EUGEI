setwd("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Redox EU-GEI/NAPLS/")

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(glmnet)
library(caret)
library(pROC)
library(dcurves)
library(probably)
library(ggplot2)
library(predRupdate)
library(tibble)
library(Metrics)
library(rms)
library(predtools)
library(Hmisc)
library(readxl)
library(survival)
library(predRupdate)
library(gtsummary)

df <- read_csv("NAPLS.csv")

summary(factor(df$`GROUP (UC = CTRL group)`))
df <- df %>% subset(select = c(
  demo_age_ym, demo_sex, demo_racial, `miR-9`, `miR-34`, `miR-132`, `miR-137`, `miR-941`, `GROUP (UC = CTRL group)`, SiteNumber,
  P1_SOPS, P2_SOPS, P3_SOPS, P4_SOPS, P5_SOPS, GlobalAssessmentFunction
))
df <- df %>%
  mutate(
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
    CAARMS_total = P1_CAARMS + P2_CAARMS + P3_CAARMS + P4_CAARMS,
    PI_CHR = 2.74417757 + (-2.685152 * `miR-9`) + (6.99712917 * `miR-34`) + (0.60454662 * `miR-132`) + (-0.447368 * `miR-137`) + (-15.028186 * `miR-941`),
    PI_HC = -2.0418329 + (-0.2277153 * `miR-9`) + (0.751332827 * `miR-34`),
    Transition = ifelse(`GROUP (UC = CTRL group)` == "CHR-C", 1, 0),
    ethnicity = case_when(
      demo_racial == "European" ~ "White",
      demo_racial == "African" ~ "Black",
      demo_racial == "East Asian" | demo_racial == "South Asian" ~ "Asian",
      demo_racial == "Interracial" ~ "Mixed",
      TRUE ~ "Other"
    )
  )
df <- df[complete.cases(df), ]
df <- df %>% filter(`miR-941` < 100)

tbl_summary(df,
  by = `GROUP (UC = CTRL group)`, include = c(demo_age_ym, demo_sex, ethnicity, CAARMS_total, GlobalAssessmentFunction), missing = "no", type = list(
    demo_age_ym ~ "continuous2",
    CAARMS_total ~ "continuous2",
    GlobalAssessmentFunction ~ "continuous2"
  ),
  statistic = list(all_continuous() ~ c("{mean}", "{sd}")),
  digits = list(
    all_continuous() ~ c(1, 1),
    all_categorical() ~ c(0, 1)
  )
) %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Group**") %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_flex_table()

##### CHR-T v CHR-NT #####

df_chr <- df %>% filter(`GROUP (UC = CTRL group)` != "UC")

concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))

df_chr$risk <- 1 / (1 + exp(-df_chr$PI_CHR))

predicted_labels <- ifelse(df_chr$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_chr$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

results_NAPLS <- data.frame(
  C = paste0(
    round(concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))$concordance, 3), " (",
    round(concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))$concordance - 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))$cvar), 3), "-",
    round(concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))$concordance + 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit")))$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm$byClass[4] * 100, 1), "%"
  )
)

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_chr$risk
)

calibration <- calibration %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)

cal_plot_breaks(calibration, truth = observed, estimate = predicted)

recal_model <- glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

results_NAPLS_recal <- data.frame(
  C = paste0(
    round(concordance(recal_model)$concordance, 3), " (",
    round(concordance(recal_model)$concordance - 1.96 * (concordance(recal_model)$cvar), 3), "-",
    round(concordance(recal_model)$concordance + 1.96 * (concordance(recal_model)$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm_recal$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm_recal$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm_recal$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm_recal$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm_recal$byClass[4] * 100, 1), "%"
  )
)
# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = predicted_labels_recal
)

calibration_recal <- calibration_recal %>%
  rowwise() %>%
  mutate(
    predicted = case_when(
      predicted > 0.999 ~ sample(c(0.999, 0.9991, 0.9992, 0.9993, 0.9994, 0.9995, 0.9996, 0.9997, 0.9998, 0.9999), size = 1),
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  ) %>%
  ungroup()
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted, cal_plot = FALSE)

cal_plot_breaks(calibration_recal, truth = observed, estimate = predicted)

# Mean Offset
### Compute Global Means ###
global_mean <- colMeans(df_chr[, 4:8, drop = FALSE], na.rm = TRUE)

### Mean Offset Correction ###
batch_test <- as.factor(df_chr$SiteNumber)
df_chr_corrected <- df_chr # Start with the original data

for (b in levels(batch_test)) {
  batch_indices <- which(batch_test == b) # Indices for samples in batch `b`

  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }

  # Extract the subset of test for the current batch
  batch_data <- df_chr[batch_indices, 4:8, drop = FALSE] # Exclude site and outcome

  # Compute means for each predictor in this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }
  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  # Subtract the offset to align the batch with the global mean
  df_chr_corrected[batch_indices, 4:8] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
  paste("Batch", b, "corrected.")
}

df_chr_corrected <- df_chr_corrected %>%
  mutate(
    PI_CHR = 2.74417757 + (-2.685152 * `miR-9`) + (6.99712917 * `miR-34`) + (0.60454662 * `miR-132`) + (-0.447368 * `miR-137`) + (-15.028186 * `miR-941`)
  )
concordance(glm(Transition ~ PI_CHR, data = df_chr_corrected, family = binomial(link = "logit")))

df_chr_corrected$risk <- 1 / (1 + exp(-df_chr_corrected$PI_CHR))

predicted_labels <- ifelse(df_chr_corrected$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_chr_corrected$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_chr_corrected$risk
)

calibration <- calibration %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)

recal_model <- glm(Transition ~ PI_CHR, data = df_chr_corrected, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = predicted_labels_recal
)

calibration_recal <- calibration_recal %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted, cal_plot = FALSE)

##### CHR-T v HC #####

df_hc <- df %>% filter(`GROUP (UC = CTRL group)` != "CHR-NC")

concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))

df_hc$risk <- 1 / (1 + exp(-df_hc$PI_HC))

predicted_labels <- ifelse(df_hc$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_hc$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

results_NAPLS_hc <- data.frame(
  C = paste0(
    round(concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))$concordance, 3), " (",
    round(concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))$concordance - 1.96 * (concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))$cvar), 3), "-",
    round(concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))$concordance + 1.96 * (concordance(glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit")))$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm$byClass[4] * 100, 1), "%"
  )
)

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_hc$risk
)

calibration <- calibration %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)
cal_plot_breaks(calibration, truth = observed, estimate = predicted)

recal_model <- glm(Transition ~ PI_HC, data = df_hc, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

results_NAPLS_recal <- data.frame(
  C = paste0(
    round(concordance(recal_model)$concordance, 3), " (",
    round(concordance(recal_model)$concordance - 1.96 * (concordance(recal_model)$cvar), 3), "-",
    round(concordance(recal_model)$concordance + 1.96 * (concordance(recal_model)$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm_recal$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm_recal$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm_recal$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm_recal$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm_recal$byClass[4] * 100, 1), "%"
  )
)

# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = predicted_labels_recal
)

calibration_recal <- calibration_recal %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
calibration_recal <- calibration_recal %>% mutate(truth = case_when(observed == 1 ~ "Transition", observed == 0 ~ "No Transition"))
calibration_recal$truth <- factor(calibration_recal$truth, levels = c("No Transition", "Transition"))
calibration_recal$truth <- factor(calibration_recal$observed, levels = c(0, 1))
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted)
cal_plot_breaks(calibration_recal, truth = observed, estimate = predicted)


##### Decision Curve Analysis #####
dca_recal_chr <- dca(observed ~ predicted,
  data = calibration_recal,
  prevalence = 0.22,
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

dca_recal_chr <- dca_recal_chr %>% filter(label == "predicted")
dca_recal_chr$label <- "NAPLS"
dca_recal_chr$variable <- "NAPLS"

dca_recal_hc <- dca(observed ~ predicted,
  data = calibration_recal,
  prevalence = 0.000266,
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

dca_recal_hc <- dca_recal_hc %>% filter(label == "predicted")
dca_recal_hc$label <- "NAPLS"
dca_recal_hc$variable <- "NAPLS"

##### Mean Offset #####
batch <- as.factor(df_chr$SiteNumber)
df_cc_corrected <- df_chr # Start with the original data

### Compute Global Means ###
global_mean <- colMeans(df_chr[, 4:8, drop = FALSE], na.rm = TRUE)

# Iterate over each batch
for (b in levels(batch)) {
  batch_indices <- which(batch == b) # Indices for samples in batch `b`

  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }

  # Extract the subset of test for the current batch
  batch_data <- df_chr[batch_indices, 4:8, drop = FALSE]

  # Compute row-wise means for this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }

  # Subtract batch means
  df_cc_corrected[batch_indices, 4:8] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
}

df_cc_corrected <- df_cc_corrected %>%
  mutate(
    PI_CHR = -2.6723787 + (-0.2158028 * `miR-9`) + (1.1518603 * `miR-34`) + (-0.1884227 * `miR-132`) + (-0.1594170 * `miR-137`) + (-0.9113908 * `miR-941`)
  )

concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))

df_cc_corrected$risk <- 1 / (1 + exp(-df_cc_corrected$PI_CHR))

predicted_labels <- ifelse(df_cc_corrected$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_cc_corrected$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

results_NAPLS <- data.frame(
  C = paste0(
    round(concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))$concordance, 3), " (",
    round(concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))$concordance - 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))$cvar), 3), "-",
    round(concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))$concordance + 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit")))$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm$byClass[4] * 100, 1), "%"
  )
)

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_cc_corrected$risk
)

calibration <- calibration %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)

cal_plot_breaks(calibration, truth = observed, estimate = predicted)

recal_model <- glm(Transition ~ PI_CHR, data = df_cc_corrected, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = predicted_labels_recal
)

calibration_recal <- calibration_recal %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted, cal_plot = FALSE)

##### Mean Offset HC #####
batch <- as.factor(df_hc$SiteNumber)
df_hc_corrected <- df_hc # Start with the original data

### Compute Global Means ###
global_mean <- colMeans(df_hc[, 4:8, drop = FALSE], na.rm = TRUE)

# Iterate over each batch
for (b in levels(batch)) {
  batch_indices <- which(batch == b) # Indices for samples in batch `b`

  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }

  # Extract the subset of test for the current batch
  batch_data <- df_hc[batch_indices, 4:8, drop = FALSE]

  # Compute row-wise means for this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }

  # Subtract batch means
  df_hc_corrected[batch_indices, 4:8] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
}

df_hc_corrected <- df_hc_corrected %>%
  mutate(
    PI_HC = -2.4362817 + (-0.7213687 * `miR-9`) + (1.0082008 * `miR-34`) + (0.1832857 * `miR-132`) + (0.7240845 * `miR-137`) + (-0.9344316 * `miR-941`)
  )

concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))

df_hc_corrected$risk <- 1 / (1 + exp(-df_hc_corrected$PI_HC))

predicted_labels <- ifelse(df_hc_corrected$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_hc_corrected$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

results_NAPLS <- data.frame(
  C = paste0(
    round(concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))$concordance, 3), " (",
    round(concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))$concordance - 1.96 * (concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))$cvar), 3), "-",
    round(concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))$concordance + 1.96 * (concordance(glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit")))$cvar), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(cm$byClass[11] * 100, 1), "%"
  ),
  sensitivity = paste0(
    round(cm$byClass[1] * 100, 1), "%"
  ),
  specificity = paste0(
    round(cm$byClass[2] * 100, 1), "%"
  ),
  ppv = paste0(
    round(cm$byClass[3] * 100, 1), "%"
  ),
  npv = paste0(
    round(cm$byClass[4] * 100, 1), "%"
  )
)

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_hc_corrected$risk
)

calibration <- calibration %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)

cal_plot_breaks(calibration, truth = observed, estimate = predicted)

recal_model <- glm(Transition ~ PI_HC, data = df_hc_corrected, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = predicted_labels_recal
)

calibration_recal <- calibration_recal %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  )
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted, cal_plot = FALSE)

##### Likelihood Ratio Plot #####
classification_chr <- data.frame(
  threshold = numeric(),
  balanced_accuracy = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  ppv = numeric(),
  npv = numeric(),
  LR_pos = numeric(),
  LR_neg = numeric()
)
recal_model <- glm(Transition ~ PI_CHR, data = df_chr, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")

for (i in seq(0.01, 1, by = 0.01)) {
  predicted_labels <- ifelse(recalibrated_probs >= i, 1, 0)
  observed <- as.numeric(as.factor(df_chr$Transition)) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

  classification_chr <- rbind(classification_chr, data.frame(
    threshold = i,
    balanced_accuracy = cm$byClass[11],
    sensitivity = cm$byClass[1],
    specificity = cm$byClass[2],
    ppv = cm$byClass[3],
    npv = cm$byClass[4],
    LR_pos = posLr(obs = observed, pred = predicted_labels, pos_level = 1)$posLr,
    LR_neg = negLr(obs = observed, pred = predicted_labels, pos_level = 1)$negLr
  ))
}

ggplot(data = classification_chr, aes(x = threshold * 100)) +
  geom_line(aes(y = LR_pos, color = "Positive Likelihood Ratio"), size = 1.5) +
  geom_line(aes(y = LR_neg, color = "Negative Likelihood Ratio"), size = 1.5) +
  labs(
    x = "Threshold Probability (%)",
    y = "Likelihood Ratio",
    color = "Metric"
  ) +
  scale_color_manual(
    values = c("Negative Likelihood Ratio" = "#c8526a", "Positive Likelihood Ratio" = "#599ec4"),
    labels = c("Negative Likelihood Ratio", "Positive Likelihood Ratio")
  ) +
  theme_classic() +
  theme(text = element_text(family = "Roboto", face = "bold", size = 20))
ggsave("likelihood_ratio_plot_210725.png", width = 42, height = 32, units = "cm", scale = 0.65)

ggplot(data = classification_chr, aes(x = threshold * 100)) +
  geom_line(aes(y = ppv, color = "PPV"), size = 1.5) +
  geom_line(aes(y = npv, color = "NPV"), size = 1.5) +
  labs(
    x = "Threshold Probability (%)",
    y = "Positive/Negative Predictive Value",
    color = "Metric"
  ) +
  scale_color_manual(
    values = c("NPV" = "#c8526a", "PPV" = "#599ec4"),
    labels = c("NPV", "PPV")
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 20))
ggsave("ppv_npv_plot_210725.png", width = 42, height = 32, units = "cm", scale = 0.65)

ggplot(data = classification_chr, aes(x = threshold * 100)) +
  geom_line(aes(y = specificity, color = "Specificity"), size = 1.5) +
  geom_line(aes(y = sensitivity, color = "Sensitivity"), size = 1.5) +
  labs(
    x = "Threshold Probability (%)",
    y = "Sensitivity/Specificity",
    color = "Metric"
  ) +
  scale_color_manual(
    values = c("Sensitivity" = "#c8526a", "Specificity" = "#599ec4"),
    labels = c("Sensitivity", "Specificity")
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 20))
ggsave("sens_spec_plot_210725.png", width = 42, height = 32, units = "cm", scale = 0.65)
