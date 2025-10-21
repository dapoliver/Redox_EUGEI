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
library(tibble)
library(Metrics)
library(rms)
library(predtools)
library(Hmisc)
library(readxl)
library(survival)
library(survminer)
library(predRupdate)
library(metrica)

set.seed(123)
df <- read_excel("Eugei vf 0810 final BATCH 091025.xlsx")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

df_chr <- df %>% filter(Group != "Control")
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)

df_cc <- df_chr %>% subset(select = c(MIR132, MIR34A, MIR9, MIR941, MIR137, Transition_status, site))
df_cc <- df_cc[complete.cases(df_cc), ]

data <- df_cc

data$Transition_status <- as.factor((as.character(data$Transition_status)))
levels(data$Transition_status) <- c("NT", "T")
summary(as.factor(data$Transition_status))

# Define the predictor sets
predictors <- list(
  a = c("MIR132", "MIR34A", "MIR9", "MIR941", "MIR137")
)

# Create data frames to store predictions and coefficients
all_predictions <- data.frame(
  Subject_ID = integer(), True_Label = integer(),
  Predicted_Probability = numeric(),
  Fold = integer(), Repeat = integer(), stringsAsFactors = FALSE
)

coef_results <- data.frame(Variable = character(), Coefficient = numeric(), stringsAsFactors = FALSE)

# Create 5-fold cross-validation with 5 repeats for outer folds
outer_folds <- createMultiFolds(data$Transition_status, k = 5, times = 10)

temp_results <- data.frame(
  fold_num = numeric(),
  C = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  balanced_accuracy = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  LR_pos = numeric(),
  LR_neg = numeric(),
  stringsAsFactors = FALSE
)
# Loop over each outer fold
for (fold_num in seq_along(outer_folds)) {
  # set.seed(123 + fold_num) # Ensure reproducibility for each fold
  train_idx <- outer_folds[[fold_num]]
  test_idx <- setdiff(seq_len(nrow(data)), train_idx)

  train <- data[train_idx, ]
  test <- data[test_idx, ]

  ### Compute Global Means ###
  global_mean <- colMeans(test[, 1:5, drop = FALSE], na.rm = TRUE)

  ### Mean Offset Correction ###
  batch_test <- as.factor(test$site)
  test_corrected <- test # Start with the original data

  for (b in levels(batch_test)) {
    batch_indices <- which(batch_test == b) # Indices for samples in batch `b`

    if (length(batch_indices) == 0) {
      warning(paste("Batch", b, "is empty. Skipping."))
      next
    }

    # Extract the subset of test for the current batch
    batch_data <- test[batch_indices, 1:5, drop = FALSE] # Exclude site and outcome

    # Compute means for each predictor in this batch
    batch_mean <- colMeans(batch_data, na.rm = TRUE)

    if (length(batch_mean) == 0) {
      warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
      next
    }

    # Compute the offset: batch mean - global mean
    offset <- batch_mean - global_mean

    # Subtract the offset to align the batch with the global mean
    test_corrected[batch_indices, 1:5] <- sweep(
      batch_data,
      1,
      offset,
      "-"
    )
  }
  test <- test_corrected # %>% subset(select=c(-site))

  ### Mean offset correction ###
  batch_train <- as.factor(train$site)
  train_corrected <- train # Start with the original data

  ### Compute Global Means ###
  global_mean <- colMeans(train[, 1:5, drop = FALSE], na.rm = TRUE)
  # Iterate over each batch
  for (b in levels(batch_train)) {
    batch_indices <- which(batch_train == b) # Indices for samples in batch `b`

    if (length(batch_indices) == 0) {
      warning(paste("Batch", b, "is empty. Skipping."))
      next
    }

    # Extract the subset of test for the current batch
    batch_data <- train[batch_indices, 1:5, drop = FALSE]

    # Compute row-wise means for this batch
    batch_mean <- colMeans(batch_data, na.rm = TRUE)

    # Compute the offset: batch mean - global mean
    offset <- batch_mean - global_mean

    if (length(batch_mean) == 0) {
      warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
      next
    }

    # Subtract batch means
    train_corrected[batch_indices, 1:5] <- sweep(
      batch_data,
      1,
      offset,
      "-"
    )
    cat("Processed batch", b, "\n")
  }

  train <- train_corrected # %>% subset(select=c(-site))
  x <- model.matrix(~ . - 1, train[, predictors[[1]]])
  y <- train$Transition_status

  # Define fold IDs for cross-validation
  # set.seed(123 + fold_num) # Ensure reproducibility for fold IDs
  foldid <- sample(rep(1:5, length.out = length(y)))

  # Fit model and predict
  cv_model <- cv.glmnet(x, y,
    family = "binomial", alpha = 1,
    nfolds = 5, foldid = foldid
  )

  final_model <- glmnet(x, y,
    family = "binomial", alpha = 1,
    lambda = cv_model$lambda.min
  )

  x <- model.matrix(~ . - 1, data[, predictors[[1]]])
  y <- data$Transition_status

  predictions <- as.vector(predict(final_model, newx = x[test_idx, ], type = "response"))

  # Create confusion matrix for a threshold (e.g., 0.5)
  predicted_labels <- ifelse(predictions >= 0.5, 1, 0)
  observed <- as.numeric(as.factor(y[test_idx])) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")
  model_test <- glm(y[test_idx] ~ predictions, family = binomial)

  # Extract performance measures
  temp_results <- rbind(temp_results, data.frame(
    fold_num = fold_num,
    C = concordance(model_test)$concordance,
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    balanced_accuracy = cm$byClass["Balanced Accuracy"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    LR_pos = posLr(obs = observed, pred = predicted_labels, pos_level = 1)$posLr,
    LR_neg = negLr(obs = observed, pred = predicted_labels, pos_level = 1)$negLr
  ))
  # Save predictions
  all_predictions <- rbind(all_predictions, data.frame(
    Subject_ID = test_idx, True_Label = y[test_idx], Predicted_Probability = predictions,
    Fold = fold_num, Repeat = NA
  ))
}

temp_results <- temp_results %>% mutate(LR_pos = case_when(
  is.infinite(LR_pos) ~ NA,
  TRUE ~ LR_pos
))

temp_mean <- temp_results %>% aggregate(. ~ 1, FUN = "mean", na.rm = TRUE)
temp_se <- temp_results %>% aggregate(. ~ 1, FUN = function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))

# Fit the final model on the full dataset and save coefficients
foldid <- sample(rep(1:5, length.out = length(y)))

all <- rbind(train, test)
x <- model.matrix(~ . - 1, all[, predictors[[1]]])
y <- all$Transition_status

cv_model <- cv.glmnet(x, y,
  family = "binomial", alpha = 1,
  nfolds = 5, foldid = foldid
)

final_full_model_chr <- glmnet(x, y,
  family = "binomial", alpha = 1,
  lambda = cv_model$lambda.min
)

predictions_data <- all_predictions %>% mutate(
  True_Label = case_when(
    True_Label == "T" ~ 1,
    TRUE ~ 0
  ),
  Predicted_Probability = case_when(
    Predicted_Probability > 0.999 ~ 0.999,
    Predicted_Probability < 0.001 ~ 0.001,
    TRUE ~ Predicted_Probability
  )
) # Recode True_Label to binary format

# Calibration analysis
calibration <- data.frame(
  observed = predictions_data$True_Label, # Already in binary format
  predicted = predictions_data$Predicted_Probability
)

# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)
cal_plot_breaks(calibration, truth = observed, estimate = predicted)

results_new <- data.frame(
  C = paste0(
    round(temp_mean$C, 3), " (",
    round(temp_mean$C - 1.96 * (temp_se$C), 3), "-",
    round(temp_mean$C + 1.96 * (temp_se$C), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(temp_mean$balanced_accuracy * 100, 1), "% (",
    round(temp_mean$balanced_accuracy * 100 - 1.96 * (temp_se$balanced_accuracy) * 100, 1), "%-",
    round(temp_mean$balanced_accuracy * 100 + 1.96 * (temp_se$balanced_accuracy) * 100, 1), "%)"
  ),
  sensitivity = paste0(
    round(temp_mean$sensitivity * 100, 1), "% (",
    round(temp_mean$sensitivity * 100 - 1.96 * (temp_se$sensitivity) * 100, 1), "%-",
    round(temp_mean$sensitivity * 100 + 1.96 * (temp_se$sensitivity) * 100, 1), "%)"
  ),
  specificity = paste0(
    round(temp_mean$specificity * 100, 1), "% (",
    round(temp_mean$specificity * 100 - 1.96 * (temp_se$specificity) * 100, 1), "%-",
    round(temp_mean$specificity * 100 + 1.96 * (temp_se$specificity) * 100, 1), "%)"
  ),
  ppv = paste0(
    round(temp_mean$PPV * 100, 1), "% (",
    round(temp_mean$PPV * 100 - 1.96 * (temp_se$PPV) * 100, 1), "%-",
    round(temp_mean$PPV * 100 + 1.96 * (temp_se$PPV) * 100, 1), "%)"
  ),
  npv = paste0(
    round(temp_mean$NPV * 100, 1), "% (",
    round(temp_mean$NPV * 100 - 1.96 * (temp_se$NPV) * 100, 1), "%-",
    round(temp_mean$NPV * 100 + 1.96 * (temp_se$NPV) * 100, 1), "%)"
  ),
  posLR = paste0(
    round(temp_mean$LR_pos, 3), " (",
    round(temp_mean$LR_pos - 1.96 * (temp_se$LR_pos), 3), "-",
    round(temp_mean$LR_pos + 1.96 * (temp_se$LR_pos), 3), ")"
  ),
  negLR = paste0(
    round(temp_mean$LR_neg, 3), " (",
    round(temp_mean$LR_neg - 1.96 * (temp_se$LR_neg), 3), "-",
    round(temp_mean$LR_neg + 1.96 * (temp_se$LR_neg), 3), ")"
  ),
  brier = paste0(round(logistic_calibration$BrierScore[1], 2), "(", round(logistic_calibration$Brier_lower[1], 2), "-", round(logistic_calibration$Brier_upper[1], 2), ")"),
  calibration_intercept = paste0(round(logistic_calibration$CalInt[1], 2), " (", round(logistic_calibration$CalInt_lower[1], 2), "-", round(logistic_calibration$CalInt_upper[1], 2), ")"),
  calibration_slope = paste0(round(logistic_calibration$CalSlope[1], 2), " (", round(logistic_calibration$CalSlope_lower[1], 2), "-", round(logistic_calibration$CalSlope_upper[1], 2), ")")
)

write_csv(results_new, "CV_results_mean_offset_171025.csv")

##### External Validation #####

df_NAPLS <- read_excel("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/NAPLS/NAPLS BATCH info17102025.xlsx")

summary(factor(df_NAPLS$`GROUP (UC = CTRL group)`))
df_NAPLS <- df_NAPLS %>% subset(select = c(
  demo_age_ym, demo_sex, demo_racial, `miR-9`, `miR-34`, `miR-132`, `miR-137`, `miR-941`, `GROUP (UC = CTRL group)`, SiteNumber,
  P1_SOPS, P2_SOPS, P3_SOPS, P4_SOPS, P5_SOPS, GlobalAssessmentFunction
))

df_NAPLS <- df_NAPLS %>%
  dplyr::rename(MIR9 = `miR-9`, MIR34A = `miR-34`, MIR132 = `miR-132`, MIR137 = `miR-137`, MIR941 = `miR-941`, site = SiteNumber)
df_NAPLS <- df_NAPLS %>%
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
    PI_CHR = predict(final_full_model_chr, newx = as.matrix(df_NAPLS[, c("MIR9", "MIR34A", "MIR132", "MIR137", "MIR941")]), type = "link")[, 1],
    Transition = ifelse(`GROUP (UC = CTRL group)` == "CHR-C", 1, 0),
    ethnicity = case_when(
      demo_racial == "European" ~ "White",
      demo_racial == "African" ~ "Black",
      demo_racial == "East Asian" | demo_racial == "South Asian" ~ "Asian",
      demo_racial == "Interracial" ~ "Mixed",
      TRUE ~ "Other"
    )
  )
df_NAPLS <- df_NAPLS[complete.cases(df_NAPLS), ]

df_NAPLS_chr <- df_NAPLS %>% filter(`GROUP (UC = CTRL group)` != "UC")

### Compute Global Means ###
global_mean <- colMeans(df_NAPLS_chr[, 4:8, drop = FALSE], na.rm = TRUE)

### Mean Offset Correction ###
batch_test <- as.factor(df_NAPLS_chr$site)
test_corrected <- df_NAPLS_chr # Start with the original data

for (b in levels(batch_test)) {
  batch_indices <- which(batch_test == b) # Indices for samples in batch `b`

  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }

  # Extract the subset of test for the current batch
  batch_data <- df_NAPLS_chr[batch_indices, 4:8, drop = FALSE] # Exclude site and outcome

  # Compute means for each predictor in this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }

  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  # Subtract the offset to align the batch with the global mean
  test_corrected[batch_indices, 4:8] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
}
df_NAPLS_chr <- test_corrected # %>% subset(select=c(-site))

concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))

df_NAPLS_chr$risk <- 1 / (1 + exp(-df_NAPLS_chr$PI_CHR))

predicted_labels <- ifelse(df_NAPLS_chr$risk >= 0.5, 1, 0)
observed <- as.numeric(as.factor(df_NAPLS_chr$Transition)) - 1
cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

calibration <- data.frame(
  observed = observed, # Already in binary format
  predicted = df_NAPLS_chr$risk
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
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted, cal_plot = FALSE)

cal_plot_breaks(calibration, truth = observed, estimate = predicted)

results_NAPLS <- data.frame(
  C = paste0(
    round(concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))$concordance, 3), " (",
    round(concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))$concordance - 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))$cvar), 3), "-",
    round(concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))$concordance + 1.96 * (concordance(glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit")))$cvar), 3), ")"
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
  ),
  LR_pos = posLr(obs = observed, pred = predicted_labels, pos_level = 1)$posLr,
  LR_neg = negLr(obs = observed, pred = predicted_labels, pos_level = 1)$negLr,
  calibration_intercept = paste0(round(logistic_calibration$CalInt[1], 2), " (", round(logistic_calibration$CalInt_lower[1], 2), "-", round(logistic_calibration$CalInt_upper[1], 2), ")"),
  calibration_slope = paste0(round(logistic_calibration$CalSlope[1], 2), " (", round(logistic_calibration$CalSlope_lower[1], 2), "-", round(logistic_calibration$CalSlope_upper[1], 2), ")"),
  brier = paste0(round(logistic_calibration$BrierScore[1], 2), "(", round(logistic_calibration$Brier_lower[1], 2), "-", round(logistic_calibration$Brier_upper[1], 2), ")")
)
write.csv(results_NAPLS, "external_validation_results_mean_offset_171025.csv", row.names = FALSE)

recal_model <- glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit"))
recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(recalibrated_probs >= 0.5, 1, 0)
cm_recal <- confusionMatrix(as.factor(predicted_labels_recal), as.factor(observed), positive = "1")

# Recalibrate
calibration_recal <- data.frame(
  observed = observed, # Already in binary format
  predicted = recalibrated_probs
)

calibration_recal <- calibration_recal %>%
  rowwise() %>%
  mutate(
    predicted = case_when(
      predicted > 0.9999 ~ 0.9999,
      predicted < 0.0001 ~ 0.0001,
      TRUE ~ predicted
    )
  ) %>%
  ungroup()
logistic_calibration_recal <- predRupdate::pred_val_probs(binary_outcome = calibration_recal$observed, Prob = calibration_recal$predicted, cal_plot = FALSE)

cal_plot_breaks(calibration_recal, truth = observed, estimate = predicted)

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
  ),
  LR_pos = posLr(obs = observed, pred = predicted_labels_recal, pos_level = 1)$posLr,
  LR_neg = negLr(obs = observed, pred = predicted_labels_recal, pos_level = 1)$negLr,
  brier = paste0(round(logistic_calibration_recal$BrierScore[1], 2), "(", round(logistic_calibration_recal$Brier_lower[1], 2), "-", round(logistic_calibration_recal$Brier_upper[1], 2), ")"),
  calibration_intercept = paste0(round(logistic_calibration_recal$CalInt[1], 2), " (", round(logistic_calibration_recal$CalInt_lower[1], 2), "-", round(logistic_calibration_recal$CalInt_upper[1], 2), ")"),
  calibration_slope = paste0(round(logistic_calibration_recal$CalSlope[1], 2), " (", round(logistic_calibration_recal$CalSlope_lower[1], 2), "-", round(logistic_calibration_recal$CalSlope_upper[1], 2), ")")
)
write.csv(results_NAPLS_recal, "external_validation_results_recal_mean_offset_171025.csv", row.names = FALSE)
