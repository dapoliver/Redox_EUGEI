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
library(neuroCombat)

set.seed(123)
df <- read_excel("Eugei vf 0810 final BATCH 091025.xlsx")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

df_chr <- df %>% filter(Group != "AtRisk_NoTr")
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)

df_cc <- df_chr %>% subset(select = c(Group, Age, Gender.x, Ethnicity, MIR132, MIR34A, MIR9, MIR941, MIR137, site, `BATCH NUMB`, Transition_status))
df_cc <- df_cc %>% rename(batch = `BATCH NUMB`, Gender = Gender.x)
df_cc <- df_cc[complete.cases(df_cc), ]

data <- df_cc

data$Transition_status <- as.factor((as.character(data$Transition_status)))
levels(data$Transition_status) <- c("NT", "T")
summary(as.factor(data$Transition_status))

# Define the predictor sets
predictors <- list(
  a = c("MIR9", "MIR34A", "MIR132", "MIR137", "MIR941")
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

  combat_train <- neuroCombat(dat = t(train[, predictors[[1]]]), batch = train$batch, mod = NULL)
  train[, predictors[[1]]] <- t(combat_train$dat.combat)

  est <- combat_train$estimates

  gamma_star <- est$gamma.star
  delta_star <- est$delta.star
  stand_mean <- est$stand.mean
  var_pooled <- est$var.pooled
  test_mat <- t(test[, predictors[[1]]])
  batch_test <- test$batch

  # Standardize test data using TRAIN parameters
  s_data <- (test_mat - stand_mean[, 1]) / sqrt(var_pooled)

  # Apply batch adjustment
  for (b in levels(batch_test)) {
    batch_idx <- which(batch_test == b)

    s_data[, batch_idx] <-
      (s_data[, batch_idx] - gamma_star[, b]) /
        sqrt(delta_star[, b])
  }

  # Back-transform
  test_corrected <- t(s_data)

  test[, predictors[[1]]] <- test_corrected * sqrt(var_pooled) + stand_mean[, 1]

  ### Compute Global Means ###
  global_mean <- colMeans(test[, predictors[[1]], drop = FALSE], na.rm = TRUE)

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
    batch_data <- test[batch_indices, predictors[[1]], drop = FALSE] # Exclude site and outcome

    # Compute means for each predictor in this batch
    batch_mean <- colMeans(batch_data, na.rm = TRUE)

    if (length(batch_mean) == 0) {
      warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
      next
    }

    # Compute the offset: batch mean - global mean
    offset <- batch_mean - global_mean

    # Subtract the offset to align the batch with the global mean
    test_corrected[batch_indices, predictors[[1]]] <- sweep(
      batch_data,
      1,
      offset,
      "-"
    )
  }
  test <- test_corrected

  ### Mean offset correction ###
  batch_train <- as.factor(train$site)
  train_corrected <- train # Start with the original data

  ### Compute Global Means ###
  global_mean <- colMeans(train[, predictors[[1]], drop = FALSE], na.rm = TRUE)
  # Iterate over each batch
  for (b in levels(batch_train)) {
    batch_indices <- which(batch_train == b) # Indices for samples in batch `b`

    if (length(batch_indices) == 0) {
      warning(paste("Batch", b, "is empty. Skipping."))
      next
    }

    # Extract the subset of test for the current batch
    batch_data <- train[batch_indices, predictors[[1]], drop = FALSE]

    # Compute row-wise means for this batch
    batch_mean <- colMeans(batch_data, na.rm = TRUE)

    # Compute the offset: batch mean - global mean
    offset <- batch_mean - global_mean

    if (length(batch_mean) == 0) {
      warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
      next
    }

    # Subtract batch means
    train_corrected[batch_indices, predictors[[1]]] <- sweep(
      batch_data,
      1,
      offset,
      "-"
    )
    cat("Processed batch", b, "\n")
  }

  models <- list()

  # Back-transform
  test_corrected <- test

  test_corrected[, predictors[[1]]] <- t(s_data)

  test_corrected[, predictors[[1]]] <- test_corrected[, predictors[[1]]] * sqrt(var_pooled) + stand_mean[, 1]

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

  x_test <- as.matrix(test_corrected[, predictors[[1]]])
  predictions <- as.vector(predict(final_model, newx = x_test, type = "response"))

  # Create confusion matrix for a threshold (e.g., 0.5)
  predicted_labels <- ifelse(predictions >= 0.5, 1, 0)
  observed <- as.numeric(as.factor(test_corrected$Transition_status)) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")
  model_test <- glm(test_corrected$Transition_status ~ predictions, family = binomial)

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
    Subject_ID = test_idx, True_Label = test_corrected$Transition_status, Predicted_Probability = predictions,
    Fold = fold_num, Repeat = NA
  ))
}

data_corrected <- data # Start with the original data

combat <- neuroCombat(dat = t(data[, predictors[[1]]]), batch = data$batch, mod = NULL)
data_corrected[, predictors[[1]]] <- t(combat$dat.combat)

### Mean offset correction ###
batch <- as.factor(data$site)

### Compute Global Means ###
global_mean <- colMeans(data[, predictors[[1]], drop = FALSE], na.rm = TRUE)
# Iterate over each batch
for (b in levels(batch)) {
  batch_indices <- which(batch == b) # Indices for samples in batch `b`

  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }

  # Extract the subset of test for the current batch
  batch_data <- data[batch_indices, predictors[[1]], drop = FALSE]

  # Compute row-wise means for this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }

  # Subtract batch means
  data_corrected[batch_indices, predictors[[1]]] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
  cat("Processed batch", b, "\n")
}

x <- model.matrix(~ . - 1, data_corrected[, predictors[[1]]])
y <- data_corrected$Transition_status
foldid <- sample(rep(1:5, length.out = length(y)))

cv_model <- cv.glmnet(x, y,
  family = "binomial", alpha = 1,
  nfolds = 5, foldid = foldid
)

final_full_model_chr <- glmnet(x, y,
  family = "binomial", alpha = 1,
  lambda = cv_model$lambda.min
)

# Extract coefficients and store
coeffs <- data.frame(
  Variable = rownames(coef(final_full_model_chr)),
  Coefficient = as.vector(coef(final_full_model_chr))
)
coef_results <- rbind(coef_results, coeffs)

write.csv(coef_results, "Results/logistic_regression_LASSO_coefficients_110325_HC.csv", row.names = FALSE)

# ============================================================
# Bootstrap confidence intervals using subject-level resampling
# ============================================================

library(dplyr)
library(pROC)
library(caret)
library(writexl)

# collapse repeated predictions per subject
boot_df <- all_predictions %>%
  group_by(Subject_ID) %>%
  summarise(
    obs = first(factor(as.numeric(as.factor(True_Label)) - 1, levels = c(0, 1))),
    pred = mean(Predicted_Probability),
    .groups = "drop"
  )

B <- 2000
n <- nrow(boot_df)

# storage
C_boot <- numeric(B)
BA_boot <- numeric(B)
SEN_boot <- numeric(B)
SPEC_boot <- numeric(B)
PPV_boot <- numeric(B)
NPV_boot <- numeric(B)
LR_pos_boot <- numeric(B)
LR_neg_boot <- numeric(B)

Intercept_boot <- numeric(B)
Slope_boot <- numeric(B)
Brier_boot <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  d <- boot_df[idx, ]

  # ---------- C-index ----------
  C_boot[b] <- concordance(glm(obs ~ pred, data = d, family = "binomial"))$concordance

  # ---------- Classification metrics ----------
  pred_class <- ifelse(d$pred >= 0.5, 1, 0)

  cm <- caret::confusionMatrix(
    factor(pred_class, levels = c(0, 1)),
    d$obs,
    positive = "1"
  )

  BA_boot[b] <- cm$byClass["Balanced Accuracy"]
  SEN_boot[b] <- cm$byClass["Sensitivity"]
  SPEC_boot[b] <- cm$byClass["Specificity"]
  PPV_boot[b] <- cm$byClass["Pos Pred Value"]
  NPV_boot[b] <- cm$byClass["Neg Pred Value"]
  LR_pos_boot[b] <- posLr(obs = d$obs, pred = pred_class, pos_level = 1)$posLr
  LR_neg_boot[b] <- negLr(obs = d$obs, pred = pred_class, pos_level = 1)$negLr

  # ---------- Calibration ----------

  logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = as.numeric(d$obs) - 1, Prob = d$pred)

  Intercept_boot[b] <- logistic_calibration$CalInt[1]
  Slope_boot[b] <- logistic_calibration$CalSlope[1]

  Brier_boot[b] <- logistic_calibration$BrierScore[1]
}

# helper function
boot_summary <- function(x) {
  q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
  c(
    Mean = mean(x, na.rm = TRUE),
    CI_low = unname(q[1]),
    CI_high = unname(q[2])
  )
}

# classification metrics
class_boot_summary <- rbind(
  C = boot_summary(C_boot),
  BalancedAccuracy = boot_summary(BA_boot),
  Sensitivity = boot_summary(SEN_boot),
  Specificity = boot_summary(SPEC_boot),
  PPV = boot_summary(PPV_boot),
  NPV = boot_summary(NPV_boot),
  LR_pos = boot_summary(LR_pos_boot),
  LR_neg = boot_summary(LR_neg_boot)
)

# calibration metrics
calib_boot_summary <- rbind(
  Intercept = boot_summary(Intercept_boot),
  Slope = boot_summary(Slope_boot),
  Brier = boot_summary(Brier_boot)
)

class_boot_summary <- as.data.frame(class_boot_summary)
colnames(class_boot_summary) <- c("Mean", "CI_low", "CI_high")
class_boot_summary$Metric <- rownames(class_boot_summary)
rownames(class_boot_summary) <- NULL

calib_boot_summary <- as.data.frame(calib_boot_summary)
colnames(calib_boot_summary) <- c("Mean", "CI_low", "CI_high")
calib_boot_summary$Metric <- rownames(calib_boot_summary)
rownames(calib_boot_summary) <- NULL

# reorder columns
class_boot_summary <- class_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

calib_boot_summary <- calib_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

# helper to get a metric row from summary df
get_metric <- function(df, metric_name) {
  row <- df[df$Metric == metric_name, ]
  if (nrow(row) == 0) {
    return(NULL)
  }
  list(mean = row$Mean, low = row$CI_low, high = row$CI_high)
}

# helper formatters
fmt_num <- function(mean, low, high, digits = 3) {
  paste0(
    formatC(mean, format = "f", digits = digits), " (",
    formatC(low, format = "f", digits = digits), "-",
    formatC(high, format = "f", digits = digits), ")"
  )
}

fmt_pct <- function(mean, low, high, digits_mean = 1, digits_ci = 1) {
  paste0(
    formatC(mean * 100, format = "f", digits = digits_mean), "% (",
    formatC(low * 100, format = "f", digits = digits_ci), "%-",
    formatC(high * 100, format = "f", digits = digits_ci), "%)"
  )
}

# pull classification metrics (expects metrics named exactly as in your class_boot_summary)
C_m <- get_metric(class_boot_summary, "C")
BA_m <- get_metric(class_boot_summary, "BalancedAccuracy")
SEN_m <- get_metric(class_boot_summary, "Sensitivity")
SPEC_m <- get_metric(class_boot_summary, "Specificity")
PPV_m <- get_metric(class_boot_summary, "PPV")
NPV_m <- get_metric(class_boot_summary, "NPV")

# optionally LR metrics if present in your class_boot_summary
LRpos_m <- get_metric(class_boot_summary, "LR_pos") # or "posLR" depending on your Metric name
LRneg_m <- get_metric(class_boot_summary, "LR_neg") # or "negLR"

# calibration / brier (expects names in calib_boot_summary)
Intercept_m <- get_metric(calib_boot_summary, "Intercept")
Slope_m <- get_metric(calib_boot_summary, "Slope")
Brier_m <- get_metric(calib_boot_summary, "Brier")

# Build the results_new row, using fallbacks if any metric is missing
results_new <- data.frame(
  C = if (!is.null(C_m)) fmt_num(C_m$mean, C_m$low, C_m$high, digits = 3) else NA_character_,
  balanced_accuracy = if (!is.null(BA_m)) fmt_pct(BA_m$mean, BA_m$low, BA_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  sensitivity = if (!is.null(SEN_m)) fmt_pct(SEN_m$mean, SEN_m$low, SEN_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  specificity = if (!is.null(SPEC_m)) fmt_pct(SPEC_m$mean, SPEC_m$low, SPEC_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  ppv = if (!is.null(PPV_m)) fmt_pct(PPV_m$mean, PPV_m$low, PPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  npv = if (!is.null(NPV_m)) fmt_pct(NPV_m$mean, NPV_m$low, NPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  posLR = if (!is.null(LRpos_m)) fmt_num(LRpos_m$mean, LRpos_m$low, LRpos_m$high, digits = 3) else NA_character_,
  negLR = if (!is.null(LRneg_m)) fmt_num(LRneg_m$mean, LRneg_m$low, LRneg_m$high, digits = 3) else NA_character_,
  brier = if (!is.null(Brier_m)) paste0(formatC(Brier_m$mean, format = "f", digits = 2), " (", formatC(Brier_m$low, format = "f", digits = 2), "-", formatC(Brier_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_intercept = if (!is.null(Intercept_m)) paste0(formatC(Intercept_m$mean, format = "f", digits = 2), " (", formatC(Intercept_m$low, format = "f", digits = 2), "-", formatC(Intercept_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_slope = if (!is.null(Slope_m)) paste0(formatC(Slope_m$mean, format = "f", digits = 2), " (", formatC(Slope_m$low, format = "f", digits = 2), "-", formatC(Slope_m$high, format = "f", digits = 2), ")") else NA_character_,
  stringsAsFactors = FALSE
)

results_new
write_csv(results_new, "Results/CV_results_mean_offset_170326_HC.csv")

calibration <- data.frame(
  observed = as.numeric(all_predictions$True_Label) - 1, # Already in binary format
  predicted = all_predictions$Predicted_Probability
)

# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)
cal_plot_breaks(calibration, truth = observed, estimate = predicted)

##### External Validation #####

df_NAPLS <- read_excel("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/NAPLS/NAPLS BATCH info17102025.xlsx")

summary(factor(df_NAPLS$`GROUP (UC = CTRL group)`))
df_NAPLS <- df_NAPLS %>% subset(select = c(
  demo_age_ym, demo_sex, demo_racial, `miR-9`, `miR-34`, `miR-132`, `miR-137`, `miR-941`, `GROUP (UC = CTRL group)`, BATCH, SiteNumber,
  P1_SOPS, P2_SOPS, P3_SOPS, P4_SOPS, P5_SOPS, GlobalAssessmentFunction
))

df_NAPLS <- df_NAPLS %>%
  dplyr::rename(MIR9 = `miR-9`, MIR34A = `miR-34`, MIR132 = `miR-132`, MIR137 = `miR-137`, MIR941 = `miR-941`, site = SiteNumber, batch = BATCH)
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

df_NAPLS_chr <- df_NAPLS %>% filter(`GROUP (UC = CTRL group)` != "CHR-NC")

predictors <- list(
  a = c("MIR9", "MIR34A", "MIR132", "MIR137", "MIR941")
)

df_NAPLS_combat <- neuroCombat(dat = t(df_NAPLS_chr[, predictors[[1]]]), batch = df_NAPLS_chr$batch, mod = NULL)
df_NAPLS_chr[, predictors[[1]]] <- as.data.frame(t(df_NAPLS_combat$dat.combat))

### Compute Global Means ###
global_mean <- colMeans(df_NAPLS_chr[, predictors[[1]], drop = FALSE], na.rm = TRUE)

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
  batch_data <- df_NAPLS_chr[batch_indices, predictors[[1]], drop = FALSE] # Exclude site and outcome

  # Compute means for each predictor in this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)

  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }

  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean

  # Subtract the offset to align the batch with the global mean
  test_corrected[batch_indices, predictors[[1]]] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
}

df_NAPLS_chr <- test_corrected # %>% subset(select=c(-site))
df_NAPLS_chr$obs <- factor(as.numeric(as.factor(df_NAPLS_chr$Transition)) - 1, levels = c(0, 1))

df_NAPLS_chr$PI_CHR <- predict(final_full_model_chr, newx = as.matrix(df_NAPLS_chr[, predictors[[1]]]), type = "link")[, 1]
df_NAPLS_chr$pred <- predict(final_full_model_chr, newx = as.matrix(df_NAPLS_chr[, predictors[[1]]]), type = "response")[, 1]
n <- nrow(df_NAPLS_chr)

# storage
C_boot <- numeric(B)
BA_boot <- numeric(B)
SEN_boot <- numeric(B)
SPEC_boot <- numeric(B)
PPV_boot <- numeric(B)
NPV_boot <- numeric(B)
LR_pos_boot <- numeric(B)
LR_neg_boot <- numeric(B)

Intercept_boot <- numeric(B)
Slope_boot <- numeric(B)
Brier_boot <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  d <- df_NAPLS_chr[idx, ]

  # ---------- C-index ----------
  C_boot[b] <- concordance(glm(obs ~ pred, data = d, family = "binomial"))$concordance

  # ---------- Classification metrics ----------
  pred_class <- ifelse(d$pred >= 0.5, 1, 0)

  cm <- caret::confusionMatrix(
    factor(pred_class, levels = c(0, 1)),
    d$obs,
    positive = "1"
  )

  BA_boot[b] <- cm$byClass["Balanced Accuracy"]
  SEN_boot[b] <- cm$byClass["Sensitivity"]
  SPEC_boot[b] <- cm$byClass["Specificity"]
  PPV_boot[b] <- cm$byClass["Pos Pred Value"]
  NPV_boot[b] <- cm$byClass["Neg Pred Value"]
  LR_pos_boot[b] <- posLr(obs = d$obs, pred = pred_class, pos_level = 1)$posLr
  LR_neg_boot[b] <- negLr(obs = d$obs, pred = pred_class, pos_level = 1)$negLr

  # ---------- Calibration ----------

  logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = as.numeric(d$obs) - 1, Prob = d$pred)

  Intercept_boot[b] <- logistic_calibration$CalInt[1]
  Slope_boot[b] <- logistic_calibration$CalSlope[1]

  Brier_boot[b] <- logistic_calibration$BrierScore[1]
}

# classification metrics
class_boot_summary <- rbind(
  C = boot_summary(C_boot),
  BalancedAccuracy = boot_summary(BA_boot),
  Sensitivity = boot_summary(SEN_boot),
  Specificity = boot_summary(SPEC_boot),
  PPV = boot_summary(PPV_boot),
  NPV = boot_summary(NPV_boot),
  LR_pos = boot_summary(LR_pos_boot),
  LR_neg = boot_summary(LR_neg_boot)
)

# calibration metrics
calib_boot_summary <- rbind(
  Intercept = boot_summary(Intercept_boot),
  Slope = boot_summary(Slope_boot),
  Brier = boot_summary(Brier_boot)
)

class_boot_summary <- as.data.frame(class_boot_summary)
colnames(class_boot_summary) <- c("Mean", "CI_low", "CI_high")
class_boot_summary$Metric <- rownames(class_boot_summary)
rownames(class_boot_summary) <- NULL

calib_boot_summary <- as.data.frame(calib_boot_summary)
colnames(calib_boot_summary) <- c("Mean", "CI_low", "CI_high")
calib_boot_summary$Metric <- rownames(calib_boot_summary)
rownames(calib_boot_summary) <- NULL

# reorder columns
class_boot_summary <- class_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

calib_boot_summary <- calib_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

# pull classification metrics (expects metrics named exactly as in your class_boot_summary)
C_m <- get_metric(class_boot_summary, "C")
BA_m <- get_metric(class_boot_summary, "BalancedAccuracy")
SEN_m <- get_metric(class_boot_summary, "Sensitivity")
SPEC_m <- get_metric(class_boot_summary, "Specificity")
PPV_m <- get_metric(class_boot_summary, "PPV")
NPV_m <- get_metric(class_boot_summary, "NPV")

# optionally LR metrics if present in your class_boot_summary
LRpos_m <- get_metric(class_boot_summary, "LR_pos") # or "posLR" depending on your Metric name
LRneg_m <- get_metric(class_boot_summary, "LR_neg") # or "negLR"

# calibration / brier (expects names in calib_boot_summary)
Intercept_m <- get_metric(calib_boot_summary, "Intercept")
Slope_m <- get_metric(calib_boot_summary, "Slope")
Brier_m <- get_metric(calib_boot_summary, "Brier")

# Build the results_new row, using fallbacks if any metric is missing
results_NAPLS <- data.frame(
  C = if (!is.null(C_m)) fmt_num(C_m$mean, C_m$low, C_m$high, digits = 3) else NA_character_,
  balanced_accuracy = if (!is.null(BA_m)) fmt_pct(BA_m$mean, BA_m$low, BA_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  sensitivity = if (!is.null(SEN_m)) fmt_pct(SEN_m$mean, SEN_m$low, SEN_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  specificity = if (!is.null(SPEC_m)) fmt_pct(SPEC_m$mean, SPEC_m$low, SPEC_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  ppv = if (!is.null(PPV_m)) fmt_pct(PPV_m$mean, PPV_m$low, PPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  npv = if (!is.null(NPV_m)) fmt_pct(NPV_m$mean, NPV_m$low, NPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  posLR = if (!is.null(LRpos_m)) fmt_num(LRpos_m$mean, LRpos_m$low, LRpos_m$high, digits = 3) else NA_character_,
  negLR = if (!is.null(LRneg_m)) fmt_num(LRneg_m$mean, LRneg_m$low, LRneg_m$high, digits = 3) else NA_character_,
  brier = if (!is.null(Brier_m)) paste0(formatC(Brier_m$mean, format = "f", digits = 2), " (", formatC(Brier_m$low, format = "f", digits = 2), "-", formatC(Brier_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_intercept = if (!is.null(Intercept_m)) paste0(formatC(Intercept_m$mean, format = "f", digits = 2), " (", formatC(Intercept_m$low, format = "f", digits = 2), "-", formatC(Intercept_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_slope = if (!is.null(Slope_m)) paste0(formatC(Slope_m$mean, format = "f", digits = 2), " (", formatC(Slope_m$low, format = "f", digits = 2), "-", formatC(Slope_m$high, format = "f", digits = 2), ")") else NA_character_,
  stringsAsFactors = FALSE
)

results_NAPLS
write.csv(results_NAPLS, "Results/external_validation_results_mean_offset_170326_HC.csv", row.names = FALSE)

logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = df_NAPLS_chr$obs, Prob = df_NAPLS_chr$pred)
cal_plot_breaks(df_NAPLS_chr, truth = Transition, estimate = pred)

recal_model <- glm(Transition ~ PI_CHR, data = df_NAPLS_chr, family = binomial(link = "logit"))
df_NAPLS_chr$recalibrated_probs <- predict(recal_model, type = "response")
predicted_labels_recal <- ifelse(df_NAPLS_chr$recalibrated_probs >= 0.5, 1, 0)

# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = df_NAPLS_chr$obs, Prob = df_NAPLS_chr$recalibrated_probs)
cal_plot_breaks(df_NAPLS_chr, truth = Transition, estimate = recalibrated_probs)

# storage
C_boot <- numeric(B)
BA_boot <- numeric(B)
SEN_boot <- numeric(B)
SPEC_boot <- numeric(B)
PPV_boot <- numeric(B)
NPV_boot <- numeric(B)
LR_pos_boot <- numeric(B)
LR_neg_boot <- numeric(B)

Intercept_boot <- numeric(B)
Slope_boot <- numeric(B)
Brier_boot <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  d <- df_NAPLS_chr[idx, ]

  # ---------- C-index ----------
  C_boot[b] <- concordance(glm(obs ~ recalibrated_probs, data = d, family = "binomial"))$concordance

  # ---------- Classification metrics ----------
  pred_class <- ifelse(d$recalibrated_probs >= 0.5, 1, 0)

  cm <- caret::confusionMatrix(
    factor(pred_class, levels = c(0, 1)),
    d$obs,
    positive = "1"
  )

  BA_boot[b] <- cm$byClass["Balanced Accuracy"]
  SEN_boot[b] <- cm$byClass["Sensitivity"]
  SPEC_boot[b] <- cm$byClass["Specificity"]
  PPV_boot[b] <- cm$byClass["Pos Pred Value"]
  NPV_boot[b] <- cm$byClass["Neg Pred Value"]
  LR_pos_boot[b] <- posLr(obs = d$obs, pred = pred_class, pos_level = 1)$posLr
  LR_neg_boot[b] <- negLr(obs = d$obs, pred = pred_class, pos_level = 1)$negLr

  # ---------- Calibration ----------

  logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = as.numeric(d$obs) - 1, Prob = d$recalibrated_probs)

  Intercept_boot[b] <- logistic_calibration$CalInt[1]
  Slope_boot[b] <- logistic_calibration$CalSlope[1]

  Brier_boot[b] <- logistic_calibration$BrierScore[1]
}

# classification metrics
class_boot_summary <- rbind(
  C = boot_summary(C_boot),
  BalancedAccuracy = boot_summary(BA_boot),
  Sensitivity = boot_summary(SEN_boot),
  Specificity = boot_summary(SPEC_boot),
  PPV = boot_summary(PPV_boot),
  NPV = boot_summary(NPV_boot),
  LR_pos = boot_summary(LR_pos_boot),
  LR_neg = boot_summary(LR_neg_boot)
)

# calibration metrics
calib_boot_summary <- rbind(
  Intercept = boot_summary(Intercept_boot),
  Slope = boot_summary(Slope_boot),
  Brier = boot_summary(Brier_boot)
)

class_boot_summary <- as.data.frame(class_boot_summary)
colnames(class_boot_summary) <- c("Mean", "CI_low", "CI_high")
class_boot_summary$Metric <- rownames(class_boot_summary)
rownames(class_boot_summary) <- NULL

calib_boot_summary <- as.data.frame(calib_boot_summary)
colnames(calib_boot_summary) <- c("Mean", "CI_low", "CI_high")
calib_boot_summary$Metric <- rownames(calib_boot_summary)
rownames(calib_boot_summary) <- NULL

# reorder columns
class_boot_summary <- class_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

calib_boot_summary <- calib_boot_summary %>%
  select(Metric, Mean, CI_low, CI_high)

# pull classification metrics (expects metrics named exactly as in your class_boot_summary)
C_m <- get_metric(class_boot_summary, "C")
BA_m <- get_metric(class_boot_summary, "BalancedAccuracy")
SEN_m <- get_metric(class_boot_summary, "Sensitivity")
SPEC_m <- get_metric(class_boot_summary, "Specificity")
PPV_m <- get_metric(class_boot_summary, "PPV")
NPV_m <- get_metric(class_boot_summary, "NPV")

# optionally LR metrics if present in your class_boot_summary
LRpos_m <- get_metric(class_boot_summary, "LR_pos") # or "posLR" depending on your Metric name
LRneg_m <- get_metric(class_boot_summary, "LR_neg") # or "negLR"

# calibration / brier (expects names in calib_boot_summary)
Intercept_m <- get_metric(calib_boot_summary, "Intercept")
Slope_m <- get_metric(calib_boot_summary, "Slope")
Brier_m <- get_metric(calib_boot_summary, "Brier")

# Build the results_new row, using fallbacks if any metric is missing
results_NAPLS_recal <- data.frame(
  C = if (!is.null(C_m)) fmt_num(C_m$mean, C_m$low, C_m$high, digits = 3) else NA_character_,
  balanced_accuracy = if (!is.null(BA_m)) fmt_pct(BA_m$mean, BA_m$low, BA_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  sensitivity = if (!is.null(SEN_m)) fmt_pct(SEN_m$mean, SEN_m$low, SEN_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  specificity = if (!is.null(SPEC_m)) fmt_pct(SPEC_m$mean, SPEC_m$low, SPEC_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  ppv = if (!is.null(PPV_m)) fmt_pct(PPV_m$mean, PPV_m$low, PPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  npv = if (!is.null(NPV_m)) fmt_pct(NPV_m$mean, NPV_m$low, NPV_m$high, digits_mean = 1, digits_ci = 1) else NA_character_,
  posLR = if (!is.null(LRpos_m)) fmt_num(LRpos_m$mean, LRpos_m$low, LRpos_m$high, digits = 3) else NA_character_,
  negLR = if (!is.null(LRneg_m)) fmt_num(LRneg_m$mean, LRneg_m$low, LRneg_m$high, digits = 3) else NA_character_,
  brier = if (!is.null(Brier_m)) paste0(formatC(Brier_m$mean, format = "f", digits = 2), " (", formatC(Brier_m$low, format = "f", digits = 2), "-", formatC(Brier_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_intercept = if (!is.null(Intercept_m)) paste0(formatC(Intercept_m$mean, format = "f", digits = 2), " (", formatC(Intercept_m$low, format = "f", digits = 2), "-", formatC(Intercept_m$high, format = "f", digits = 2), ")") else NA_character_,
  calibration_slope = if (!is.null(Slope_m)) paste0(formatC(Slope_m$mean, format = "f", digits = 2), " (", formatC(Slope_m$low, format = "f", digits = 2), "-", formatC(Slope_m$high, format = "f", digits = 2), ")") else NA_character_,
  stringsAsFactors = FALSE
)
write.csv(results_NAPLS_recal, "Results/external_validation_results_recal_mean_offset_121225_HC.csv", row.names = FALSE)


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

for (i in seq(0.01, 0.99, by = 0.01)) {
  predicted_labels <- ifelse(df_NAPLS_chr$recalibrated_probs >= i, 1, 0)
  observed <- as.numeric(as.factor(df_NAPLS_chr$Transition)) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

  classification_chr <- rbind(classification_chr, data.frame(
    threshold = i,
    balanced_accuracy = cm$byClass[11],
    sensitivity = cm$byClass[1],
    specificity = cm$byClass[2],
    ppv = cm$byClass[3],
    npv = cm$byClass[4],
    LR_pos = NA,
    LR_neg = NA
  ))

  classification_chr$LR_pos[classification_chr$threshold == i] <- case_when(
    sum(predicted_labels) != 0 & sum(predicted_labels) != length(predicted_labels) ~ posLr(obs = observed, pred = predicted_labels, pos_level = 1)$posLr,
    TRUE ~ NA
  )
  classification_chr$LR_neg[classification_chr$threshold == i] <- case_when(
    sum(predicted_labels) != 0 & sum(predicted_labels) != length(predicted_labels) ~ negLr(obs = observed, pred = predicted_labels, pos_level = 1)$negLr,
    TRUE ~ NA
  )
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
ggsave("likelihood_ratio_plot_210326.png", width = 42, height = 32, units = "cm", scale = 0.65)

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
ggsave("ppv_npv_plot_210325.png", width = 42, height = 32, units = "cm", scale = 0.65)

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
ggsave("sens_spec_plot_210326.png", width = 42, height = 32, units = "cm", scale = 0.65)

##### DCA #####
dca_EUGEI <- dca(True_Label ~ Predicted_Probability,
  data = all_predictions,
  prevalence = 0.00027,
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

dca_all <- dca_EUGEI %>% filter(label != "Predicted_Probability")
dca_EUGEI <- dca_EUGEI %>% filter(label == "Predicted_Probability")

dca_EUGEI$variable <- "EUGEI"
dca_EUGEI$label <- "EUGEI"

dca_all <- rbind(dca_all, dca_EUGEI)

dca_NAPLS <- dca(Transition ~ recalibrated_probs,
  data = df_NAPLS_chr,
  prevalence = 0.00027,
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

dca_NAPLS <- dca_NAPLS %>% filter(label == "recalibrated_probs")
dca_NAPLS$label <- "NAPLS"
dca_NAPLS$variable <- "NAPLS"

dca_all <- rbind(dca_all, dca_NAPLS)
summary_table <- dca_all %>% subset(select = c(variable, threshold, net_benefit))
summary_table_wide <- summary_table %>%
  pivot_wider(names_from = variable, values_from = net_benefit)
summary_table_wide <- summary_table_wide %>% mutate(
  EUGEI = case_when(
    all > 0 ~ EUGEI - all,
    TRUE ~ EUGEI
  ),
  EUGEI_snb = EUGEI / 0.22,
  NAPLS = case_when(
    all > 0 ~ NAPLS - all,
    TRUE ~ NAPLS
  ),
  NAPLS_snb = NAPLS / 0.22
)
write.csv(summary_table_wide, "net_benefit_summary_table_210326_HC.csv", row.names = FALSE)

dca_all$label <- factor(dca_all$label, levels = c("Treat All", "Treat None", "EUGEI", "NAPLS"))
ggplot(data = dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) +
  coord_cartesian(ylim = c(-0.15, 0.075)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(labels = c("Treat All", "Treat None", "EU-GEI", "NAPLS-3"), values = c("gray80", "#000000", "#599ec4", "#c8526a")) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 40), legend.title = element_text(size = 23), legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("dca_all_summary_210326_HC.png", width = 20, height = 15, scale = 0.3)
