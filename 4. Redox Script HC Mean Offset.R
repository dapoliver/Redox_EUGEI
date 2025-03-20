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
library(predRupdate)
library(tibble)
library(Metrics)
library(rms)
library(predtools)
library(Hmisc)
library(readxl)

set.seed(123)
df <- read_csv("data_survival.csv")
clinical <- read.csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/PPS EU-GEI/Databases/PPS_processed.csv")

df_chr <- df %>% filter(Group != "AtRisk_NoTr")
df_chr <- merge(df_chr, clinical, by.x = "st_subjid", by.y = "ID", all.x = TRUE)

df_chr <- df_chr %>% rename(Gender = Gender.x)
df_cc <- df_chr %>% subset(select = c(MIR132, MIR34A, MIR9, MIR941, MIR137, Transition_status))
df_cc <- df_cc[complete.cases(df_cc), ]
df_cc <- df_cc %>% filter(MIR137 < 75)

data <- df_cc

data$Transition_status <- as.factor((as.character(data$Transition_status)))
levels(data$Transition_status) <- c("Ctrl", "T")

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
outer_folds <- createMultiFolds(data$Transition_status, k = 5, times = 5)
temp_results <- data.frame(
  fold_num = numeric(),
  C = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  balanced_accuracy = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each predictor set
# Loop over each outer fold
for (fold_num in seq_along(outer_folds)) {
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
  }


  train <- train_corrected # %>% subset(select=c(-site))

  x <- model.matrix(~ . - 1, train[, predictors[[1]]])
  y <- train$Transition_status

  # Define fold IDs for cross-validation
  foldid <- sample(rep(1:5, length.out = length(y[train_idx])))

  # Fit model and predict
  cv_model <- cv.glmnet(x, y,
    family = "binomial", alpha = 1,
    nfolds = 5, foldid = foldid
  )

  final_model <- glmnet(x, y,
    family = "binomial", alpha = 1,
    lambda = cv_model$lambda.min
  )

  x_test <- model.matrix(~ . - 1, test[, predictors[[1]]])
  y <- test$Transition_status

  predictions <- as.vector(predict(final_model, newx = x_test, type = "response", s = cv_model$lambda.min))
  PI <- predict(final_model, newx = x_test, type = "link", s = cv_model$lambda.min)

  # Create confusion matrix for a threshold (e.g., 0.5)
  predicted_labels <- ifelse(predictions >= 0.5, 1, 0)
  observed <- as.numeric(as.factor(test$Transition_status)) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")
  model_test <- glm(y ~ PI, family = binomial)

  # Extract performance measures
  temp_results <- rbind(temp_results, data.frame(
    fold_num = fold_num,
    C = concordance(model_test)$concordance,
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    balanced_accuracy = cm$byClass["Balanced Accuracy"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"]
  ))
  # Save predictions
  all_predictions <- rbind(all_predictions, data.frame(
    Subject_ID = test_idx, True_Label = test$Transition_status, Predicted_Probability = predictions,
    Fold = fold_num, Repeat = NA
  ))
}

temp_mean <- temp_results %>% aggregate(. ~ 1, FUN = "mean")
temp_sd <- temp_results %>% aggregate(. ~ 1, FUN = "sd")
results_new <- data.frame(
  C = paste0(
    round(temp_mean$C, 3), " (",
    round(temp_mean$C - 1.96 * (temp_sd$C), 3), "-",
    round(temp_mean$C + 1.96 * (temp_sd$C), 3), ")"
  ),
  balanced_accuracy = paste0(
    round(temp_mean$balanced_accuracy * 100, 1), "% (",
    round(temp_mean$balanced_accuracy * 100 - 1.96 * (temp_sd$balanced_accuracy) * 100, 1), "%-",
    round(temp_mean$balanced_accuracy * 100 + 1.96 * (temp_sd$balanced_accuracy) * 100, 1), "%)"
  ),
  sensitivity = paste0(
    round(temp_mean$sensitivity * 100, 1), "% (",
    round(temp_mean$sensitivity * 100 - 1.96 * (temp_sd$sensitivity) * 100, 1), "%-",
    round(temp_mean$sensitivity * 100 + 1.96 * (temp_sd$sensitivity) * 100, 1), "%)"
  ),
  specificity = paste0(
    round(temp_mean$specificity * 100, 1), "% (",
    round(temp_mean$specificity * 100 - 1.96 * (temp_sd$specificity) * 100, 1), "%-",
    round(temp_mean$specificity * 100 + 1.96 * (temp_sd$specificity) * 100, 1), "%)"
  ),
  ppv = paste0(
    round(temp_mean$PPV * 100, 1), "% (",
    round(temp_mean$PPV * 100 - 1.96 * (temp_sd$PPV) * 100, 1), "%-",
    round(temp_mean$PPV * 100 + 1.96 * (temp_sd$PPV) * 100, 1), "%)"
  ),
  npv = paste0(
    round(temp_mean$NPV * 100, 1), "% (",
    round(temp_mean$NPV * 100 - 1.96 * (temp_sd$NPV) * 100, 1), "%-",
    round(temp_mean$NPV * 100 + 1.96 * (temp_sd$NPV) * 100, 1), "%)"
  )
)

write_csv(results_new, "CV_results_corr_HC.csv")

write.csv(all_predictions, "predictions_output_corr_HC_210225.csv", row.names = FALSE)

predictions_data <- read.csv("predictions_output_corr_HC_210225.csv") # Read the predictions output file
predictions_data <- predictions_data %>% mutate(
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

# Initialize results data frame
results <- data.frame(
  brier = numeric(),
  calibration_intercept = numeric(),
  calibration_slope = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor set
subset_data <- predictions_data

# Average predictions for each subject
averaged_data <- subset_data %>%
  group_by(Subject_ID) %>%
  dplyr::summarize(
    observed = first(True_Label), # True_Label should be the same for each subject
    predicted = mean(Predicted_Probability, na.rm = TRUE),
    .groups = "drop"
  )
# Calibration analysis
calibration <- data.frame(
  observed = averaged_data$observed, # Already in binary format
  predicted = averaged_data$predicted
)

# Fit logistic calibration
logistic_calibration <- predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted)

# Store calibration results
results <- rbind(results, data.frame(
  brier = round(logistic_calibration$BrierScore[1], 2),
  calibration_intercept = round(logistic_calibration$CalInt[1], 2),
  calibration_slope = round(logistic_calibration$CalSlope[1], 2)
))

# Save calibration results
write.csv(results, "calibration_results_corr_HC_210225.csv", row.names = FALSE)
