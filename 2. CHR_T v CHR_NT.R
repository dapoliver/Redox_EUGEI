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

df_chr <- df %>% filter(Group != "Control")

df_cc <- df_chr %>% subset(select = c(Group, Age, Gender, Ethnicity, MIR132, MIR34A, MIR9, MIR941, MIR137, Transition_status))
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

  x <- model.matrix(~ . - 1, data[, predictors[[1]]])
  y <- data$Transition_status

  # Define fold IDs for cross-validation
  # set.seed(123 + fold_num) # Ensure reproducibility for fold IDs
  foldid <- sample(rep(1:5, length.out = length(y[train_idx])))

  # Fit model and predict
  cv_model <- cv.glmnet(x[train_idx, ], y[train_idx],
    family = "binomial", alpha = 1,
    nfolds = 5, foldid = foldid
  )

  final_model <- glmnet(x[train_idx, ], y[train_idx],
    family = "binomial", alpha = 1,
    lambda = cv_model$lambda.min
  )
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

write.csv(coef_results, "logistic_regression_LASSO_coefficients_111025.csv", row.names = FALSE)
write.csv(all_predictions, "predictions_output_111025.csv", row.names = FALSE)

predictions_data <- read.csv("predictions_output_111025.csv") # Read the predictions output file
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

# Initialize summary table for net benefits
summary_table <- data.frame(
  Threshold = numeric(),
  Net_Benefit = numeric(),
  stringsAsFactors = FALSE
)

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
  brier = paste0(
    logistic_calibration$BrierScore[1], " (",
    round(logistic_calibration$Brier_lower[1], 2), "-",
    round(logistic_calibration$Brier_upper[1], 2), ")"
  ),
  calibration_intercept = paste0(
    round(logistic_calibration$CalInt[1], 2), " (",
    round(logistic_calibration$CalInt_lower[1], 2), "-",
    round(logistic_calibration$CalInt_upper[1], 2), ")"
  ),
  calibration_slope = paste0(
    round(logistic_calibration$CalSlope[1], 2), " (",
    round(logistic_calibration$CalSlope_lower[1], 2), "-",
    round(logistic_calibration$CalSlope_upper[1], 2), ")"
  )
)

write_csv(results_new, "CV_results_111025.csv")

# Save calibration plot
png(paste0("calibration_plot_171025.png"), width = 600, height = 500)
print(cal_plot_logistic(calibration, truth = observed, estimate = predicted, smooth = FALSE))
dev.off()

# Save calibration results
write.csv(results, "calibration_results_111025.csv", row.names = FALSE)

##### External Validation #####

df_NAPLS <- read_excel("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Papers/Submitted/Redox EU-GEI/NAPLS/NAPLS BATCH info17102025.xlsx")

summary(factor(df_NAPLS$`GROUP (UC = CTRL group)`))
df_NAPLS <- df_NAPLS %>% subset(select = c(
  demo_age_ym, demo_sex, demo_racial, `miR-9`, `miR-34`, `miR-132`, `miR-137`, `miR-941`, `GROUP (UC = CTRL group)`, SiteNumber,
  P1_SOPS, P2_SOPS, P3_SOPS, P4_SOPS, P5_SOPS, GlobalAssessmentFunction
))

df_NAPLS <- df_NAPLS %>%
  dplyr::rename(MIR9 = `miR-9`, MIR34A = `miR-34`, MIR132 = `miR-132`, MIR137 = `miR-137`, MIR941 = `miR-941`)
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
write.csv(results_NAPLS, "external_validation_results_171025.csv", row.names = FALSE)

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
write.csv(results_NAPLS_recal, "external_validation_results_recal_171025.csv", row.names = FALSE)

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

for (i in seq(0.01, 1, by = 0.01)) {
  predicted_labels <- ifelse(recalibrated_probs >= i, 1, 0)
  observed <- as.numeric(as.factor(df_NAPLS_chr$Transition)) - 1
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
ggsave("likelihood_ratio_plot_171025.png", width = 42, height = 32, units = "cm", scale = 0.65)

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
ggsave("ppv_npv_plot_171025.png", width = 42, height = 32, units = "cm", scale = 0.65)

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
ggsave("sens_spec_plot_171025.png", width = 42, height = 32, units = "cm", scale = 0.65)

##### DCA summary #####
dca <- read_csv("dca_summary_111025.csv")

dca_all <- dca %>% filter(label != "pred")
dca <- dca %>% filter(label == "pred")

dca$variable <- "EUGEI"
dca$label <- "EUGEI"

dca_all <- rbind(dca_all, dca)

dca_recal_chr <- dca(observed ~ predicted,
  data = calibration_recal,
  prevalence = 0.22,
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

dca_recal_chr <- dca_recal_chr %>% filter(label == "predicted")
dca_recal_chr$label <- "NAPLS"
dca_recal_chr$variable <- "NAPLS"

dca_all <- rbind(dca_all, dca_recal_chr)
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
write.csv(summary_table_wide, "net_benefit_summary_table_111025.csv", row.names = FALSE)

dca_all$label <- factor(dca_all$label, levels = c("Treat All", "Treat None", "EUGEI", "NAPLS"))
ggplot(data = dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) +
  coord_cartesian(ylim = c(-0.005, 0.25)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(labels = c("Treat All", "Treat None", "EU-GEI", "NAPLS-3"), values = c("gray80", "#000000", "#599ec4", "#c8526a")) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 40), legend.title = element_text(size = 23), legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("dca_all_summary_171025.png", width = 20, height = 15, scale = 0.3)
