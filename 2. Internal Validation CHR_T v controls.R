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

set.seed(123)

df <- read_csv("data_survival.csv")

df_chr <- df %>% filter(Group != "AtRisk_NoTr")
df_cc <- df_chr %>% subset(select = c(MIR132, MIR34A, MIR9, MIR941, MIR137, Group, Transition_status, day_exit))
df_cc <- df_cc[complete.cases(df_cc), ]
data <- df_cc

data$Transition_status <- as.factor((as.character(data$Transition_status)))
levels(data$Transition_status) <- c("Ctrl", "T")

summary(as.factor(data$Transition_status))
df2 <- df %>% filter(!is.na(MIR132))
summary(as.factor(df2$Transition_status))

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
set.seed(123)

outer_folds <- createMultiFolds(data$Transition_status, k = 5, times = 5)

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

# Loop over each predictor set
x <- model.matrix(~ . - 1, data[, predictors[[1]]])
y <- data$Transition_status

# Loop over each outer fold
for (fold_num in seq_along(outer_folds)) {
  train_idx <- outer_folds[[fold_num]]
  test_idx <- setdiff(seq_len(nrow(data)), train_idx)

  # Define fold IDs for cross-validation
  foldid <- sample(rep(1:5, length.out = length(y[train_idx])))

  set.seed(123)
  # Fit model and predict
  cv_model <- cv.glmnet(x[train_idx, ], y[train_idx],
    family = "binomial", alpha = 1,
    nfolds = 5, foldid = foldid
  )

  set.seed(123)
  final_model <- glmnet(x[train_idx, ], y[train_idx],
    family = "binomial", alpha = 1,
    lambda = cv_model$lambda.min
  )
  predictions <- as.vector(predict(final_model, newx = x[test_idx, ], type = "response"))

  # Create confusion matrix for a threshold (e.g., 0.5)
  predicted_labels <- ifelse(predictions >= 0.5, 1, 0)
  observed <- as.numeric(as.factor(y[test_idx])) - 1
  cm <- confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

  # Extract performance measures
  temp_results <- rbind(temp_results, data.frame(
    fold_num = fold_num,
    C = pROC::auc(roc(observed, predictions)),
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
temp_sd <- temp_results %>% aggregate(. ~ 1, FUN = "sd", na.rm = TRUE)

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
  ),
  posLR = paste0(
    round(temp_mean$LR_pos, 3), " (",
    round(temp_mean$LR_pos - 1.96 * (temp_sd$LR_pos), 3), "-",
    round(temp_mean$LR_pos + 1.96 * (temp_sd$LR_pos), 3), ")"
  ),
  negLR = paste0(
    round(temp_mean$LR_neg, 3), " (",
    round(temp_mean$LR_neg - 1.96 * (temp_sd$LR_neg), 3), "-",
    round(temp_mean$LR_neg + 1.96 * (temp_sd$LR_neg), 3), ")"
  )
)


write_csv(results_new, "CV_results_HC.csv")

# Fit the final model on the full dataset and generate confusion matrix
final_full_model <- glmnet(x, y, family = "binomial", alpha = 1, lambda = cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)$lambda.min)
final_predictions <- as.vector(predict(final_full_model, newx = x, type = "response"))
predicted_class <- ifelse(final_predictions > 0.5, "T", "NT")
confusion_matrix <- table(Predicted = predicted_class, Actual = factor(y, levels = c("NT", "T")))

# Extract coefficients and store
coeffs <- data.frame(
  Variable = rownames(coef(final_full_model)),
  Coefficient = as.vector(coef(final_full_model))
)
coef_results <- rbind(coef_results, coeffs)

write.csv(coef_results, "logistic_regression_LASSO_coefficients_hc_210225.csv", row.names = FALSE)
write.csv(all_predictions, "predictions_output_hc_210225.csv", row.names = FALSE)

predictions_data <- read.csv("predictions_output_hc_210225.csv") # Read the predictions output file
predictions_data <- predictions_data %>% mutate(
  True_Label = case_when(
    True_Label == "T" ~ 1,
    TRUE ~ 0
  ),
  Predicted_Probability = case_when(
    Predicted_Probability > 0.99 ~ 0.99,
    Predicted_Probability < 0.01 ~ 0.01,
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

# Initialize summary table for net benefits
summary_table <- data.frame(
  Threshold = numeric(),
  Net_Benefit = numeric(),
  stringsAsFactors = FALSE
)

# Average predictions for each subject
averaged_data <- predictions_data %>%
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
cal_plot_breaks(calibration, truth = observed, estimate = predicted, smooth = FALSE)

# Store calibration results
results <- rbind(results, data.frame(
  brier = round(logistic_calibration$BrierScore[1], 3),
  calibration_intercept = round(logistic_calibration$CalInt[1], 3),
  calibration_slope = round(logistic_calibration$CalSlope[1], 3)
))

# Decision curve analysis
dca <- data.frame(
  obs = as.numeric(as.factor(averaged_data$observed)) - 1, # Already binary
  pred = averaged_data$predicted
)

dca_assessment <- dca(obs ~ pred,
  data = dca,
  prevalence = 0.000266, # (0.017*0.22),
  thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble()

# Summarize net benefit
dca_assessment <- dca_assessment %>%
  group_by(variable, label, threshold)

write.csv(dca_assessment, paste0("dca_summary_hc_210225.csv"), row.names = FALSE)

# Save calibration results
write.csv(results, "calibration_results_HC_210225.csv", row.names = FALSE)

# Save summary table for net benefits
write.csv(summary_table, "net_benefit_summary_table_HC_210225.csv", row.names = FALSE)

##### DCA summary #####
dca_a <- read_csv("dca_summary_hc_210225.csv")

dca_all <- dca %>% filter(label != "pred")
dca <- dca %>% filter(label == "pred")

dca$variable <- "miRNA"
dca$label <- "miRNA"

dca_all <- rbind(dca_all, dca)

summary_table <- dca_all %>% subset(select = c(variable, threshold, net_benefit))
summary_table_wide <- summary_table %>%
  pivot_wider(names_from = variable, values_from = net_benefit)
summary_table_wide <- summary_table_wide %>% mutate(
  miRNA = case_when(
    all > 0 ~ miRNA - all,
    TRUE ~ miRNA
  ),
  miRNA_snb = miRNA / 0.000266,
)
write.csv(summary_table_wide, "net_benefit_summary_table_hc.csv", row.names = FALSE)

dca_all$label <- factor(dca_all$label, levels = c("Treat All", "Treat None", "miRNA"))
ggplot(data = dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) +
  coord_cartesian(ylim = c(-0.15, 0.075)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values = c("gray80", "#000000", "#599ec4")) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 30), legend.title = element_text(size = 23), legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("dca_all_summary_hc_210225.png", width = 20, height = 15, scale = 0.3)

dca <- read_csv("/Users/domoliver/Library/CloudStorage/Dropbox/Work/Redox EU-GEI/dca_summarya_hc_210225.csv")

dca_all <- dca %>% filter(label != "pred")
dca <- dca %>% filter(label == "pred")

dca$variable <- "EUGEI"
dca$label <- "EUGEI"

dca_all <- rbind(dca_all, dca)
dca_all <- rbind(dca_all, dca_recal_hc)
summary_table <- dca_all %>% subset(select = c(variable, threshold, net_benefit))
summary_table_wide <- summary_table %>%
  pivot_wider(names_from = variable, values_from = net_benefit)
summary_table_wide <- summary_table_wide %>% mutate(
  EUGEI = case_when(
    all > 0 ~ EUGEI - all,
    TRUE ~ EUGEI
  ),
  EUGEI_snb = EUGEI / 0.000266,
  NAPLS = case_when(
    all > 0 ~ NAPLS - all,
    TRUE ~ NAPLS
  ),
  NAPLS_snb = NAPLS / 0.000266
)
write.csv(summary_table_wide, "net_benefit_summary_table_hc_140525.csv", row.names = FALSE)

dca_all$label <- factor(dca_all$label, levels = c("Treat All", "Treat None", "EUGEI", "NAPLS"))
ggplot(data = dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) +
  coord_cartesian(ylim = c(-0.15, 0.075)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(labels = c("Treat All", "Treat None", "EU-GEI", "NAPLS-3"), values = c("gray80", "#000000", "#599ec4", "#c8526a")) +
  theme(text = element_text(family = "Roboto", face = "bold", size = 40), legend.title = element_text(size = 23), legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("dca_all_summary_hc_210525.png", width = 20, height = 15, scale = 0.3)
