setwd("~/Dropbox/Work/Redox EU-GEI/")

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

df = read_csv("data_survival.csv")

df_chr = df %>% filter(Group!="Ctrl")
df_chr = df_chr %>%
  mutate(Transition_status = ifelse(is.na(Transition_status) & Group == "control", 0, Transition_status))
df_cc = df_chr %>% subset(select=c(MIR132, MIR34A, MIR9, MIR941, MIR137, Age, Gender, Ethnicity,Transition_status, day_exit))
df_cc = df_cc[complete.cases(df_cc),]
df_cc = df_cc %>% mutate(Gender = as.factor(Gender),
                         Ethnicity = as.factor(Ethnicity))
data = df_cc

data$Transition_status = as.factor((as.character(data$Transition_status)))
levels(data$Transition_status) = c("NT", "T")


# Define the predictor sets
predictors_cox = list(a=c("MIR132", "MIR34A", "MIR9", "MIR941"),
                      #b=c("MIR132", "MIR34A", "MIR9", "MIR941","MIR137"),
                      c=c("MIR132", "MIR34A", "MIR9", "MIR941","Age","Gender","Ethnicity")
                      #d=c("MIR132", "MIR34A", "MIR9", "MIR941","MIR137","Age","Gender","Ethnicity")
)

# Custom function to calculate Brier Score
brier_score = function(observed, predicted) {
  mean((observed - predicted) ^ 2)
}

# Create data frames to store predictions and coefficients
all_predictions = data.frame(Subject_ID = integer(), True_Label = integer(),
                             Predicted_Probability = numeric(), Predictor_Set = character(),
                             Fold = integer(), Repeat = integer(), stringsAsFactors = FALSE)

coef_results = data.frame(Predictor_Set = character(), Variable = character(), Coefficient = numeric(), stringsAsFactors = FALSE)

# Create 5-fold cross-validation with 5 repeats for outer folds
outer_folds = createMultiFolds(data$Transition_status, k = 5, times = 5)
temp_results = data.frame(predictors = character(), 
                          fold_num = numeric(),
                          C = numeric(), 
                          sensitivity = numeric(),
                          specificity = numeric(),
                          balanced_accuracy = numeric(),
                          PPV = numeric(),
                          NPV = numeric(),
                          precision = numeric(),
                          recall = numeric(),
                          F1_Score = numeric(),
                          Brier = numeric(),
                          Calibration_Intercept = numeric(),
                          Calibration_Slope = numeric(),
                          stringsAsFactors = FALSE)
# Loop over each predictor set
for (i in 1:2) {
  x = model.matrix(~.-1,data[, predictors_cox[[i]]])
  y = data$Transition_status

  # Loop over each outer fold
  for (fold_num in seq_along(outer_folds)) {
    train_idx = outer_folds[[fold_num]]
    test_idx = setdiff(seq_len(nrow(data)), train_idx)
    
    # Fit model and predict
    final_model = glmnet(x[train_idx, ], y[train_idx], family = "binomial", alpha = 1,
                         lambda = cv.glmnet(x[train_idx, ], y[train_idx], family = "binomial", alpha = 1, nfolds = 5)$lambda.min)
    predictions = as.vector(predict(final_model, newx = x[test_idx, ], type = "response"))

    # Create confusion matrix for a threshold (e.g., 0.5)
    predicted_labels = ifelse(predictions >= 0.5, 1, 0)
    observed = as.numeric(as.factor(y[test_idx]))-1
    cm = confusionMatrix(as.factor(predicted_labels), as.factor(observed), positive = "1")

    # Extract performance measures
    temp_results <- rbind(temp_results, data.frame(predictors = i, 
                                                   fold_num = fold_num,
                                                   C = pROC::auc(roc(observed, predictions)),
                                                   sensitivity = cm$byClass["Sensitivity"],
                                                   specificity = cm$byClass["Specificity"],
                                                   balanced_accuracy = cm$byClass["Balanced Accuracy"],
                                                   PPV = cm$byClass["Pos Pred Value"],
                                                   NPV = cm$byClass["Neg Pred Value"],
                                                   precision = cm$byClass["Precision"],
                                                   recall = cm$byClass["Recall"],
                                                   F1_Score = cm$byClass["F1"],
                                                   Brier = brier_score(observed, predictions)))
    # Save predictions
    all_predictions = rbind(all_predictions, data.frame(
      Subject_ID = test_idx, True_Label = y[test_idx], Predicted_Probability = predictions,
      Predictor_Set = names(predictors_cox)[i], Fold = fold_num, Repeat = NA
    ))
  }

    temp_mean <- temp_results %>% aggregate(.~predictors, FUN="mean")
    temp_sd <- temp_results %>% aggregate(.~predictors, FUN="sd")
    results_new <- data.frame(predictors=c(1:4),
                              C=paste0(round(temp_mean$C,3)," (",
                                       round(temp_mean$C-1.96*(temp_sd$C),3),"-",
                                       round(temp_mean$C+1.96*(temp_sd$C),3), ")"),
                              balanced_accuracy=paste0(round(temp_mean$balanced_accuracy*100,1),"% (",
                                                       round(temp_mean$balanced_accuracy*100-1.96*(temp_sd$balanced_accuracy)*100,1),"%-",
                                                       round(temp_mean$balanced_accuracy*100+1.96*(temp_sd$balanced_accuracy)*100,1), "%)"),
                              sensitivity=paste0(round(temp_mean$sensitivity*100,1),"% (",
                                                       round(temp_mean$sensitivity*100-1.96*(temp_sd$sensitivity)*100,1),"%-",
                                                       round(temp_mean$sensitivity*100+1.96*(temp_sd$sensitivity)*100,1), "%)"),
                              specificity=paste0(round(temp_mean$specificity*100,1),"% (",
                                                 round(temp_mean$specificity*100-1.96*(temp_sd$specificity)*100,1),"%-",
                                                 round(temp_mean$specificity*100+1.96*(temp_sd$specificity)*100,1), "%)"),
                              ppv=paste0(round(temp_mean$PPV*100,1),"% (",
                                                 round(temp_mean$PPV*100-1.96*(temp_sd$PPV)*100,1),"%-",
                                                 round(temp_mean$PPV*100+1.96*(temp_sd$PPV)*100,1), "%)"),
                              npv=paste0(round(temp_mean$NPV*100,1),"% (",
                                                 round(temp_mean$NPV*100-1.96*(temp_sd$NPV)*100,1),"%-",
                                                 round(temp_mean$NPV*100+1.96*(temp_sd$NPV)*100,1), "%)"),
                              precision=paste0(round(temp_mean$precision*100,1),"% (",
                                                 round(temp_mean$precision*100-1.96*(temp_sd$precision)*100,1),"%-",
                                                 round(temp_mean$precision*100+1.96*(temp_sd$precision)*100,1), "%)"),
                              recall=paste0(round(temp_mean$recall*100,1),"% (",
                                                 round(temp_mean$recall*100-1.96*(temp_sd$recall)*100,1),"%-",
                                                 round(temp_mean$recall*100+1.96*(temp_sd$recall)*100,1), "%)"),
                              F1=paste0(round(temp_mean$F1_Score*100,1),"% (",
                                                 round(temp_mean$F1_Score*100-1.96*(temp_sd$F1_Score)*100,1),"%-",
                                                 round(temp_mean$F1_Score*100+1.96*(temp_sd$F1_Score)*100,1), "%)"),
                              Brier=paste0(round(temp_mean$Brier*100,1),"% (",
                                                 round(temp_mean$Brier*100-1.96*(temp_sd$Brier)*100,1),"%-",
                                                 round(temp_mean$Brier*100+1.96*(temp_sd$Brier)*100,1), "%)")
                              )

    write_csv(results_new, "CV_results.csv")
    
  # Fit the final model on the full dataset and generate confusion matrix
  final_full_model = glmnet(x, y, family = "binomial", alpha = 1, lambda = cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)$lambda.min)
  final_predictions = as.vector(predict(final_full_model, newx = x, type = "response"))
  predicted_class = ifelse(final_predictions > 0.5, "T", "NT")
  confusion_matrix = table(Predicted = predicted_class, Actual = factor(y, levels = c("NT", "T")))
  
  # Save confusion matrix to CSV
  write.csv(as.data.frame(confusion_matrix), paste0("confusion_matrix_", names(predictors_cox)[i], ".csv"), row.names = TRUE)
  
  # Extract coefficients and store
  coeffs = data.frame(Predictor_Set = names(predictors_cox)[i], 
                      Variable = rownames(coef(final_full_model)), 
                      Coefficient = as.vector(coef(final_full_model)))
  coef_results = rbind(coef_results, coeffs)
}
write.csv(coef_results, "logistic_regression_LASSO_coefficients.csv", row.names = FALSE)
write.csv(all_predictions, "predictions_output.csv", row.names = FALSE)

predictions_data = read.csv("predictions_output.csv") # Read the predictions output file
predictions_data$True_Label = ifelse(predictions_data$True_Label == "T", 1, 0) # Recode True_Label to binary format
predictor_sets = unique(predictions_data$Predictor_Set) # Get unique predictor sets

# Initialize results data frame
results = data.frame(Predictor_Set = character(), 
                     calibration_intercept = numeric(), 
                     calibration_slope = numeric(),
                     stringsAsFactors = FALSE)

# Initialize summary table for net benefits
summary_table = data.frame(Predictor_Set = character(),
                           Threshold = numeric(),
                           Net_Benefit = numeric(),
                           stringsAsFactors = FALSE)

# Loop through each predictor set
for (set in predictor_sets) {
  
  # Filter data for the current predictor set
  subset_data = predictions_data %>% filter(Predictor_Set == set)
  
  # Average predictions for each subject
  averaged_data = subset_data %>%
    group_by(Subject_ID) %>%
    dplyr::summarize(observed = first(True_Label),  # True_Label should be the same for each subject
                     predicted = mean(Predicted_Probability, na.rm = TRUE),
                     .groups = 'drop')
  averaged_data <- averaged_data %>% mutate(predicted=case_when(predicted>0.99 ~ 0.99,
                                                                predicted<0.001 ~ 0.001,
                                            TRUE ~ predicted))
  # Calibration analysis
  calibration = data.frame(observed = as.numeric(as.factor(averaged_data$observed))-1,  # Already in binary format
                           predicted = averaged_data$predicted)
  
  # Fit logistic calibration
  logistic_calibration = predRupdate::pred_val_probs(binary_outcome = calibration$observed, Prob = calibration$predicted-0.0001)
  
  # Store calibration results
  results = rbind(results, data.frame(Predictor_Set = set,
                                      calibration_intercept = round(logistic_calibration$CalInt[1], 2),
                                      calibration_slope = round(logistic_calibration$CalSlope[1], 2)))
  
  # Generate calibration plot
  cal_plot = ggplot(calibration, aes(x = predicted, y = observed)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Calibration Plot for", set),
         x = "Predicted Probability",
         y = "Observed Probability") +
    theme_minimal()
  
  # Save calibration plot
  ggsave(paste0("calibration_plot_", set, ".png"), plot = cal_plot, width = 6, height = 5, dpi = 300)
  
  # Decision curve analysis
  dca = data.frame(obs = as.numeric(as.factor(averaged_data$observed))-1,  # Already binary
                   pred = averaged_data$predicted)
  
  dca_assessment = dca(obs ~ pred,
                       data = dca,
                       prevalence = 0.017*0.22, 
                       thresholds = seq(0, 0.5, 0.01)) %>%
    as_tibble()
  
  # Summarize net benefit
  dca_assessment = dca_assessment %>%
    group_by(variable, label, threshold) 
  
  write.csv(dca_assessment, paste0("dca_summary", set,".csv"), row.names = FALSE)
  
}

# Save calibration results
write.csv(results, "calibration_results.csv", row.names = FALSE)

# Save summary table for net benefits
write.csv(summary_table, "net_benefit_summary_table.csv", row.names = FALSE)

##### DCA summary #####
dca_b <- read_csv("dca_summarya.csv")
dca_d <- read_csv("dca_summaryc.csv")

dca_all <- dca_b %>% filter(label!="pred")
dca_b <- dca_b %>% filter(label=="pred")
dca_d <- dca_d %>% filter(label=="pred")

dca_b$variable <- "miRNA"
dca_d$variable <- "miRNA+demo"
dca_b$label <- "miRNA"
dca_d$label <- "miRNA+demo"

dca_all <- rbind(dca_all, dca_b)
dca_all <- rbind(dca_all, dca_d)

dca_all$label <- factor(dca_all$label,levels = c("Treat All", "Treat None", "miRNA", "miRNA+demo"))
ggplot(data=dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.5) +
  coord_cartesian(ylim = c(-0.005, 0.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values=c("gray80","#000000","#599ec4","#c8526a")) +
  theme(text = element_text(family="Roboto", face="bold", size=30),legend.title = element_text(size=23),legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("dca_all_summary_091024.png", width=20, height=15, scale=0.5)
