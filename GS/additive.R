setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

library(rrBLUP)
library(tidyverse)


# Set seed for reproducibility
set.seed(123)
df_index = read.csv('index_boild_mapping.csv', header = T, na.strings = c("", "."))
df_dosage = read.csv('dosage_boild_converted.csv', na.strings = c(" ", "", "."), header = T, 
                     row.names = 1)
df_index$MGIDI <- as.numeric(scale(df_index$MGIDI))



# performing additive genomic relationship matrix (G) estimation for GS analysis


df_dosage <- as.matrix(df_dosage)  # Exclude the first column (geno)   # convert the dataframe to a matrix

if (!all(df_index$geno %in% rownames(df_dosage))) {
  stop("Error: Some genotypes in index file are not in dosage file")
}
# Ensure df_index is ordered the same as A matrix rows
# df_index <- df_index[match(rownames(df_dosage), df_index$geno), ]

# perform the additive genomic relationship matrix (G) estimation
A <- A.mat(df_dosage)

# Define 10 folds
fold <- 10
n <- dim(df_index)[1]
k <- fold
folds <- sample(rep(1:k, length.out = n))


# i want to extract all the predicted  and actual values for each fold 

actual_values <- list()
predicted_values <- list()
true_values <- list()


# Initialize metrics data frame
matrics <- data.frame(folds = c(1:fold), r = NA, rmse = NA, mae = NA)

# Perform 10-fold cross-validation
for(i in 1:k) {
  # Define test and train indices
  test_idx <- which(folds == i)
  train_idx <- setdiff(1:n, test_idx)
  
  # Create y_train with test set masked as NA
  y_train <- df_index
  y_train$MGIDI[test_idx] <- NA

  # append the actual values to the list
  actual_values[[i]] <- df_index$MGIDI[test_idx]
  
  # Fit GBLUP model using training data only (test set masked as NA)
  # mixed.solve() will automatically handle NA values and predict them
  fit <- mixed.solve(y = y_train$MGIDI, K = A)
  
  # Extract predictions for test set (add intercept)
  pred_test <- fit$u[test_idx] + as.numeric(fit$beta)
  
  # actual values for test set
  true_test <- df_index$MGIDI[test_idx]

  # append the predicted values to the list
  predicted_values[[i]] <- pred_test

  # append the true values to the list
  true_values[[i]] <- true_test
  
  # Calculate metrics
  r <- cor(true_test, pred_test)
  rmse <- sqrt(mean((true_test - pred_test)^2))
  mae <- mean(abs(true_test - pred_test))
  
  # Store metrics
  matrics[i, ] <- c(i, r, rmse, mae)
  print(paste("Fold", i, "completed"))
  
}


# calculate overall metrics using all folds combined
# unlist the lists of vectors into single numeric vectors
all_actual <- unlist(actual_values)
all_predicted <- unlist(predicted_values)
all_true_value <- unlist(true_values)

# Calculate overall metrics across all folds
r <- cor(all_actual, all_predicted)
rmse <- sqrt(mean((all_actual - all_predicted)^2))
mae <- mean(abs(all_actual - all_predicted))

print(paste("Overall R:", round(r, 4)))
print(paste("Overall RMSE:", round(rmse, 4)))
print(paste("Overall MAE:", round(mae, 4)))


# Fit model on ALL data to get population-level (true) heritability
# This uses the complete dataset to estimate variance components
print("\n=== Fitting model on ALL data for population heritability ===")
fit_full <- mixed.solve(y = df_index$MGIDI, K = A)

# Extract variance components from full model
Vu_full <- fit_full$Vu  # Genetic variance (σ²g)
Ve_full <- fit_full$Ve  # Residual variance (σ²e)

# Calculate population heritability
h2_full <- Vu_full / (Vu_full + Ve_full)

print(paste("Genetic Variance (Vu):", round(Vu_full, 6)))
print(paste("Residual Variance (Ve):", round(Ve_full, 6)))
print(paste("Population Heritability (h²):", round(h2_full, 4)))
print(paste("Heritability as percentage:", round(h2_full * 100, 2), "%"))

# Create data frame with only Vu and Ve, using model name as row name
variance_components <- data.frame(
  Vu = Vu_full,
  Ve = Ve_full,
  row.names = "Additive_GBLUP"
)

# Save variance components to CSV (only Vu and Ve)
write.csv(variance_components, "Results/variance_components_Additive_GBLUP.csv", row.names = TRUE)
print("\nVariance components saved to: Results/variance_components_Additive_GBLUP.csv")

# plot the actual vs predicted values make a scatter plot and a line of best fit both dot should have different colors
plot(all_actual, all_predicted, col = "blue", pch = 16, xlab = "Actual", ylab = "Predicted")
fit_line <- lm(all_predicted ~ all_actual)
abline(fit_line, col = "red", lwd = 2, lty = 2) 
title("Actual vs Predicted")
legend("topleft", legend = c("Actual", "Predicted"), col = c("blue", "red"), pch = 16)








# Create Results folder if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}
View(matrics)
write.csv(matrics, paste("Results/additivegblup_", fold, "folds.csv", sep = ""), row.names = FALSE) # save the metrics to a csv file
