setwd('/Users/mac/Documents/Cursor_Project/GS/Tuber_quality_data/boild_data_folder')

library(rrBLUP)
library(tidyverse)
library(BGLR)      # for BGLR() and RKHS model
library(AGHmatrix) # for Gmatrix() with method = "Vitezica"



df_index = read.csv('index_boild_mapping.csv')
df_dosage = read.csv('dosage_boild_converted.csv')



# define the function to calculate the R² score
r2_score <- function(y_true, y_pred) {
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  return(1 - (ss_res / ss_tot))
}




# performing additive genomic relationship matrix (G) estimation for GS analysis

rownames(df_dosage) <- df_dosage$geno
df_dosage <- as.matrix(df_dosage[, -c(1)])  # Exclude the first column (geno)   # convert the dataframe to a matrix

if (!all(df_index$geno %in% rownames(df_dosage))) {
  stop("Error: Some genotypes in index file are not in dosage file")
}



# perform the additive genomic relationship matrix (G) estimation
A <- A.mat(df_dosage)

# perform the dominance genomic relationship matrix (D) estimation
D <- Gmatrix(df_dosage, method = "Vitezica", ploidy = 2)

# define the ETA list
ETA <- list(Add = list(K = A, model = "RKHS"),
            Dom = list(K = D, model = "RKHS"))
# set the seed for reproducibility
set.seed(123)
 predict <- 
# define the number of folds
n <- dim(df_index)[1]
k <- 10 # 10 folds
folds <- sample(rep(1:k, length.out = n))


# initialize the metrics data frame
matrics <- data.frame( folds = c(1:10) , r = NA, r2 = NA, rmse = NA, mae = NA)


# perform the 10-fold cross-validation
for (i in 1:k) {
  
  test_indices <- which(folds == i)
  # define the train indices
  train_indices <- which(folds != i)
  
  y <- df_index$MGIDI
  y[test_indices] <- NA  # Mask the test set phenotypes
  
# nIter = 50000, burnIn = 10000 for publication-quality results
# nIter = 10000, burnIn = 2000 for quick testing


  fm <- BGLR(y = y, ETA = ETA, nIter = 10000, burnIn = 2000, verbose = FALSE)
  predictions <- fm$yHat[test_indices]
  actuals <- df_index$MGIDI[test_indices]
  
  # Calculate prediction accuracy 
  r <- cor(actuals, predictions)
  # calculate the R² score 
  r2 <- r2_score(actuals, predictions)
  # calculate the RMSE
  rmse <- sqrt(mean((predictions - actuals)^2))
  # calculate the MAE
  mae <- mean(abs(predictions - actuals))
  
  matrics[i, ] <- c(i, r, r2, rmse, mae)
  
  print(paste("Fold", i, "completed"))
}

# save the metrics to a csv file
# Create Results folder if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

plot(df_index$MGIDI,EBB)
write.csv(matrics, "Results/10foldcv_additive_dominance_gblup.csv", row.names = FALSE) # save the metrics to a csv file
