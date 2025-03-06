#' @title Repeated k-Fold Cross-Validation with Custom Test Metrics
#'
#' @description This function performs repeated k-fold cross-validation and computes a selected performance metric across all repetitions and folds.
#' It allows for different types of performance tests, such as MCC, sensitivity, specificity, R2, F1, and more.
#'
#' @param data A data frame or matrix of features (predictor variables).
#' @param labels A vector of class labels corresponding to the rows of \code{data}.
#' @param k_folds An integer specifying the number of cross-validation folds (default = 5).
#' @param repeats An integer specifying the number of times the cross-validation is repeated (default = 3).
#' @param A number of score components
#' @param test_type A character string specifying the type of test to use. Options include:
#'   \itemize{
#'     \item 'mccTest' for Matthews Correlation Coefficient (MCC),
#'     \item 'sensitivityTest' for Sensitivity,
#'     \item 'specificityTest' for Specificity,
#'     \item 'R2Test' for R-squared,
#'     \item 'scoreTest' for Score,
#'     \item 'F1Test' for F1 Score,
#'     \item 'FMTest' for Fowlkes-Mallows Index (FM),
#'     \item 'AUCTest' for Area Under the Curve (AUC),
#'     \item 'dQ2Test' for dQ2.
#'   }
#'   Default is 'mccTest'.
#' @param seed An integer for setting the random seed to ensure reproducibility (default = 1234).
#'
#' @return A numeric value representing the average performance metric across the outer folds.
#' @import caret
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(15,15),m = 6,nvar_rel = 5,A = 1)
#' data <- datas$X
#' labels <- datas$Y
#' mean_mcc <- repeatedCV_test(data, labels, A = 1, test_type = 'mccTest')
#' cat('Mean MCC:', mean_mcc, '\n')
#'
#' mean_score <- repeatedCV_test(data, labels, A = 1, test_type = 'scoreTest')
#' cat('Mean Sensitivity:', mean_score, '\n')
#'
#' @export

repeatedCV_test <- function(data, labels, k_folds = 5, repeats = 3, A = 1, test_type = "mccTest",
                            seed = 1234) {
  set.seed(seed)

  labels <- as.factor(labels)

  # Function to compute the selected test metric
  compute_test <- function(X, Y, A, test_type) {
    if (test_type == "mccTest") {
      return(mccTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "sensitivityTest") {
      return(sensitivityTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "specificityTest") {
      return(specificityTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "R2Test") {
      return(R2Test(X = X, Y = Y, A = A)$test)
    } else if (test_type == "scoreTest") {
      return(scoreTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "F1Test") {
      return(F1Test(X = X, Y = Y, A = A)$test)
    } else if (test_type == "FMTest") {
      return(FMTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "AUCTest") {
      return(AUCTest(X = X, Y = Y, A = A)$test)
    } else if (test_type == "dQ2Test") {
      return(dQ2Test(X = X, Y = Y, A = A)$test)
    } else {
      stop("Test type not recognized. Please choose a valid test type.")
    }
  }

  # Initialize results storage
  all_results <- numeric(k_folds * repeats)
  result_index <- 1

  for (rep in 1:repeats) {
    cat("Repetition", rep, "of", repeats, "\n")

    # Create folds for this repetition
    fold_indices <- createFolds(labels, k = k_folds, list = TRUE)

    for (fold in seq_along(fold_indices)) {
      cat("  Fold", fold, "of", k_folds, "\n")

      test_indices <- fold_indices[[fold]]
      train_indices <- setdiff(seq_len(nrow(data)), test_indices)

      train_data <- data[train_indices, ]
      test_data <- data[test_indices, ]

      train_labels <- labels[train_indices]
      test_labels <- labels[test_indices]

      # Compute test metric on the test set
      test_value <- compute_test(X = test_data, Y = test_labels, A = A, test_type = test_type)

      # Store the result
      all_results[result_index] <- test_value
      result_index <- result_index + 1
    }
  }

  # Compute the average performance metric across all folds and repetitions
  mean_test_value <- mean(all_results)
  return(mean_test_value)
}
