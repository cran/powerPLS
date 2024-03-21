#' @title MCC test
#' @description Performs randomization test based on Matthews Correlation Coefficient
#' @usage mccTest(X, Y, nperm = 100, A, randomization = FALSE, Y.prob = FALSE, eps = 0.01,...)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations. Default 100.
#' @param A number of score components
#' @param randomization Boolean value. Default @FALSE. If @TRUE the permutation p-value is computed
#' @param Y.prob Boolean value. Default @FALSE. IF @TRUE \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a list with the corresponding statistical tests,
#' raw and adjusted p-values
#' @importFrom compositions ilr
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom stats cor
#' @export
#' @seealso The type of tests implemented: \code{\link{scoreTest}} \code{\link{R2Test}}.
#' @author Angela Andreella
#' @return List with the following objects: \code{pv}: raw p-value, \code{pv_adj}: adjusted p-value, \code{test} estimated statistical test.
#' @export
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- mccTest(X = datas$X, Y = datas$Y, A = 1)
#' out




mccTest <- function(X, Y, nperm = 100, A, randomization = FALSE, Y.prob = FALSE, eps = 0.01,...){

  out <- PLSc(X = X, Y = Y, A = A, ...)

  if(!is.null(dim(out$Y_fitted))){
    Y_fitted <- as.factor(out$Y_fitted[,2])
  }else{
    Y_fitted <- as.factor(out$Y_fitted)
  }
  Yf <- as.factor(Y)
  levels(Yf) <-  c(0,1)
  levels(Y_fitted) <-  c(0,1)
  confMatrix <- table(Yf, Y_fitted)
  mcc_obs <- mcc(confMatrix = confMatrix)

  if(randomization){

    null_distr <- foreach(j=seq(nperm-1))%dopar%{

      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      out <- PLSc(X = Xkp, Y = Y, A = A, ...)

      if(!is.null(dim(out$Y_fitted))){
        Y_fitted <- as.factor(out$Y_fitted[,2])
      }else{
        Y_fitted <- as.factor(out$Y_fitted)
      }
      Yf <- as.factor(Y)
      levels(Yf) <- c(0,1)
      levels(Y_fitted) <- c(0,1)
      confMatrix <- table(Yf, Y_fitted)
      mcc_p <- mcc(confMatrix = confMatrix)

    }

    null_distr <-c(mcc_obs, unlist(null_distr))
    pv <- mean(null_distr >= mcc_obs)


  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = mcc_obs))
}
