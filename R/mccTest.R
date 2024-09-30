#' @title MCC test
#' @description Performs permutation-based test based on Matthews Correlation Coefficient
#' @usage mccTest(X, Y, nperm = 200, A, randomization = FALSE,
#' Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
#' post.transformation = TRUE)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations. Default to 200.
#' @param A number of score components
#' @param randomization Boolean value. Default to \code{FALSE}. If \code{TRUE} the permutation p-value is computed
#' @param Y.prob Boolean value. Default \code{FALSE}. IF \code{TRUE} \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector
#' @param scaling Type of scaling, one of
#' \code{c("auto-scaling", "pareto-scaling", "mean-centering")}. Default "auto-scaling".
#' @param post.transformation Boolean value. \code{TRUE} if you want to apply post transformation. Default \code{TRUE}
#' @author Angela Andreella
#' @importFrom compositions ilr
#' @importFrom stats cor
#' @export
#' @seealso Other test statistics implemented: \code{\link{scoreTest}} \code{\link{R2Test}}.
#' @author Angela Andreella
#' @return List with the following objects:
#' \describe{
#'   \item{pv}{raw p-value. It equals \code{NA} if \code{randomization = FALSE}}
#'   \item{pv_adj}{adjusted p-value. It equals \code{NA} if \code{randomization = FALSE}}
#'   \item{test}{estimated test statistic}
#' }
#' @export
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' Andreella, A., Fino, L., Scarpa, B., & Stocchero, M. (2024). Towards a power analysis for PLS-based methods. arXiv preprint \url{https://arxiv.org/abs/2403.10289}.
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 1)
#' out <- mccTest(X = datas$X, Y = datas$Y, A = 1)
#' out




mccTest <- function(X, Y, nperm = 200, A, randomization = FALSE,
                    Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
                    post.transformation = TRUE){

  out <- PLSc(X = X, Y = Y, A = A,
              transformation = "clr",
              scaling = scaling,
              post.transformation = post.transformation,
              eps = eps, Y.prob = Y.prob)

  #Fitted from pilot data
  rownames(out$Y_fitted) <- NULL
  Y_fitted <- as.factor(out$Y_fitted[,2])

  #Observed one
  Yf <- as.factor(Y)

  #check levels
  levels(Yf) <-  c(0,1)
  levels(Y_fitted) <-  c(0,1)

  #confusion matrix
  confMatrix <- table(Yf, Y_fitted)

  #MCC observed
  mcc_obs <- mcc(confMatrix = confMatrix)

  if(randomization){
    null_distr <- replicate(nperm -1, {
      #Permute rows
      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      #Compute Y fitted
      out <- PLSc(X = Xkp, Y = Y, A = A,
                  transformation = "clr",
                  scaling = scaling,
                  post.transformation = post.transformation,
                  eps = eps, Y.prob = Y.prob)

      rownames(out$Y_fitted) <- NULL

      #Compute permuted MCC
      if(length(table(out$Y_fitted))==1){
        Y_fitted <- as.factor(out$Y_fitted)
        levels(Y_fitted) <- c(0,1)
        Y_fitted <- ifelse(Y_fitted == 0, 1, 0)
        lev_drop <- as.numeric(names(table(out$Y_fitted))) -2
        confMatrix <- table(Yf, Y_fitted)

        if(lev_drop== 0){
        confMatrix <- cbind(c(0,0), confMatrix)
        colnames(confMatrix) <- c(0,1)
        mcc(confMatrix = confMatrix)
        }else{
        confMatrix <- cbind(confMatrix,c(0,0))
        colnames(confMatrix) <- c(0,1)
        mcc(confMatrix = confMatrix)
        }

      }else{
        Y_fitted <- as.factor(out$Y_fitted[,2])
        levels(Y_fitted) <- c(0,1)
        confMatrix <- table(Yf, Y_fitted)
        mcc(confMatrix = confMatrix)
      }


    })

    null_distr <-c(mcc_obs, null_distr)
    pv <- mean(null_distr >= mcc_obs)


  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = mcc_obs))
}
