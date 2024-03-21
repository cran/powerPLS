#' @title sample size estimation
#' @description Compute optimal sample size calculation
#' @usage computeSampleSize(n, X, Y, A, alpha, beta,
#' nperm, Nsim, seed, test = "R2",...)
#' @param n vector of sample sizes to consider
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A number of score components
#' @param alpha level of type I error. Default 0.05
#' @param beta level of type II error. Default 0.2.
#' @param Nsim number of simulations. Default 100.
#' @param nperm number of permutations. Default 100.
#' @param seed seed value
#' @param test type of test, one of \code{c("score", "mcc", "R2")}.
#' Default to @R2.
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a dataframe that contains the estimated power for each
#' sample size and number of components considered
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @export
#' @examples
#' \donttest{
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- computeSampleSize(X = datas$X, Y = datas$Y, A = 2, n = 20)
#' }
#' @references
#'
#' Andreella, A., Finos, L., Scarpa, B. and Stocchero, M. "Towards a power analysis for PLS-based methods" 	arXiv:2403.10289 stat.ME.
#'


computeSampleSize <- function(n, X, Y, A, alpha = 0.05, beta = 0.2,
                              nperm = 100, Nsim = 100, seed = 123, test = "R2",...){

  x <-1
samplesize <- foreach(x = seq(length(n)), .combine=cbind) %dopar%
    {
      computePower(X = X, Y = Y, A = A,
                   n = n[x], nperm = nperm, Nsim = Nsim,
                   alpha = alpha, test = test, ...)
    }

  samplesize <- list(size = n,
                     power = samplesize,
                     A = A)

  out <- data.frame(Power = as.vector(t(samplesize$power)),
                    Size = rep(n ),
                    A =rep(seq(A), each =length(n)))


  return(out)
}
