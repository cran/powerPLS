#' @title Sample size estimation
#' @description Compute optimal sample size
#' @usage computeSampleSize(n, X, Y, A, alpha, beta,
#' nperm, Nsim, seed, test = "R2",...)
#' @param n Vector of sample sizes to consider
#' @param X Data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y Data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A Number of score components
#' @param alpha Type I error level. Default to 0.05
#' @param beta Type II error level. Default to 0.2.
#' @param Nsim Number of simulations. Default to 100.
#' @param nperm Number of permutations. Default to 100.
#' @param seed Seed value
#' @param test Type of test, one of \code{c("score", "mcc", "R2")}.
#' Default to "R2".
#' @param ... Further parameters.
#' @author Angela Andreella
#' @return Returns a data frame that contains the estimated power for each
#' sample size and number of components considered
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @export
#' @examples
#' \dontrun{
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- computeSampleSize(X = datas$X, Y = datas$Y, A = 2, A = 3, n = 20, test = "R2")
#' }
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' Andreella, A., Fino, L., Scarpa, B., & Stocchero, M. (2024). Towards a power analysis for PLS-based methods. arXiv preprint \url{https://arxiv.org/abs/2403.10289}.
#' @seealso \code{\link{computePower}}


computeSampleSize <- function(n, X, Y, A, alpha = 0.05, beta = 0.2,
                              nperm = 100, Nsim = 100, seed = 123, test = "R2",...){

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
