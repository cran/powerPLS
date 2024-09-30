#' @title Score test
#' @description Performs permutation-based test based on predictive score vector
#' @usage scoreTest(X, Y, nperm = 200, A, randomization = FALSE,
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
#' @importFrom stats t.test
#' @importFrom stats var
#' @export
#' @seealso Other test statistics implemented: \code{\link{mccTest}} \code{\link{R2Test}}.
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
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- scoreTest(X = datas$X, Y = datas$Y, A = 1)
#' out





scoreTest <- function(X, Y, nperm = 200, A, randomization = FALSE,
                      Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
                      post.transformation = TRUE){

  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling,
              post.transformation = post.transformation,
              eps = eps, Y.prob = Y.prob, transformation = "clr")


  T_score <- out$T_score

  if(!is.na(out$M)){
    M <- out$M
  }else{
    M <- A
  }

  if(A!=M){
    Tp <- T_score[,(M+1):A]
  }else{
    Tp <- T_score
  }

  lev <- unique(as.vector(Y))

  if(length(lev)!=2){stop("Y must be a binary variable")}

  if(is.null(dim(Tp))){
    Tp1 <- as.vector(Tp[Y == lev[1]])
    Tp2 <- as.vector(Tp[Y == lev[2]])
  }else{
    Tp1 <- as.vector(Tp[Y == lev[1],])
    Tp2 <- as.vector(Tp[Y == lev[2],])
  }

  if(length(Tp1)==1){
    var1 <- 0
  }else{
    var1 <- var(Tp1)
  }

  if(length(Tp2)==1){
    var2 <- 0
  }else{
    var2 <- var(Tp2)
  }

  effect_obs <- abs(mean(Tp1) - mean(Tp2)) / sqrt(var1/length(Tp1) + var2/length(Tp2))

  if(randomization){

    null_distr <- replicate(nperm-1,{

    idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
    Xkp <- X[idx,]

    out <- PLSc(X = Xkp, Y = Y, A = A, scaling = scaling,
                post.transformation = post.transformation,
                eps = eps, Y.prob = Y.prob, transformation = "clr")


      T_score <- out$T_score
      T_score

      if(!is.na(out$M)){
        M <- out$M
      }else{
        M <- A
      }

      if(A!=M){
        Tp <- T_score[,(M+1):A]
      }else{
        Tp <- T_score
      }


      if(is.null(dim(Tp))){
        Tp1 <- as.vector(Tp[Y == lev[1]])
        Tp2 <- as.vector(Tp[Y == lev[2]])
      }else{
        Tp1 <- as.vector(Tp[Y == lev[1],])
        Tp2 <- as.vector(Tp[Y == lev[2],])
      }

      if(length(Tp1)==1){
        var1 <- 0
      }else{
        var1 <- var(Tp1)
      }

      if(length(Tp2)==1){
        var2 <- 0
      }else{
        var2 <- var(Tp2)
      }

      abs(mean(Tp1) - mean(Tp2)) / sqrt(var1/length(Tp1) + var2/length(Tp2))

    })

    null_distr <-c(effect_obs, null_distr)
    pv <- mean(null_distr >= effect_obs)
  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = effect_obs))
}
