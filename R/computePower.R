#' @title Power estimation
#' @description estimate power for a given sample size, alpha level and number
#' of score components.
#' @usage computePower(X, Y, A, n, seed = 123,
#' Nsim = 100, nperm = 200, alpha = 0.05,
#' test = "R2", Y.prob = FALSE, eps = 0.01, ...)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A number of score components
#' @param n sample size
#' @param seed seed value
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha type I error
#' @param test type of test, one of \code{c("score", "mcc", "R2")}.
#' @param Y.prob Boolean value. Default @FALSE. IF @TRUE \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector.
#' Default to "R2".
#' @param ... Futher parameters see \code{\link{PLSc}}
#' @author Angela Andreella
#' @return Returns the corresponding estimated power
#' @export
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @examples
#' \donttest{
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- computePower(X = datas$X, Y = datas$Y, A = 3, n = 20)
#' }
#' @references
#'
#' Andreella, A., Finos, L., Scarpa, B. and Stocchero, M. "Towards a power analysis for PLS-based methods" 	arXiv:2403.10289 stat.ME.
#'



computePower <- function(X, Y, A, n, seed = 123,
                         Nsim = 100, nperm = 200, alpha = 0.05,
                         test = "R2", Y.prob = FALSE, eps = 0.01, ...){

  if(any(!(test %in% c("R2", "mcc", "score")))){
    stop("available tessts are R2, mcc and score")
  }

  #Build the reference model PLS2c

  outPLS <- PLSc(X = X, Y = Y, A = A, Y.prob = Y.prob, eps = eps, ...)

  pw <- matrix(0, ncol = length(test), nrow = A)
  colnames(pw)<- test

i <- NULL
pw <-  foreach(i = c(1:Nsim)) %dopar% {

    #Model the distribution of the X-data
    outsim <- sim_XY(out = outPLS, n = n, seed = 1234+i, A = A, ...)
    #Model the distribution of the Y-data
    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1


    #Apply one test
    if(length(test) == 1){

      if(test == "mcc"){

        pv <- foreach(x = seq(A), .combine=rbind) %dopar%{
          mccTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                  randomization = TRUE, ...)
        }

      }
      if(test == "score"){
        pv <- foreach(x = seq(A), .combine=rbind) %dopar%{
          scoreTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                  randomization = TRUE, ...)
        }

      }
      if(test == "R2"){
        pv <- foreach(x = seq(A), .combine=rbind) %dopar%{
          R2Test(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    randomization = TRUE, ...)
        }

      }

      pv <- as.data.frame(pv)
      pv <- data.frame(pv = unlist(pv$pv),
                       pv_adjust = unlist(pv$pv_adj))
      for(x in seq(A)){
        if(pv$pv_adj[x] <= alpha){pw[x] <- pw[x] + 1}
      }
    }else{

      #Apply more than one test.

      if("mcc" %in% test){
        pv_mcc <- foreach(x = seq(A), .combine=cbind) %dopar%{
          mccTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    randomization = TRUE, ...)
        }

      }
      if("score" %in% test){
        pv_score <- foreach(x = seq(A), .combine=cbind) %dopar%{
          scoreTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                  randomization = TRUE, ...)
        }

      }
      if("R2" %in% test){
        pv_R2 <- foreach(x = seq(A), .combine=cbind) %dopar%{
          R2Test(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                  randomization = TRUE, ...)
        }

      }

      names_test <- ls()[ls() %in% c("pv_mcc", "pv_score", "pv_R2")]

      pv_out <- t(sapply(seq(length(names_test)),
                         function(x) eval(as.name(names_test[x]))[2,]))

      colnames(pv_out) <- gsub("pv_", "", names_test)
      rownames(pv_out) <- seq(A)

      for(x in seq(A)){
        for(y in seq(length(test))){
          if(pv_out[x,y] <= alpha){
            pw[x,y] <- pw[x,y] + 1
          }
        }

      }
    }


    pw


  }


pw <- pw[[Nsim]]
pw <- pw/Nsim


  return(pw)
}
