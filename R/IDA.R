#' @title Iteration Deflation Algorithm
#' @description Performs Iteration Deflation Algorithm
#' @usage IDA(X, Y, W)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y vector of class probabilities
#' @param W weight matrix where columns represent the \eqn{A} components and rows the \eqn{k} X variables.
#' @author Angela Andreella
#' @return Returns a matrix of scores vectors \code{Tscore}.
#' @export
#' @keywords internal


IDA <- function(X, Y, W){

  #TODO::put check if Y is transformed.


  n <- nrow(X)
  A <- ncol(W)

  E <- R <- Q <- list()

  E[[1]] <- X
  R[[1]] <- Y
  Tscore <- matrix(NA, nrow = n, ncol = A)

  for(a in seq(A)){

    #score vector t_i
    Tscore[,a] <- E[[a]] %*% W[,a]

    #Orthogonal projection matrix of scores
    Q[[a]] <- diag(n) - Tscore[,a] %*% solve(t(Tscore[,a]) %*% Tscore[,a]) %*% t(Tscore[,a])

    #Deflation step

    #X-deflation step
    E[[a+1]] <- Q[[a]] %*% E[[a]] #residual matrix X
    #Y-deflation step
    R[[a+1]] <- Q[[a]] %*% R[[a]] #residual matrix Y
  }

  return(Tscore)
}
