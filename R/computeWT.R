#' @title compute weight and score matrices from PLSc
#' @description compute weight and score matrices for Partial Least Squares classification
#' @usage computeWT(X, Y, A)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A number of score components
#' @importFrom nipals nipals
#' @author Angela Andreella
#' @return Returns a list with the following objects:
#' - \code{W}: matrix of weights
#' - \code{T_score}: matrix of Y scores
#' - \code{R}: matrix of Y residuals
#' @export
#' @keywords internal



computeWT <- function(X, Y, A){



  E <- R <- w <- r <- Q <- list()

  E[[1]] <- X
  R[[1]] <- Y
  #number of observations
  N <- nrow(Y)

  #dim(W) = nX times A

  for(i in seq(A)){

    AA <- t(E[[i]]) %*% R[[i]] %*% t(R[[i]]) %*% E[[i]]

    #Compute weight matrix
    #out <- eigen(AA) #well defined eigenvalue problem
    out <- nipals(AA, center =FALSE, scale = FALSE)
   # w[[i+1]] <- Re(out$vectors[,1])
    w[[i+1]] <- Re(out$loadings[,1])

    #score vector t_i
    r[[i+1]] <- E[[i]] %*% w[[i+1]]

    #Orthogonal projection matrix of scores
    #Q[[i+1]] <- diag(N) - r[[i+1]] %*% t(r[[i+1]])/ (t(r[[i+1]]) %*% r[[i+1]])[1]
    Q[[i+1]] <- diag(N) - r[[i+1]] %*% solve(t(r[[i+1]]) %*% r[[i+1]]) %*% t(r[[i+1]])

    #Deflation step
    #X-deflation step
    E[[i+1]] <- Q[[i+1]] %*% E[[i]] #residual matrix X
    #Y-deflation step
    R[[i+1]] <- Q[[i+1]] %*% R[[i]] #residual matrix Y
  }

  #Rearrange weight matrix
  W <- NULL
  for (i in c(2:(length(w)))) W <- cbind(W, w[[i]]) #number of obs times A


  #Rearrange scores matrix
  T_score <- NULL
  for (i in c(2:(length(r)))) T_score <- cbind(T_score, r[[i]]) #number of obs times ncomponent


  return(list(W = W,
              T_score = T_score,
              R = R))
}
