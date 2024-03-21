#' @title post transformed PLS
#' @description Performs post transformed Partial Least Squares
#' @usage ptPLSc(X, Y, W)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y vector of class probabilities
#' @param W weight matrix where columns represent the \eqn{A} components and rows the \eqn{k} X variables.
#' @author Angela Andreella
#' @return Returns a matrix of scores vectors \code{Tscore}.
#' @export
#' @importFrom nipals nipals
#' @return Returns a list with the following objects:
#' - \code{W}: matrix of weights
#' - \code{6}: post transformation matrix
#' - \code{M}: number of orthogonal components.
#' @references
#'
#' Stocchero, M., & Paris, D. (2016). Post-transformation of PLS2 (ptPLS2) by orthogonal matrix: a new approach for generating predictive and orthogonal latent variables. Journal of Chemometrics, 30(5), 242-251.
#' @keywords internal



ptPLSc <- function(X, Y, W){


  #TODO::put check if Y is transformed.

  A <- ncol(W)
 # out <- svd(t(Y) %*% X %*% W)
  rownames(W) <-NULL
  colnames(X) <- NULL
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  W <- as.matrix(W)
  out <- nipals::nipals(t(Y) %*% X %*% W, center =FALSE, scale = FALSE)
#  V <- out$v
  V <- out$loadings
 # M <- sum(Re(out$d) > 10^-8)
  M <- sum(Re(out$eig) > 10^-8)
  V <- as.matrix(Re(V[,1:M]))
#  d2<-eigen((diag(nrow(V))-V%*%t(V))%*%(t(W)%*%t(X)%*%X%*%W))
  d2 <- nipals::nipals((diag(nrow(V))-V%*%t(V))%*%(t(W)%*%t(X)%*%X%*%W), center =FALSE, scale = FALSE)
  #Compute orthogonal part
 # Go<-Re(d2$vectors[,1:M])
  Go<-Re(d2$loadings[,1:M])
  if(A ==1){Go<- 1}
 # lp <- c()
  G<-Go
  if(A !=M){
  for (i in 1:(A-M)) {
    Xp<-X-X%*%W%*%G%*%solve(t(G)%*%t(W)%*%t(X)%*%X%*%W%*%G)%*%t(G)%*%t(W)%*%t(X)%*%X
   # d3<-eigen(t(diag(A)-Go%*%t(Go))%*%(t(W)%*%t(Xp)%*%Y%*%t(Y)%*%Xp%*%W))
    d3 <- nipals::nipals(t(diag(A)-Go%*%t(Go))%*%(t(W)%*%t(Xp)%*%Y%*%t(Y)%*%Xp%*%W), center =FALSE, scale = FALSE)
   # Gp<-Re(d3$vectors[,1])
    Gp<-Re(d3$loadings[,1])
   # lp[i]<-d3$values[1]
    G<-cbind(G,Gp)
  }
  }

  if((A-M) <= min(qr(Y)$rank, A)){
    warning("The minimum number of predictive latent variables is greater than min(rank(Y), A)")
  }


  #apply G to weight matrix
  Wtilde <- W %*% G #nX times ncomp

  return(list(Wtilde = Wtilde,
              G = G,
              M = M))

}
