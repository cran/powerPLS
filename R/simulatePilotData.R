#' @title simulate pilot data
#' @description simulate cluster pilot data
#' @usage simulatePilotData(seed = 123, nvar, clus.size, nvar_rel,m, A = 2, S1 = NULL, S2 = NULL)
#' @param seed seed value
#' @param nvar number of variables
#' @param clus.size Vector of two elements, specifying the size of classes (only two classes are considered)
#' @param nvar_rel number of variables relevant to predict \code{Y}
#' @param m separation between classes
#' @param A oracle number of score components
#' @param S1 covariance matrix for the first class. Default @NULL i.e., the identity is considered.
#' @param S2 covariance matrix for the second class. Default @NULL i.e., the identity is considered.
#' @importFrom MASS mvrnorm
#' @importFrom stats princomp
#' @importFrom stats rnorm
#' @author Angela Andreella
#' @return Returns list of X and Y simulated data
#' @export
#' @examples
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)

simulatePilotData <- function(seed = 123, nvar, clus.size, nvar_rel,m, A = 2, S1 = NULL, S2 = NULL){

  set.seed(seed)
  n <- clus.size[1] + clus.size[2]
  if(is.null(S1)){S1 = diag(A)}
  if(is.null(S2)){S2 = diag(A)}

  X<- rbind(mvrnorm(n = clus.size[1], mu = rep(0, A), Sigma = S1),
            mvrnorm(n = clus.size[2], mu = rep(m, A), Sigma = S2))




  X <- scale(X, scale = FALSE)
  #Predictive latent variable n x A
  out <- princomp(X)$scores #rank

  #loading matrix nvar_rel x
  out1 <- princomp(data.frame(matrix(runif(A*nvar_rel), ncol = A)),
                   scores = TRUE)$scores


  X <- cbind(out %*% t(out1),matrix(rnorm(n*(nvar-nvar_rel)), ncol = nvar-nvar_rel))

  # X <- X + c(rep(m,clus.size[1]), rep(0, clus.size[2]))

  Y<- c(rep(0, clus.size[1]), rep(1, clus.size[2]))

  return(simData = list(X = X, Y = Y))
}
