#' @title PLS classification
#' @description Performs Partial Least Squares classification
#' @usage PLSc(X, Y, A, scaling = "auto-scaling", post.transformation = TRUE,
#' eps = 0.01, Y.prob = FALSE, transformation = "ilr")
#' @param X Data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y Data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A Number of score components
#' @param scaling Type of scaling, one of
#' \code{c("auto-scaling", "pareto-scaling", "mean-centering")}. Default to "auto-scaling"
#' @param post.transformation Boolean value. \code{TRUE} if you want to apply post transformation. Default \code{TRUE}
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector
#' @param Y.prob Boolean value. Default \code{FALSE}. IF \code{TRUE} \code{Y} is a probability vector
#' @param transformation Transformation used to map \code{Y} in probability data vector. The options are "ilr" and "clr".
#' Default @ilr.
#' @author Angela Andreella
#' @return List with the following objects:
#' \describe{
#'   \item{W}{Matrix of weights}
#'   \item{X_loading}{Matrix of \code{X} loading}
#'   \item{Y_loading}{Matrix of \code{Y} loading}
#'   \item{X}{Matrix of \code{X} data (predictor variables)}
#'   \item{Y}{Matrix of \code{Y} data (dependent variable)}
#'   \item{T_score}{Matrix of scores}
#'   \item{Y_fitted}{Fitted \code{Y} matrix}
#'   \item{B}{Matrix regression coefficients}
#'   \item{M}{Number of orthogonal components if \code{post.transformation=TRUE} is applied.}
#'   }
#' @importFrom compositions ilrInv
#' @importFrom compositions ilr
#' @importFrom compositions clr
#' @importFrom compositions clrInv
#' @importFrom stats sd
#' @importFrom stats model.matrix
#'
#' @export
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- PLSc(X = datas$X, Y = datas$Y, A = 3)
#'
#' @references Stocchero, M., De Nardi, M., & Scarpa, B. (2021). PLS for classification. Chemometrics and Intelligent Laboratory Systems, 216, 104374.


PLSc <- function(X, Y, A, scaling = "auto-scaling", post.transformation = TRUE,
                 eps = 0.01, Y.prob = FALSE, transformation = "ilr"){

  nY <- ifelse(is.null(dim(Y)), length(Y), dim(Y)[1])
  if(dim(X)[1] != nY){
    stop("X and Y must have the same number of observations!")
    }

  if(!(scaling %in% c("auto-scaling", "pareto-scaling", "mean-centering"))){
    stop("available scaling are auto-scaling, pareto-scaling and mean-centering")
  }
  X<-as.matrix(X)
  if(scaling == "auto-scaling"){
    Mm <- apply(X, 2, mean)
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "pareto-scaling"){
    Mm <- 0
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "mean-centering"){
    Mm <- apply(X, 2, mean)
    s <- 1
    X <- (X - Mm)/s
  }

  #If Y is not probability but a vector of classes
  if(!Y.prob){

    if(is.null(dim(Y)) | ncol(as.matrix(Y))==1){
      Y <- as.matrix(Y)

      if(!is.factor(Y)){
        Y <- as.factor(Y)

      }
      levels(Y) <- c(0,1)
      Y <- model.matrix(~0+Y)
    }

    #Transform to probability matrix
    Y[which(Y==0)]<-eps
    Y[which(Y==1)]<-1-(ncol(Y)-1)*eps

    #Centered log ratio transform transformation
    if(transformation == "clr"){
      P <- matrix(clr(Y), ncol = ncol(Y))
    }
    if(transformation == "ilr"){
    P <- matrix(ilr(Y), ncol = 1)
    }
  }else{
    P <- Y
  }


  #scaling Y

  if(transformation == "clr"){
    Mm <- apply(P, 2, mean)
    s <- apply(P, 2, sd)
    P <- (P - Mm)/s  }
  if(transformation == "ilr"){
    Mm <- mean(P)
    s <- sd(P)
    P <- (P - Mm)/s
    }


  n <- nrow(X)

  P <- as.matrix(P)
  X <- as.matrix(X)
  out <- computeWT(X = X, Y = P, A = A)

  W <- out$W
  T_score <- out$T_score

  R <- out$R
  if(post.transformation){

    if(transformation == "clr"){
      out <- ptPLSc(X = X, Y = matrix(clr(Y), ncol = ncol(Y)), W = W)
    }
    if(transformation == "ilr"){
      out <- ptPLSc(X = X, Y = matrix(ilr(Y), ncol = 1), W = W)
    }

    Wtilde <- out$Wtilde
    M <- out$M

    #apply G to weight matrix
    N <- dim(X)[1]
    E <- r <- Q <- list()

    E[[1]] <- X

    T_score <- IDA(X = X, Y = Y, W = Wtilde)
  }else{
    M <- NULL
  }

  #Compute loadings matrix
  X_loading = t(X) %*% T_score %*% solve(t(T_score) %*% T_score)
  Y_loading = t(Y) %*% T_score %*% solve(t(T_score) %*% T_score)

  #matrix coefficients
  Wstar = W %*% solve(t(X_loading) %*%W)
  B = Wstar %*% t(Y_loading)


  Y_fitted <- fitY(X = X, B = B, Mm = Mm, s = s, transformation = transformation)



  return(list(X_loading = X_loading,
              Y_loading = Y_loading,
              X = X,
              Y = Y,
              B = B,
              M = M,
              T_score = T_score,
              Y_fitted = Y_fitted))
}
