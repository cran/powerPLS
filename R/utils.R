
#Utils

#Fit Y
fitY <- function(X, B, Mm, s, transformation = "clr"){

  if(transformation == "clr"){
    Y.fitted = matrix(apply(compositions::clrInv(s*(X %*% B) + Mm), 1, function(x)
      which.max(x)),ncol = 1)
    }
  if(transformation == "ilr"){
    Y.fitted = matrix(apply(compositions::ilrInv(s*(X %*% B) + Mm), 1, function(x)
      which.max(x)),ncol = 1)
    }


  Y.fitted = as.factor(Y.fitted)
  lev <- levels(Y.fitted)
  if(length(lev)>1){
    Y.fitted = model.matrix(~0+Y.fitted)
  }


  return(Y.fitted)
}

permuteIndex <- function(Y, by.row = TRUE, times, replace = FALSE){

  if(by.row){
    idx <- sample(seq(nrow(Y)), size = times, replace = replace)
    Y <- Y[idx,]
  }else{
    idx <- sample(seq(ncol(Y)), size =times, replace = replace)
    Y <- Y[,idx]
  }
return(Y)
}

similarityMatrix <- function(X, Y){

  tr <- function(X){

    return(sum(diag(X)))
  }
  n <- nrow(X)
  XX <- tcrossprod(X)
  YY <- tcrossprod(Y)
  XY <- X %*% t(Y)
  YX <- Y %*% t(X)


  #RV index (Escoufier, 1973; Robert and Escoufier, 1976)
  RV <- tr(XY %*% YX)/sqrt(tr(XX)^2 * tr(YY)^2)
  #RLS index (Gower 1971; Lingoes and Schonemann (1974))
  RLS <- sqrt(tr(crossprod(X) %*% crossprod(Y)))/ sqrt(tr(XX) %*% tr(YY))

  out <- data.frame(RV = RV, RLS = RLS)

  return(out)
}

mcc <- function(confMatrix){

  confMatrix <-stats::addmargins(confMatrix)
  n11 <- confMatrix[1,1]
  n22 <- confMatrix[2,2]
  n12 <- confMatrix[1,2]
  n21 <- confMatrix[2,1]

  out <- (n11*n22 - n12*n21)/sqrt(confMatrix[1,3]*confMatrix[2,3]*confMatrix[3,1]*confMatrix[3,2])
  if(is.nan(out)){out <- 0}
  return(out)

}
