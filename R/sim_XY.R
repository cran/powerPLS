#' @title simulate data matrix under the alternative hypothesis
#' @description simulate data matrix under the alternative hypothesis
#' @usage sim_XY(out, n, seed = 123, post.transformation = TRUE, A)
#' @param out output from \code{PLSc}
#' @param n number of observations to simulate
#' @param seed seed value
#' @param post.transformation Boolean value. Default @TRUE i.e., post transformation is applied.
#' @param A number of score components used in \code{PLSc}.
#' @author Angela Andreella
#' @return Returns a simulated matrix under the alternative hypothesis.
#' @export
#' @importFrom simukde simulate_kde
#' @examples
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- PLSc(X = datas$X, Y = datas$Y, A = 2)
#' out_sim <- sim_XY(out = out, n = 10, A = 2)

sim_XY <- function(out, n, seed = 123, post.transformation = TRUE, A){


  set.seed(seed)
#  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling, post.transformation = post.transformation)
  if(post.transformation){
    M <- out$M
  }else{
    M <- A
  }

  T_score <- out$T_score
  X_loading <- out$X_loading
  Y_loading <- out$Y_loading
  B <- out$B
  X <- out$X

  if(A == M){
    T_scoreO <- T_score
    sim_TO <- sapply(seq(A), function(x) simulate_kde(x = T_scoreO[,x], n = n)$random.values)
    T_sim <- sim_TO #T target

  }else{
    T_scoreO <- T_score[,1:M]
    if(is.null(ncol(T_score[,((M+1):A)]))){
      ncp <- 1
    }else{
      ncp <- ncol(T_score[,((M+1):A)])
    }

    if(is.null(ncol(T_score[,1:M]))){
      nco <- 1
    }else{
      nco <- ncol(T_score[,1:M])
    }
    T_scoreO <- matrix(T_score[,(1:M)], ncol = nco)
    T_scoreP <- matrix(T_score[,((M+1):A)], ncol = ncp)


    sim_TP <- unlist(sapply(seq(ncol(T_scoreP)), function(x) simulate_kde(x = T_scoreP[,x],
                                                                          n = n)$random.values))
    sim_TO <- unlist(sapply(seq(ncol(T_scoreO)), function(x) simulate_kde(x = T_scoreO[,x],
                                                                          n = n)$random.values))
  #  sim_TP <- scale(sim_TP, center = FALSE)
  #  sim_TO <- scale(sim_TO, center = FALSE)


    T_sim <- cbind(sim_TO, sim_TP) #T target

  }


  out1 <- svd(T_sim %*% t(T_score) %*% T_score)

  S <- (t(T_score) %*% T_score)
  if(length(S)==1){
    T_new <- out1$u %*% t(out1$v) %*% diag(S)^(1/2)
  }else{
    T_new <- out1$u %*% t(out1$v) %*% diag(diag(S))^(1/2)
  }


  E_pilot <- X - T_score %*% t(X_loading)

  E <- permuteIndex(E_pilot, times = n, by.row = TRUE, replace = TRUE)
  E <- as.matrix(E)

  E_X <- (diag(dim(T_new)[1]) - T_new %*% solve(t(T_new) %*% T_new) %*% t(T_new)) %*% E

  X_H1 <- T_new %*% t(X_loading) + E_X


  Y_H1 <- fitY(X = X_H1, B = B, Mm = 0, s = 1)



  return(list(Y_H1 = Y_H1, X_H1 = X_H1))
}
