#' @title Simulate pilot data
#' @description Simulate data matrix under the alternative hypothesis with \code{n} observations by kernel density estimation
#' @usage sim_XY(out, n, seed = 123, post.transformation = TRUE, A, fast = FALSE)
#' @param out Output from \code{PLSc}
#' @param n Number of observations to simulate
#' @param seed Seed value
#' @param post.transformation Boolean value. Default to \code{TRUE}, i.e., post transformation is applied in \code{PLSc}
#' @param A Number of score components used in \code{PLSc}.
#' @param fast Use the function \code{fk_density} from the \code{FKSUM} \code{R} package for kernel density estimation. Default to \code{FALSE}.
#' @author Angela Andreella
#' @return Returns a list:
#' \describe{
#'   \item{Y_H1}{dependent variable, matrix with 2 columns and \code{n} rows (observations)}
#'   \item{X_H1}{predictor variables, matrix with \code{n} rows (observations) and number of columns equal to \code{out$X} (i.e., original dataset)}
#' }
#' @export
#' @seealso \code{\link{PLSc}}, \code{\link{ptPLSc}}
#' @importFrom FKSUM fk_density
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom simukde simulate_kde
#' @examples
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- PLSc(X = datas$X, Y = datas$Y, A = 3)
#' out_sim <- sim_XY(out = out, n = 10, A = 3)
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' Andreella, A., Fino, L., Scarpa, B., & Stocchero, M. (2024). Towards a power analysis for PLS-based methods. arXiv preprint \url{https://arxiv.org/abs/2403.10289}.

sim_XY <- function(out, n, seed = 123, post.transformation = TRUE, A, fast = FALSE){


  set.seed(seed)

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
   # sim_TO <- sapply(seq(A), function(x) {
   #   out_kde <- fk_density(x = T_scoreO[,x])
    #  sample(out_kde$x, size = n, prob = out_kde$y)}
   # )
    sim_TO <- sapply(seq(ncol(T_scoreO)),function(x){
      if(fast){
      out_kde <- fk_density(x = T_scoreO[,x])
      sample(out_kde$x, size = n, prob = out_kde$y)
     # kde_transf <- ks::kde(x = T_scoreO[,x])
     # ks::rkde(n = n, fhat = kde_transf)
        }
      else{
      simulate_kde(x = T_scoreO[,x], n = n)$random.values
      }
    })

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

    sim_TP <- sapply(seq(ncol(T_scoreP)),function(x){
      if(fast){
      out_kde <- fk_density(x = T_scoreP[,x])
      sample(out_kde$x, size = n, prob = out_kde$y)
     # kde_transf <- ks::kde(x = T_scoreP[,x])
     # ks::rkde(n = n, fhat = kde_transf)
      }else{
        simulate_kde(x = T_scoreP[,x], n = n)$random.values
      }
      })

    sim_TO <- sapply(seq(ncol(T_scoreO)),function(x){
      if(fast){
      out_kde <- fk_density(x = T_scoreO[,x])
      sample(out_kde$x, size = n, prob = out_kde$y)
    #  kde_transf <- ks::kde(x = T_scoreO[,x])
    #  ks::rkde(n = n, fhat = kde_transf)
      }else{
        simulate_kde(x = T_scoreO[,x], n = n)$random.values
      }
      })
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
