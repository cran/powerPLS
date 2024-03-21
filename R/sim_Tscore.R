#' @title simulate the T scores matrix
#' @description simulate the T scores matrix from multivariate continuous
#' distribution by using kernel density estimation and the accept-reject method.
#' @usage sim_Tscore(Tscore, n, seed)
#' @param Tscore score matrix
#' @param n number of observations
#' @param seed fix seed
#' @author Angela Andreella
#' @return Returns a simulated T score matrix.
#' @importFrom ks kde
#' @importFrom mvtnorm dmvnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cov
#' @importFrom stats runif
#' @keywords internal


sim_Tscore <- function(Tscore, n, seed){

  if(is.vector(Tscore)){Tscore <- as.matrix(Tscore, ncol = 1)}
  #### kernel density estimation (KDE) ####
  kd <- kde(x = Tscore)
  #for each eval points it estimated density having as information Tscore
  #so
  ####Instrumental density normal
  mean_kd <- colMeans(kd$x)
  cov_kd <- cov(kd$x)
  eval.points <- as.matrix(expand.grid(kd$eval.points))
  density_norm = dmvnorm(eval.points, mean = mean_kd, sigma = cov_kd)
  kd_est <- as.vector(kd$estimate)
  ########ACCEPT-REJECTION method

  ## constant of accept reject method

  const <- max(kd_est / density_norm)

  set.seed(seed)

  n.accepts     <- 0
  Tscore_sim <- matrix(NA, nrow = n, ncol = dim(Tscore)[2])
  m <- ncol(Tscore)
  while (n.accepts < m) {
    y_sim <- rmvnorm(n = 1, mean = mean_kd, sigma = cov_kd)
    f <- kde(x = kd$x, eval.points = y_sim)$estimate
    g <- dmvnorm(y_sim, mean = mean_kd, sigma = cov_kd)
    u <- runif(1,0,1)
    if (u < f/(const*g)) {
      n.accepts <- n.accepts+1
      Tscore_sim[n.accepts,] = y_sim
    }
  }

  return(Tscore_sim)

}
