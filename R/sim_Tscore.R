sim_Tscore <- function(Tscore, n, seed){

  if(is.vector(Tscore)){Tscore <- as.matrix(Tscore, ncol = 1)}
  #### kernel density estimation (KDE) ####
  kd <- ks::kde(x = Tscore)
  #for each eval points it estimated density having as information Tscore
  #so
  ####Instrumental density normal
  mean_kd <- colMeans(kd$x)
  cov_kd <- stats::cov(kd$x)
  eval.points <- as.matrix(expand.grid(kd$eval.points))
  density_norm = mvtnorm::dmvnorm(eval.points, mean = mean_kd, sigma = cov_kd)
  kd_est <- as.vector(kd$estimate)
  ########ACCEPT-REJECTION method

  ## constant of accept reject method

  const <- max(kd_est / density_norm)

  set.seed(seed)

  n.accepts     <- 0
  Tscore_sim <- matrix(NA, nrow = n, ncol = dim(Tscore)[2])
  m <- ncol(Tscore)
  while (n.accepts < m) {
    y_sim <- mvtnorm::rmvnorm(n = 1, mean = mean_kd, sigma = cov_kd)
    f <- ks::kde(x = kd$x, eval.points = y_sim)$estimate
    g <- mvtnorm::dmvnorm(y_sim, mean = mean_kd, sigma = cov_kd)
    u <- stats::runif(1,0,1)
    if (u < f/(const*g)) {
      n.accepts <- n.accepts+1
      Tscore_sim[n.accepts,] = y_sim
    }
  }

  return(Tscore_sim)

}
