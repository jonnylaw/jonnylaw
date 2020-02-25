#' Title
#'
#' @param A 
#'
#' @return
#' @export
#'
#' @examples
diag_inverse <- function(A) {
  diag(1 / diag(A))
}

#' Title
#'
#' @param mu 
#' @param prec 
#'
#' @return
#' @export
#'
#' @examples
mvrnorm_prec <- function(mu, prec) {
  p <- length(mu)
  z <- rnorm(p)
  root <- chol(prec)
  mu + solve(root, z)
}

#' Sample beta
#'
#' @param ys p x n matrix of observation
#' @param xs p x m matrix of time-independent covariates
#' @param theta latent state
#' @param mod dlm model
#' @param beta_mu prior mean of the coefficients
#' @param beta_sigma prior standard deviation of the coefficients
#'
#' @return
#' @export
#'
#' @examples
sample_beta <- function(ys, xs, theta, mod, beta_mu, beta_sigma) {
  n <- ncol(ys) # number of observations
  p <- nrow(ys) # number of observed locations
  m <- nrow(xs) # dimension of the covariates
  y_center <- center_y(ys, theta, mod) %>% t()
  beta_Sigma <- beta_sigma * diag(m) # create a matrix
  beta_Mu <- rep(beta_mu, m) # Create a vector
  
  posterior_precision <- beta_Sigma + n * xs %*% diag_inverse(mod$V) %*% t(xs)
  posterior_mean <- solve(posterior_precision, beta_Sigma %*% beta_Mu + n * xs %*% diag_inverse(mod$V) %*% rowMeans(y_center))
  
  mvrnorm_prec(posterior_mean, posterior_precision)
}

#' Title
#'
#' @param ys 
#' @param xs 
#' @param mod 
#' @param shape_v 
#' @param rate_v 
#' @param shape_w 
#' @param rate_w 
#' @param beta_mu 
#' @param beta_sigma 
#' @param iters 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples
spatial_dlm <- function(ys, xs, mod, shape_v, rate_v, shape_w, rate_w, beta_mu, beta_sigma, iters, progress = TRUE) {
  p <- nrow(mod$FF)
  d <- nrow(mod$GG)
  m <- nrow(xs)
  v <- matrix(NA_real_, nrow = iters, ncol = p)
  w <- matrix(NA_real_, nrow = iters, ncol = d)
  beta <- matrix(NA_real_, nrow = iters, ncol = m)
  beta[1, ] <- rnorm(m, beta_mu, beta_sigma)
  v[1, ] <- 1 / rgamma(p, shape = shape_v, rate = rate_v)
  w[1, ] <- 1 / rgamma(d, shape = shape_w, rate = rate_w)
  theta <- c(w[1,], v[1,])
  if (progress) {
    pb <- progress::progress_bar$new(
      format = "sampling [:bar] :current / :total eta: :eta",
      total = iters, clear = FALSE, width= 60)  
  }
  for (i in seq_len(iters)[-1]) {
    ys_beta <- apply(ys, 2, function(y) y - t(xs) %*% beta[i-1, ]) # take beta away
    sampled <- ffbs(ys_beta, model) # Perform ffbs to sample theta
    w[i, ] <- sample_w(sampled, mod, shape_w, rate_w)
    mod$W <- diag(w[i, ])
    v[i, ] <- sample_v(ys_beta, sampled, mod, shape_v, rate_v)
    mod$V <- diag(v[i, ])
    beta[i, ] <- sample_beta(ys, xs, sampled, mod, beta_mu, beta_sigma)
    if (progress) pb$tick()
  }
  out <- as_tibble(cbind(w, v, beta))
  colnames(out) <- c(paste0("w", seq_len(d)), paste0("v", seq_len(p)), paste0("beta", seq_len(m)))
  add_column(out, iteration = seq_len(iters))
}

#' Gibbs sample spatial DLM in parallel
#' 
#' Run 4 chains in parallel for DLM with diagonal state and processes noise 
#' covariance matrices
#'
#' @param ys 
#' @param mod 
#' @param shape_v 
#' @param rate_v 
#' @param shape_w 
#' @param rate_w 
#' @param iters 
#' @param chains 
#'
#' @return
#' @export
#'
#' @examples
spatial_dlm_parallel <- function(ys, xs, mod, shape_v, rate_v, shape_w, rate_w, beta_mu, beta_sigma, iters = 1e4, chains = 4) {
  parallel::mclapply(seq_len(chains), function(x)
    spatial_dlm(ys, xs, mod, shape_v, rate_v, shape_w, rate_w, beta_mu, beta_sigma, iters, progress = FALSE)) %>% 
    dplyr::bind_rows(.id = "chain")
}