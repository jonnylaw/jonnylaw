#' Perform Forward Filtering Backward Sampling
#'
#' @param ys p x n matrix of observation
#' @param mod the dlm model
#'
#' @return
#' @export
#'
#' @examples
dlm_ffbs <- function(ys, mod) {
  f <- mod$FF
  g <- mod$GG
  v <- mod$V
  w <- mod$W
  m0 <- mod$m0
  c0 <- mod$C0
  filtered <- dlm_kalman_filter(ys, f, g, v, w, m0, c0)
  dlm_backward_sample(mod$GG, mod$W, filtered$mt, filtered$ct, filtered$at, filtered$rt)
}

#' Sample the system noise matrix
#' 
#' Sample the system noise matrix from a d-inverse Gamma distribution
#'
#' @param theta a d x n matrix containing the sampled latent state
#' @param mod the dlm model
#' @param shape_w
#' @param rate_w
#'
#' @return
#' @export
#'
#' @examples
dlm_sample_w <- function(theta, mod, shape_w, rate_w) {
  n <- ncol(theta) - 1
  d <- nrow(theta)
  theta_center <- t(theta[, -1, drop = FALSE]) - crossprod(theta[, -(n + 1), drop = FALSE], mod$GG)
  SStheta <- drop(sapply(1:d, function(i) crossprod(theta_center[i, ])))

  1 / rgamma(d,
    shape = shape_w + 0.5 * n,
    rate = rate_w + 0.5 * SStheta
  )
}

#' Helper function for Kalman Filter and Gibbs Sampling
#'
#' Centre the observations
#'
#' @param ys p x n matrix of observation
#' @param theta d x n matrix of sampled latent state
#' @param mod the dlm model
#'
#' @return
#'
#' @examples
dlm_center_y <- function(ys, theta, mod) {
  t(ys) - crossprod(theta[, -1, drop = FALSE], t(mod$FF))
}

#' Helper function for Kalman Filter
#'
#' Sum the non-missing observations in each column of the observation matrix
#'
#' @param y
#'
#' @return
#'
#' @examples
dlm_non_missing <- function(y) {
  rowSums(!is.na(y))
}

#' Sample the observation noise matrix
#' 
#' Sample the observation noise matrix from a d-Inverse Gamma distribution
#'
#' @param ys p x n matrix of observation
#' @param theta d x n matrix of the sampled latent state
#' @param mod the dlm model
#' @param shape_v the prior shape of the IG distribution
#' @param rate_v the prior rate of the IG distribution
#'
#' @return
#' @export
#'
#' @examples
dlm_sample_v <- function(y, theta, mod, shape_v, rate_v) {
  p <- nrow(y)
  ns <- dlm_non_missing(y)
  y_center <- dlm_center_y(y, theta, mod)
  SSy <- crossprod(na.omit(y_center))

  purrr::map2_dbl(ns, diag(SSy), ~ 1 / rgamma(1, shape = shape_v + 0.5 * .x, rate = rate_v + 0.5 * .y))
}

#' Gibbs sampling DLM
#'
#' Sample diagonal system and noise covariance matrices using a Gibbs sampler
#'
#' @param ys
#' @param mod
#' @param shape_v
#' @param rate_v
#' @param shape_w
#' @param rate_w
#' @param iters
#'
#' @return
#' @export
#'
#' @examples
dlm_gibbs <- function(ys, mod, shape_v, rate_v, shape_w, rate_w, iters, progress = FALSE) {
  p <- nrow(mod$FF)
  d <- nrow(mod$GG)
  v <- matrix(NA_real_, nrow = iters, ncol = p)
  w <- matrix(NA_real_, nrow = iters, ncol = d)
  v[1, ] <- 1 / rgamma(p, shape = shape_v, rate = rate_v)
  w[1, ] <- 1 / rgamma(d, shape = shape_w, rate = rate_w)
  mod$V <- diag(v[1, ])
  mod$W <- diag(w[1, ])
  if (progress) {
    pb <- progress::progress_bar$new(
      format = "sampling [:bar] :current / :total eta: :eta",
      total = iters, clear = FALSE, width = 60
    )
  }
  for (i in seq_len(iters)) {
    sampled <- dlm_ffbs(ys, mod)
    v[i, ] <- dlm_sample_v(ys, sampled, mod, shape_v, rate_v)
    mod$V <- diag(v[i, ])
    w[i, ] <- dlm_sample_w(sampled, mod, shape_w, rate_w)
    mod$W <- diag(w[i, ])
    if (progress) pb$tick()
  }
  out <- tibble::as_tibble(cbind(w, v))
  colnames(out) <- c(paste0("w", seq_len(d)), paste0("v", seq_len(p)))
  tibble::add_column(out, iteration = seq_len(iters))
}

#' Gibbs sample DLM in parallel
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
dlm_gibbs_parallel <- function(ys, mod, shape_v, rate_v, shape_w, rate_w, iters, chains = 4) {
  chains <- seq_len(chains)
  names(chains) <- chains
  furrr::future_map_dfr(chains, dlm_gibbs(ys, mod, shape_v, rate_v, shape_w, rate_w, iters), .id = "chain")
}

#' Sample the system noise matrix
#' 
#' Sample the system noise matrix from a d-inverse Gamma distribution
#'
#' @param theta a d x n matrix containing the sampled latent state
#' @param mod the dlm model
#' @param shape_w
#' @param rate_w
#'
#' @return
#' @export
#'
#' @examples
dlm_sample_w_wishart <- function(theta, mod, nu, psi) {
  n <- ncol(theta) - 1
  d <- nrow(theta)
  theta_center <- t(theta[, -1, drop = FALSE]) - crossprod(theta[, -(n + 1), drop = FALSE], mod$GG)
  SStheta <- crossprod(theta_center)
  
  MCMCpack::riwish(v = nu + n, S = psi + SStheta)
}

#' Sample the observation noise matrix
#' 
#' Sample the observation noise matrix from a d-Inverse Gamma distribution
#'
#' @param ys p x n matrix of observation
#' @param theta d x n matrix of the sampled latent state
#' @param mod the dlm model
#' @param shape_v the prior shape of the IG distribution
#' @param rate_v the prior rate of the IG distribution
#'
#' @return
#' @export
#'
#' @examples
dlm_sample_v_wishart <- function(y, theta, mod, shape_v, rate_v) {
  p <- nrow(y)
  ns <- dlm_non_missing(y)
  y_center <- dlm_center_y(y, theta, mod)
  SSy <- crossprod(y_center)
  
  MCMCpack::riwish(ncol(theta) - 1, SSy)
}

#' Title
#'
#' @param ys 
#' @param mod 
#' @param nu_v 
#' @param psi_v 
#' @param nu_w 
#' @param psi_w 
#' @param iters 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples
dlm_gibbs_wishart <- function(ys, mod, nu_v, psi_v, nu_w, psi_w, iters, progress = FALSE) {
  p <- nrow(mod$FF)
  d <- nrow(mod$GG)
  v <- array(NA_real_, dim = c(p, p, iters))
  w <- array(NA_real_, dim = c(d, d, iters))
  v[ , , 1] <- MCMCpack::riwish(nu_v, psi_v)
  w[ , , 1] <- MCMCpack::riwish(nu_w, psi_w)
  mod$V <- v[ , , 1]
  mod$W <- w[ , , 1]
  if (progress) {
    pb <- progress::progress_bar$new(
      format = "sampling [:bar] :current / :total eta: :eta",
      total = iters, clear = FALSE, width = 60
    )
  }
  for (i in seq_len(iters)) {
    sampled <- dlm_ffbs(ys, mod)
    v[ , , i] <- dlm_sample_v_wishart(ys, sampled, mod, nu_v, psi_v)
    mod$V <- v[ , , i]
    w[ , , i] <- dlm_sample_w_wishart(sampled, mod, nu_w, psi_w)
    mod$W <- w[ , , i]
    if (progress) pb$tick()
  }
  out <- tibble::as_tibble(cbind(plyr::aaply(v, 3, c), plyr::aaply(w, 3, c)))
  colnames(out) <- c(paste0("v", seq_len(p*p)), paste0("w", seq_len(d*d)))
  tibble::add_column(out, iteration = seq_len(iters))
}

#' Transform a dlm to a KFAS State Space Model
#'
#' @param dlm_model a model created using the dlm package
#' @param ys a matrix containing observations
#' @param theta a vector of parameters containing the diagonal values of v and w
#'
#' @return
#' @export
#'
#' @examples
dlm_to_kfas <- function(dlm_model, ys, theta) {
  f <- dlm_model$FF
  g <- dlm_model$GG
  w <- dlm_model$W
  v <- dlm_model$V

  SSModel(ys ~ SSMcustom(Z = f, T = g, Q = diag(w)), H = diag(v))
}



plot.filtered <- function(filtered) {
  d <- nrow(filtered$mt)
  filtered_state <-
    filtered$mt[,-1] %>% t() %>% 
    as_tibble() 
  names(filtered_state) <- paste0(rep("state", time = d), seq_len(d))
  
  filtered_state <- filtered_state %>% 
    mutate(time = row_number()) %>% 
    pivot_longer(-time, names_to = "variable", values_to = "filtered")
  
  state_df <- state %>% 
    tibble::as_tibble()
  names(state_df) <- paste0(rep("state", time = d), seq_len(d))
  
  state_df <- state_df %>% 
    mutate(time = row_number()) %>% 
    pivot_longer(-time, names_to = "variable", values_to = "state")
  
  filtered_state %>% 
    inner_join(state_df, by = c("time", "variable")) %>% 
    pivot_longer(c(state, filtered), names_to = "state", values_to = "value") %>% 
    ggplot() +
    geom_line(aes(x = time, y = value, colour = state)) +
    facet_wrap(~variable, ncol = 1)
}