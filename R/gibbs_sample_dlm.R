#' Title
#'
#' @param ys 
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
ffbs <- function(ys, model) {
  f <- model$FF; g <- model$GG; v <- model$V; w <- model$W
  m0 <- model$m0; c0 <- model$C0
  filtered <- kalmanFilter(ys, f, g, v, w, m0, c0)
  backwardSample(model$GG, model$W, filtered$mt, filtered$ct, filtered$at, filtered$rt)
}

#' Title
#'
#' @param theta 
#' @param mod 
#' @param shape_w 
#' @param rate_w 
#'
#' @return
#' @export
#'
#' @examples
sample_w <- function(theta, mod, shape_w, rate_w) {
  n <- ncol(theta) - 1
  d <- nrow(theta)
  theta_center <- t(theta[, -1, drop = FALSE]) - crossprod(theta[, -(n + 1), drop = FALSE], model$GG)
  SStheta <- drop(sapply(1:d, function(i) crossprod(theta_center[i, ])))
  
  1 / rgamma(d, shape = shape_w + 0.5 * n,
             rate = rate_w + 0.5 * SStheta)
}


#' Title
#' 
#' Centre the observations
#'
#' @param y 
#' @param theta 
#' @param mod 
#'
#' @return
#' @export
#'
#' @examples
center_y <- function(ys, theta, mod) {
  t(ys) - crossprod(theta[, -1, drop = FALSE], t(mod$FF))
}

#' Title
#' 
#' Sum the non-missing observations in each column of the observation matrix
#'
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
non_missing <- function(y) {
  rowSums(!is.na(y))
}

#' Title
#'
#' @param y 
#' @param theta 
#' @param mod 
#' @param shape_v 
#' @param rate_v 
#'
#' @return
#' @export
#'
#' @examples
sample_v <- function(y, theta, mod, shape_v, rate_v) {
  p <- nrow(y)
  ns <- non_missing(y)
  y_center <- center_y(y, theta, mod)
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
  v[1,] <- 1 / rgamma(p, shape = shape_v, rate = rate_v)
  w[1,] <- 1 / rgamma(d, shape = shape_w, rate = rate_w)
  theta <- c(w[1,], v[1,])
  if (progress) {
    pb <- progress::progress_bar$new(
      format = "sampling [:bar] :current / :total eta: :eta",
      total = iters, clear = FALSE, width= 60)  
  }
  for (i in seq_len(iters)) {
    sampled <- ffbs(ys, model)
    v[i, ] <- sample_v(ys, sampled, mod, shape_v, rate_v)
    mod$V <- v[i, ]
    w[i, ] <- sample_w(sampled, mod, shape_w, rate_w)
    mod$W <- w[i, ]
    if (progress) pb$tick()
  }
  out <- as_tibble(cbind(w, v))
  colnames(out) <- c(paste0("w", seq_len(d)), paste0("v", seq_len(p)))
  add_column(out, iteration = seq_len(iters))
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
  parallel::mclapply(seq_len(chains), function(x)
    dlm_gibbs(ys, mod, shape_v, rate_v, shape_w, rate_w, iters)) %>%
    dplyr::bind_rows(.id = "chain")
}