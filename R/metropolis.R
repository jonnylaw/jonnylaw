#' Metropolis Step
#'
#' This implements a single step of the Metropolis Algorithm
#'
#' @param theta
#' @param log_posterior
#' @param proposal
#'
#' @return
#' @export
#'
#' @examples
metropolis_step <- function(theta, log_posterior, proposal) {
  propTheta <- proposal(theta)
  a <- log_posterior(propTheta) - log_posterior(theta)
  u <- runif(1)
  if (log(u) < a) {
    propTheta
  } else {
    theta
  }
}

#' Metropolis Helper
#'
#' This implements the metropolis algorithm given an initial parameter vector,
#' the unnormalised log-posterior and a symmetric proposal distribution
#'
#' @param theta
#' @param log_posterior
#' @param proposal
#' @param m
#'
#' @return
#'
#' @examples
metropolis_helper <- function(theta, log_posterior, proposal, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(theta))
  out[1, ] <- theta
  for (i in 2:m) {
    out[i, ] <- metropolis_step(out[i - 1, ], log_posterior, proposal)
  }
  colnames(out) <- names(theta)
  tibble::as_tibble(out)
}

#' Metropolis Algorithm
#'
#' @param theta
#' @param log_posterior
#' @param proposal
#' @param m
#' @param chains
#' @param parallel
#'
#' @return
#' @export
#'
#' @examples
metropolis <- function(theta, log_posterior, proposal, m, chains = 1, parallel = FALSE) {
  if (parallel) {
    future::plan(future::multiprocess)
    mh_samples <- furrr::future_map_dfr(
      .x = 1:chains,
      .f = function(x) metropolis_helper(theta, log_posterior, proposal, m),
      .id = "chain"
    )
  } else {
    purrr::map_dfr(
      .x = 1:chains,
      .f = function(x) metropolis_helper(theta, log_posterior, proposal, m),
      .id = "chain"
    )
  }
}

#' Title
#'
#' @param theta 
#' @param ll the pseudo-marginal log-likelihood from the previous iteration
#' @param log_posterior 
#' @param proposal 
#'
#' @return
#' @export
#'
#' @examples
pmmh_step <- function(theta, ll, log_posterior, proposal) {
  propTheta <- proposal(theta)
  newll <- log_posterior(propTheta) 
  a <- newll$ll - ll
  u <- runif(1)
  if (log(u) < a) {
    c(newll$ll, propTheta)
  } else {
    c(ll, theta)
  }
}

#' Metropolis Helper
#'
#' This implements the metropolis algorithm given an initial parameter vector,
#' the unnormalised log-posterior and a symmetric proposal distribution
#'
#' @param theta
#' @param log_posterior
#' @param proposal
#' @param m
#'
#' @return
#'
#' @examples
pmmh_helper <- function(theta, log_posterior, proposal, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(theta) + 1)
  out[1, ] <- c(-1e99, theta) # Initialise loglikelihood and state
  for (i in 2:m) {
    out[i, ] <- pmmh_step(out[i - 1, -1], out[i - 1, 1], log_posterior, proposal)
  }
  colnames(out) <- c("log_likelihood", names(theta))
  tibble::as_tibble(out)
}

#' Pseudo-Marginal Metropolis Hastings
#'
#' @param theta 
#' @param log_posterior 
#' @param proposal 
#' @param m 
#' @param chains 
#'
#' @return
#' @export
#'
#' @examples
pmmh <- function(theta, log_posterior, proposal, m, chains = 2) {
  furrr::future_map_dfr(
    .x = seq_len(chains),
    .f = function(x) jonnylaw:::pmmh_helper(theta, log_posterior, proposal, m),
    .id = "chain"
  )
}
