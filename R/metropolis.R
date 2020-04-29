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
    c(TRUE, propTheta)
  } else {
    c(FALSE, theta)
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
  out <- matrix(NA_real_, nrow = m, ncol = length(theta) + 1)
  out[1, ] <- c(TRUE, theta)
  for (i in 2:m) {
    out[i, ] <- metropolis_step(out[i - 1, -1], log_posterior, proposal)
  }
  colnames(out) <- c("accepted", names(theta))
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
