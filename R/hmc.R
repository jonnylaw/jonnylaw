#' Leapfrog Step
#'
#' @param gradient
#' @param step_size
#' @param position
#' @param momentum
#'
#' @return
#'
#' @examples
hmc_leapfrog_step <- function(gradient, step_size, position, momentum, d) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum1 + gradient(position1) * 0.5 * step_size

  matrix(c(position1, momentum2), ncol = d * 2)
}

#' Perform n_steps Leapfrog steps
#'
#' @param gradient the gradient of the target
#' @param step_size the step size of the leapfrog step
#' @param n_steps the total number of steps
#' @param position the current value of the parameter vector, termed position in Hamiltonian
#' dynamics
#' @param momentum the current "momentum" in the HMC sampler
#' @param d dimension of parmaeter vector
#'
#' @return
#'
#' @examples
hmc_leapfrogs <- function(gradient, step_size, n_steps, position, momentum, d) {
  for (i in seq_len(n_steps)) {
    pos_mom <- hmc_leapfrog_step(gradient, step_size, position, momentum, d)
    position <- pos_mom[seq_len(d)]
    momentum <- pos_mom[-seq_len(d)]
  }
  pos_mom
}

#' Log Acceptance for the HMC algorithm
#'
#' @param prop_position
#' @param prop_momentum
#' @param position
#' @param momentum
#' @param log_posterior
#'
#' @return
#'
#' @examples
hmc_log_acceptance <- function(prop_position,
                               prop_momentum,
                               position,
                               momentum,
                               log_posterior) {
  a <- log_posterior(prop_position) + sum(dnorm(prop_momentum, log = T)) -
    log_posterior(position) - sum(dnorm(momentum, log = T))
  dplyr::if_else(is.na(a), -Inf, min(a, 0.0))
}

#' HMC Step
#'
#' @param log_posterior a function for the unnormalised log-posterior
#' from parameters -> log-likelihood
#' @param gradient the gradient of the log_posterior
#' @param step_size the step size of the leapfrog steps
#' @param n_steps the number of leapfrog steps
#' @param position the current value of the parameter vector
#' @param d the dimension of the parameter vector
#'
#' @return
#'
#' @examples
hmc_step <- function(log_posterior, gradient, step_size, n_steps, position, d) {
  momentum <- rnorm(d) # Sample initial momentum from standard normal 
  pos_mom <- hmc_leapfrogs(gradient, step_size, n_steps, position, momentum, d)
  prop_position <- pos_mom[seq_len(d)]
  prop_momentum <- pos_mom[-seq_len(d)]
  a <- hmc_log_acceptance(prop_position, prop_momentum, position, momentum, log_posterior)
  if (log(runif(1)) < a) {
    prop_position
  } else {
    position
  }
}

#' Hamiltonian Monte Carlo
#'
#' Perform Hamiltonian Monte Carlo
#'
#' @param log_posterior a function for the unnormalised log-posterior
#' from parameters -> log-likelihood
#' @param gradient the gradient of the log_posterior
#' @param step_size the step size of the leapfrog steps
#' @param n_steps the number of leapfrog steps
#' @param init_parameters a named vector containing the initial parameters
#' @param iters the number of iterations
#'
#' @return a matrix of iterations from the HMC algorithm representing draws
#' from the posterior distribution
#'
#' @examples
hmc_helper <- function(log_posterior, gradient, step_size, n_steps, init_parameters, iters) {
  d <- length(init_parameters)
  out <- matrix(NA_real_, nrow = iters, ncol = d)
  out[1, ] <- init_parameters
  for (i in 2:iters) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, n_steps, out[i - 1, ], d)
  }
  colnames(out) <- names(init_parameters)
  tibble::as_tibble(out) %>%
    dplyr::mutate(iteration = dplyr::row_number())
}

#' Hamiltonian Monte Carlo
#'
#' Run multiple HMC chains and return a dataframe of parameters
#' 
#' To run this in parallel call
#' `future::plan(future::multiprocess())` before this function
#'
#' @param log_posterior a function for the unnormalised log-posterior
#' from parameters -> log-likelihood
#' @param gradient the gradient of the log_posterior
#' @param step_size the step size of the leapfrog steps
#' @param n_steps the number of leapfrog steps
#' @param init_parameters a named vector of initial parameters
#' @param iters the number of iterations 
#' @param chains the number of chains
#'
#' @return
#' @export
#'
#' @examples
hmc <- function(log_posterior, gradient, step_size, n_steps, init_parameters, iters, chains = 2) {
  furrr::future_map_dfr(
    .x = seq_len(chains),
    .f = function(x)
      jonnylaw:::hmc_helper(
        log_posterior,
        gradient,
        step_size,
        n_steps,
        init_parameters,
        iters
      ),
    .id = "chain"
  )
}
