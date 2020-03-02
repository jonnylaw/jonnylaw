#' Leapfrog Step
#'
#' @param gradient
#' @param step_size
#' @param position
#' @param momentum
#'
#' @return
#' @export
#'
#' @examples
leapfrog_step <- function(gradient, step_size, position, momentum, d) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum1 + gradient(position1) * 0.5 * step_size

  matrix(c(position1, momentum2), ncol = d * 2)
}

#' Perform l Leapfrog steps
#'
#' @param gradient
#' @param step_size
#' @param l
#' @param position
#' @param momentum
#'
#' @return
#' @export
#'
#' @examples
leapfrogs <- function(gradient, step_size, l, position, momentum, d) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum, d)
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
#' @export
#'
#' @examples
log_acceptance <- function(prop_position,
                           prop_momentum,
                           position,
                           momentum,
                           log_posterior) {
  a <- log_posterior(prop_position) + sum(dnorm(prop_momentum, log = T)) -
    log_posterior(position) - sum(dnorm(momentum, log = T))
  if_else(is.na(a), -Inf, min(a, 0.0))
}

#' HMC Step
#'
#' @param log_posterior
#' @param gradient
#' @param step_size
#' @param l
#' @param position
#'
#' @return
#' @export
#'
#' @examples
hmc_step <- function(log_posterior, gradient, step_size, l, position, d) {
  momentum <- rnorm(d)
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum, d)
  prop_position <- pos_mom[seq_len(d)]
  prop_momentum <- pos_mom[-seq_len(d)]
  a <- log_acceptance(prop_position, prop_momentum, position, momentum, log_posterior)
  if (log(runif(1)) < a) {
    prop_position
  } else {
    position
  }
}

#' Hamiltonian Monte Carlo
#'
#' @param log_posterior
#' @param gradient
#' @param step_size
#' @param l
#' @param theta_0
#' @param m
#'
#' @return
#' @export
#'
#' @examples
hmc <- function(log_posterior, gradient, step_size, l, theta_0, m) {
  d <- length(theta_0)
  out <- matrix(NA_real_, nrow = m, ncol = d)
  out[1, ] <- theta_0
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i - 1, ], d)
  }
  out
}

#' Hamiltonian Monte Carlo
#'
#' Perform HMC and return a dataframe
#'
#' @param log_posterior
#' @param gradient
#' @param step_size
#' @param l
#' @param theta_0
#' @param m
#' @param parameter_names
#'
#' @return
#' @export
#'
#' @examples
hmc_df <- function(log_posterior, gradient, step_size, l, theta_0, m, parameter_names) {
  mat <- hmc(log_posterior, gradient, step_size, l, theta_0, m)
  colnames(mat) <- parameter_names
  as.data.frame(mat) %>%
    mutate(iteration = row_number())
}
