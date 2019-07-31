#' Title
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
leapfrog_step <- function(gradient, step_size, position, momentum) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum + gradient(position1) * 0.5 * step_size

  list(position = position1, momentum = momentum2)
}

#' Title
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
leapfrogs <- function(gradient, step_size, l, position, momentum) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum)
  }
  list(position = pos_mom[["position"]], momentum = pos_mom[["momentum"]])
}

#' Title
#'
#' @param propPosition
#' @param propMomentum
#' @param position
#' @param momentum
#' @param log_posterior
#'
#' @return
#' @export
#'
#' @examples
log_acceptance <- function(propPosition,
                           propMomentum,
                           position,
                           momentum,
                           log_posterior) {
  log_posterior(propPosition) + sum(dnorm(propMomentum, log = T)) -
    log_posterior(position) - sum(dnorm(momentum, log = T))
}

#' Title
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
hmc_step <- function(log_posterior, gradient, step_size, l, position) {
  d <- length(position)
  momentum <- rnorm(d)
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum)
  propPosition <- pos_mom[["position"]]
  propMomentum <- pos_mom[["momentum"]]
  a <- log_acceptance(propPosition, propMomentum, position, momentum, log_posterior)
  if (log(runif(1)) < a) {
    propPosition
  } else {
    position
  }
}

hmc <- function(log_posterior, gradient, step_size, l, initP, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(initP))
  out[1, ] <- initP
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i-1,])
  }
  out
}
