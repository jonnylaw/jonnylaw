#' Forward step for HMM
#' 
#' A single forward step of the forward filtering algorithm for a hidden Markov model
#'
#' @param alpha The previous posterior distribution over the discrete states of the HMM
#' @param x0 the initial state distribution 
#' @param transition_matrix the transition matrix
#' @param observation the observation distribution
#'
#' @return
#'
#' @examples
hmm_forward_step <- function(alpha, y, transition_matrix, observation) {
  normalise(observation(y) * t(transition_matrix) %*% alpha)
}

#' Forward Filtering an HMM
#'
#' @param ys a sequence of observations
#' @param x0 the initial state distribution 
#' @param transition_matrix the transition matrix
#' @param observation the observation distribution
#'
#' @return
#' @export
#'
#' @examples
hmm_forward <- function(ys, x0, transition_matrix, observation) {
  alpha <- normalise(observation(ys[1]) * x0)
  purrr::accumulate(
    ys[-1],
    hmm_forward_step,
    observation = observation,
    transition_matrix = transition_matrix,
    .init = alpha
  )
}

#' Title
#'
#' @param state 
#' @param y 
#' @param transition_matrix 
#' @param observation 
#'
#' @return
#'
#' @examples
ll_step <- function(state, y, transition_matrix, observation) {
  unnorm_state <- observation(y) * t(transition_matrix) %*% state[[2]]
  list(
    state[[1]] + sum(log(unnorm_state)),
    normalise(unnorm_state)
  )
}

#' Marginal log-likelihood HMM
#' 
#' Calculate the marginal log-likelihood for an HMM using the forward filtering algorithm
#'
#' @param ys a sequence of observations
#' @param x0 the initial state distribution 
#' @param transition_matrix the transition matrix
#' @param observation the observation distribution
#'
#' @return
#' @export
#'
#' @examples
hmm_log_likelihood <- function(ys, x0, transition_matrix, observation) {
  alpha <- normalise(observation(ys[1]) * x0)
  init <- list(0, alpha)
  purrr::reduce(ys, function(x, y) ll_step(x, y, observation, transition_matrix), .init = init)[[1]]
}

#' HMM Backward Smoothing
#'
#' @param beta 
#' @param y 
#' @param transition_matrix the transition matrix
#' @param observation the observation distribution
#'
#' @return
#'
#' @examples
hmm_backward_step <- function(beta, y, transition_matrix, observation) {
  normalise(transition_matrix %*% (observation(y) * beta))
}

#' Backward smoothing for the HMM
#'
#' @param ys 
#' @param transition_matrix 
#' @param observation 
#'
#' @return
#' @export
#'
#' @examples
hmm_backward <- function(ys, transition_matrix, observation) {
  purrr::accumulate(
    rev(ys),
    hmm_backward_step,
    observation = observation, 
    transition_matrix = transition_matrix,
    .init = c(1, 1)
  ) %>% rev()
}

#' Perform smoothing using the forward-backward algorithm
#'
#' @param ys The observed values
#' @param x0 the distribution of the initial state
#' @param transition_matrix the transition matrix
#' @param observation the observation distribution
#'
#' @return
#' @export
#'
#' @examples
hmm_forward_backward <- function(ys, x0, transition_matrix, observation) {
  n <- length(ys)
  f <- hmm_forward(ys, x0, transition_matrix, observation)
  b <- hmm_backward(ys, transition_matrix, observation)
  purrr::map2(f, b[-(n+1)], ~ normalise(.x * .y))
}

#' Normalise such that the values of the vector sum to one
#'
#' @param a a vector of numbers
#'
#' @return
#' @export
#'
#' @examples
normalise <- function(a) {
  a / sum(a)
}

#' Simulate from a hidden Markov model
#'
#' @param n the number of observations
#' @param x0 the initial distribution
#' @param transition_matrix the transition matrix
#' @param sim_observation the observation distribution p(y_t|x_t)
#'
#' @return
#' @export
#'
#' @examples
hmm_simulate <- function(n, x0, transition_matrix, sim_observation) {
  n_states <- ncol(transition_matrix)
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- sample(n_states, size = 1, prob = x0)
  y[1] <- sim_observation(x[1])
  for (i in seq_len(n - 1)) {
    x[i + 1] <- sample(n_states, size = 1, prob = transition_matrix[x[i], ])
    y[i + 1] <- sim_observation(x[i + 1])
  }
  tibble::tibble(time = seq_len(n), x, y)
}
