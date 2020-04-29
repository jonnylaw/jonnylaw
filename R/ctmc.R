#' Simulate from a continuous time markov chain
#'
#' @param x an integer, the initial state
#' @param Q a transition rate matrix
#' @param n the number of transitions to simulates
#'
#' @return
#' @export
#'
#' @examples
ctmc_sim_exact <- function(x, Q, n) {
  xs = numeric(n + 1)
  ts = numeric(n)
  r = nrow(Q)
  t = 0
  ts[1] <- t
  xs[1] = x
  for (i in seq_len(n)) {
    t <- t + rexp(1, -Q[x, x]) # Sim time to next observation
    weights <- Q[x, ] # Calculate the probability of transitioning away to another state
    weights[x] <- 0 # We can't stay in the same state
    x <- sample(r, 1, prob = weights) # Sample the next state
    xs[i + 1] <- x # add to vector of states
    ts[i + 1] <- t # add to vector of event times
    if (Q[x, x] == 0) { # If the new state is an absorbing state then return event times
      return(tibble(time = ts[seq_len(i + 1)], state = xs[seq_len(i + 1)]))
    }
  }
  tibble(time = ts, state = xs) # Return event times without absorbing
}

#' Title
#'
#' @param lambda 
#' @param n_states 
#' @param n_absorbing 
#'
#' @return
#' @export
#'
#' @examples
ctmc_build_rate_matrix <- function(lambda, n_states = 3, n_absorbing = 1) {
  q <- matrix(rep(0, times = n_states * n_states), byrow = TRUE, nrow = n_states)
  k <- 1
  for (i in 1:n_states) {
    for (j in 1:n_states) {
      if (i != j & i <= (n_states - n_absorbing)) {
        q[i, j] = lambda[k]
        k = k + 1
      }
    }
  }
  
  ## fix the diagonal to be negative the sum of the remaining elements in the row
  diag(q) = -rowSums(q)
  q
}

#' Title
#'
#' @param rate_matrix transition rate-matrix
#' @param tau a vector containing total time in each state
#' @param observed_transitions a matrix with all observed transitions
#'
#' @return
#' @export
#'
#' @examples
ctmc_log_likelihood <- function(rate_matrix, tau, observed_transitions) {
  n_states <- nrow(rate_matrix)
  qi <- - diag(rate_matrix)
  sum(log(rate_matrix - diag(diag(rate_matrix))) * observed_transitions, na.rm = TRUE) - (n_states - 1) * sum(qi * tau)
}