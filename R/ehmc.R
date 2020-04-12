#' Title
#'
#' Find the longest batch for a singe position
#'
#' @param position 
#' @param momentum 
#' @param step_size 
#' @param n_steps 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
longest_batch <- function(
  gradient,
  position,
  momentum,
  step_size,
  n_steps,
  d) {
  
  ## Initialise values
  l <- 1
  prop_position <- position
  prop_position <- momentum
  final_position <- position
  final_position <- momentum
  
  while (sum((prop_position - position) * prop_momentum) >= 0 & l < 2000) { ## Horrible dot-product
    pos_mom <-
      hmc_leapfrogs(gradient, step_size, n_steps, prop_position, prop_momentum, d)
    prop_position <- pos_mom[seq_len(d)]
    prop_momentum <- pos_mom[-seq_len(d)]
    if (l == n_steps) {
      final_position <- prop_position
      final_momentum <- prop_momentum
    }
    l = l + 1
  }
   
  ## Return a vector containing the number of leapfrog steps and 
  ## the position and momentum at the end of the steps
  if (l > n_steps) {
    c(l, final_position, final_momentum)      
  } else {
    c(l, prop_position, prop_momentum)
  }
}

longest_batch_step <- function(
  gradient,
  log_posterior,
  position, 
  step_size,
  l0,
  d) {

  ## Sample momentum from standard normal
  momentum <- rnorm(d)
  
  lpos_mom <- longest_batch(gradient, position, momentum, step_size, l0, d)
  l <- lpos_mom[1]
  if (l < l0) {
    pos_mom <- hmc_leapfrogs(gradient, step_size, l0 - l, position, momentum, d)
  } else {
    pos_mom <- lpos_mom[-1]
  }
  
  a = hmc_log_acceptance(prop_position,
                         prop_momentum,
                         position,
                         momentum,
                         log_posterior)
  if (log(runif(1)) < a) {
    c(l, prop_position)
  } else {
    c(l, position)
  }
}

## Determine the empirical distribution of step-sizes
empirical_longest_batch <- function(
  gradient,
  log_posterior,
  position,
  step_size,
  l0,
  k = 200) {

  ## Initialise values
  d <- length(position)
  ls <- numeric(k)
  
  for (i in seq_len(k)) {
    lpos <- longest_batch_step(gradient, log_posterior, position, step_size, l0, d)
    position <- lpos[-1]
    ls[i] <- lpos[1]
  }
  
  ls
}

#' Title
#'
#' @param log_posterior 
#' @param gradient 
#' @param initial_step_size 
#' @param l0 
#' @param iteration 
#' @param warmup_iterations 
#' @param position 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
ehmc_step <- function(
  log_posterior,
  gradient,
  position,
  initial_step_size,
  l0 = 10,
  iteration,
  warmup_iterations,
  d) {
  
  ## sample the momentum from a standard normal
  momentum = rnorm(d)
  
  ## Set the step-size
  step_size <- dplyr::if_else(iteration < warmup_iterations, exp(logeps), exp(logepsbar))
  
  ## draw the number of leapfrog steps from the empirical distribution of step sizes
  l <- sample(empirical_n_steps, size = 1)
  
  ## perform l leapfrog ststep_size
  pos_mom = hmc_leapfrogs(gradient, step_size, l, position, momentum, d)
  prop_position <- pos_mom[seq_len(d)]
  prop_momentum <- pos_mom[-seq_len(d)]
  
  ## calculate acceptance
  a <- hmc_log_acceptance(prop_position, prop_momentum, position, momentum, log_posterior)
  
  ## Update the step size in the adaptation phase
  if (iteration < warmup_iterations) {
    update_step_size(iteration, min(1.0, exp(a)), log(10 * initial_step_size))
  }
  
  ## Calculate the empirical distribution of step sizes at the final warmup iteration
  if (iteration == warmup_iterations)
    empirical_n_steps <<- empirical_longest_batch(gradient, log_posterior, prop_position, step_size, l0)
  
  ## Accept new parameters with probability max(1.0, exp(a))
  if (log(runif(1)) < a) {
    ## return a vector containing a boolean denoting
    ## an accepted move and the new parameters
    c(TRUE, prop_position) 
  } else {
    c(FALSE, position) 
  }
}

#' Empirical Hamiltonian Monte Carlo
#'
#' Perform Empirical Hamiltonian Monte Carlo
#'
#' @param log_posterior a function for the unnormalised log-posterior
#' from parameters -> log-likelihood
#' @param gradient the gradient of the log_posterior
#' @param n_steps the initial number of leapfrog steps
#' @param init_parameters a named vector containing the initial parameters
#' @param iters the number of iterations
#'
#' @return a matrix of iterations from the HMC algorithm representing draws
#' from the posterior distribution
#'
#' @examples
ehmc_helper <- function(
  log_posterior, 
  gradient, 
  n_steps, 
  init_parameters, 
  iters, 
  warmup_iters = iters / 2, 
  logging = TRUE) {
  
  ## Calculate the dimension of the parameter vector
  d <- length(init_parameters)
  
  ## Find a reasonable initial step size
  initial_step_size <-
    find_reasonable_epsilon(init_parameters, momentum, log_posterior, gradient, d)
  
  ## Initialise the (scary) global values for dual averaging and EHMC
  hm <<- 0
  logepsbar <<- 0
  logeps <<- log(initial_step_size)
  empirical_n_steps <<- n_steps
  
  ## Initialise a matrix to contain the sampled parameter vectors and
  ## Accepted boolean
  out <- matrix(NA_real_, nrow = iters, ncol = d + 1)
  out[1, -1] <- init_parameters
  
  ## Initialise the logged values
  if (logging) {
    hms <- numeric(iters)
    logepsbars <- numeric(iters)
    logepss <- numeric(iters)
  }
  
  for (i in 2:iters) {
    out[i, ] <- ehmc_step(log_posterior,
                          gradient,
                          out[i - 1, -1],
                          initial_step_size,
                          n_steps,
                          i, 
                          warmup_iters,
                          d)
    
    ## Log the values of hm, logepsbar and logeps
    if (logging) {
      hms[i] <- hm
      logepsbars[i] <- logepsbar
      logepss[i] <- logeps
    }
  }
  
  colnames(out) <- c("accepted", names(init_parameters))
  out <- tibble::as_tibble(out) %>%
    dplyr::mutate(iteration = dplyr::row_number())
  
  if (logging) {
    out <- out %>%
      tibble::add_column(hm = hms,
                         log_step_size = logepss,
                         log_step_size_bar = logepsbars)
    list(iters = out, ls = empirical_n_steps)
  } else {
    out
  }
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
ehmc <- function(log_posterior, gradient, n_steps, init_parameters, iters, chains = 2) {
  furrr::future_map_dfr(
    .x = seq_len(chains),
    .f = function(x)
      jonnylaw:::ehmc_helper(
        log_posterior,
        gradient,
        n_steps,
        init_parameters,
        iters
      ),
    .id = "chain"
  )
}