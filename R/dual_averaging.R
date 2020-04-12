#' Title
#'
#' Update the value of epsilon during an adaptation phase
#'
#' @param iteration the current iteration
#' @param accept_prob the acceptance probability from the previous time step
#' @param mu tuning parameter, default set to log(10 * eps0)
#' @param delta target acceptance probability, default 0.65
#' @param k tuning parameter, default set to 0.75
#' @param gamma tuning parameter, default set to 0.05
#' @param t0 tuning parameter, default set to 10.0
#' @param hm0 the previous state
#' @param logeps0 the previous state
#' @param logepsbar0 the previous state
#'
#' @return
#'
#' @examples
update_step_size <- function(
                       iteration,
                       accept_prob,
                       mu,
                       delta = 0.65,
                       k = 0.75,
                       gamma = 0.05,
                       t0 = 10.0) {
  ra <- 1 / (iteration + t0)
  hm <<- (1 - ra) * hm + ra * (delta - accept_prob)
  logeps <<- mu - ((sqrt(iteration) * hm) / gamma)
  power <- iteration^(-k)
  logepsbar <<- power * logeps + (1.0 - power) * logepsbar
}

#' Title
#'
#' @param position 
#' @param momentum 
#' @param log_posterior 
#' @param gradient 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
find_reasonable_epsilon <- function(
                                    position,
                                    momentum,
                                    log_posterior,
                                    gradient,
                                    d) {
  current_eps <- 1.0
  pos_mom <- hmc_leapfrog_step(gradient, current_eps, position, momentum, d)
  prop <- function(prop_position, prop_momentum) {
    hmc_log_acceptance(prop_position, prop_momentum, position, momentum, log_posterior)
  }
  position <- pos_mom[seq_len(d)]
  momentum <- pos_mom[-seq_len(d)]
  i <- prop(position, momentum) > log(0.5)
  a <- dplyr::if_else(i, 1.0, -1.0)

  count <- 0
  while (a * prop(position, momentum) > -a * log(2.0) & count < 100) {
    pos_mom <- hmc_leapfrog_step(gradient, current_eps, position, momentum, d)
    position <- pos_mom[seq_len(d)]
    momentum <- pos_mom[-seq_len(d)]
    current_eps <- 2^a * current_eps
    count <- count + 1
  }

  if (count == 100) {
    message("Could not find reasonable epsilon in 100 steps")
  }

  current_eps
}


#' Title
#'
#' @param log_posterior 
#' @param gradient 
#' @param position 
#' @param initial_step_size 
#' @param n_steps 
#' @param iteration 
#' @param warmup_iterations 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
da_step <-
  function(log_posterior,
           gradient,
           position,
           initial_step_size,
           n_steps,
           iteration,
           warmup_iterations,
           d) {
    # Sample initial momentum from standard normal
    momentum <- rnorm(d)
    
    ## Set the step-size
    step_size <- dplyr::if_else(iteration < warmup_iterations, exp(logeps), exp(logepsbar))
    
    ## Perform leapfrog steps
    pos_mom <-
      hmc_leapfrogs(gradient, step_size, n_steps, position, momentum, d)
    prop_position <- pos_mom[seq_len(d)]
    prop_momentum <- pos_mom[-seq_len(d)]
    
    ## Calculate acceptance probability
    a <-
      hmc_log_acceptance(prop_position,
                         prop_momentum,
                         position,
                         momentum,
                         log_posterior)
    
    ## Update the step size in the adaptation phase
    if (iteration < warmup_iterations) {
      update_step_size(iteration, min(1.0, exp(a)), log(10 * initial_step_size))
    }
    
    ## Accept new parameters with probability max(1.0, exp(a))
    if (log(runif(1)) < a) {
      ## return a vector containing a boolean denoting
      ## an accepted move and the new parameters
      c(TRUE, prop_position) 
    } else {
      c(FALSE, position)
    }
  }

#' Hamiltonian Monte Carlo
#'
#' Perform Hamiltonian Monte Carlo with dual averaging
#'
#' @param log_posterior a function for the unnormalised log-posterior
#' from parameters -> log-likelihood
#' @param gradient the gradient of the log_posterior
#' @param n_steps the number of leapfrog steps
#' @param init_parameters a named vector containing the initial parameters
#' @param iters the number of iterations
#' @param warmup_iters the total number of warmup iterations to use for the adaptation phase
#' @param logging enable logging
#'
#' @return a matrix of iterations from the HMC algorithm representing draws
#' from the posterior distribution
#'
#' @examples
hmc_da_helper <-
  function(log_posterior,
           gradient,
           n_steps,
           init_parameters,
           iters,
           warmup_iters = iters / 2, 
           logging = TRUE) {
    
    ## Calculate the dimension of the parameter vector
    d <- length(init_parameters)
    
    ## Sample momentum
    momentum <- rnorm(d)
    
    ## Find a reasonable initial step size
    initial_step_size <-
      find_reasonable_epsilon(init_parameters, momentum, log_posterior, gradient, d)
    
    ## Initialise the (scary) global values for dual averaging
    hm <<- 0
    logepsbar <<- 0
    logeps <<- log(initial_step_size)
    
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
    
    ## MCMC Loop
    for (i in 2:iters) {
      out[i,] <-
        da_step(log_posterior,
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
    
    ## Assign parameter names and create dataframe
    colnames(out) <- c("accepted", names(init_parameters))
    out <- tibble::as_tibble(out) %>%
      dplyr::mutate(iteration = dplyr::row_number())
    
    if (logging) {
      out <- out %>%
        tibble::add_column(hm = hms,
                           log_step_size = logepss,
                           log_step_size_bar = logepsbars)
    } 
    
    out
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
hmc_da <- function(log_posterior, gradient, n_steps, init_parameters, iters, chains = 2) {
  furrr::future_map_dfr(
    .x = seq_len(chains),
    .f = function(x)
      jonnylaw:::hmc_da_helper(
        log_posterior,
        gradient,
        n_steps,
        init_parameters,
        iters
      ),
    .id = "chain"
  )
}
